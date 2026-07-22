# drivers/sigma_r_niter_v2.jl — MLEM-iteration stability of the washed R50
# precision for isotope-labelled generation-2 samples.
#
# Each v2 acquisition leaf already represents its complete [t1,t2] window, and
# the stored isotope label selects the corresponding survival factor. One MLEM
# pass per realisation records R50 at every iteration checkpoint.
#
# Run:
#   julia -t auto --project=. drivers/sigma_r_niter_v2.jl
#       [--leaves del120s_ac300s_1Gy,del180s_ac300s_1Gy]
#       [--dose 1.0] [--realizations 50]
#       [--niter-max 150] [--check-every 10] [--seed 8000000]
#
# Writes <config>/washout_v2/sigma_r_vs_niter_v2.toml

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using TOML

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const PARAMS = load_run_parameters()
const CFG = PARAMS.config
const SCENARIO, TOPOLOGY, RING = CFG.scenario, CFG.topology, CFG.scanner
const DEV = Metal.functional() ? MtlArray : identity

argvalue(flag, default, parsefn) = let i = findfirst(==(flag), ARGS)
    i === nothing ? default : parsefn(ARGS[i + 1])
end

const LEAVES = argvalue(
    "--leaves",
    ["del120s_ac300s_1Gy", "del180s_ac300s_1Gy", "del300s_ac300s_1Gy"],
    value -> String.(split(value, ",")),
)
const THIN_DOSE = argvalue("--dose", 1.0, value -> parse(Float64, value))
const REALIZATIONS = argvalue("--realizations", 50, value -> parse(Int, value))
const NITER_MAX = argvalue("--niter-max", 150, value -> parse(Int, value))
const CHECK_EVERY = argvalue("--check-every", 10, value -> parse(Int, value))
const SEED_BASE = argvalue("--seed", 8_000_000, value -> parse(Int, value))

function validate_arguments()
    REALIZATIONS >= 2 || error("--realizations must be at least 2")
    0 < THIN_DOSE <= 1 || error("--dose must be in (0,1] Gy")
    CHECK_EVERY > 0 || error("--check-every must be positive")
    NITER_MAX >= CHECK_EVERY ||
        error("--niter-max must be at least --check-every")
    NITER_MAX % CHECK_EVERY == 0 ||
        error("--niter-max must be divisible by --check-every")
    isempty(LEAVES) && error("--leaves must contain at least one v2 leaf")
end

function r50_vs_iteration(ctx, xs, xe, attenuation; device=DEV)
    params = ctx.params
    sensitivity = scaled_sensitivity(ctx.base, size(xs, 2), params.n_sens)
    model = ListmodePoissonModel(
        device(xs),
        device(xe),
        device(sensitivity);
        img_origin=params.grid.img_origin,
        voxsize=params.grid.voxsize,
        mult=device(attenuation),
    )

    image = device(Float32.(sensitivity .> 0))
    iterations = collect(CHECK_EVERY:CHECK_EVERY:NITER_MAX)
    fitted_r50 = fill(NaN, length(iterations))
    completed = 0

    for (index, target_iteration) in enumerate(iterations)
        image = mlem(model, image; niter=target_iteration - completed)
        completed = target_iteration
        z, profile = depth_profile(
            Array(image);
            voxel_size_mm=params.grid.voxsize,
            beam_axis=3,
            roi_radius_mm=params.roi.radius_mm,
            z_origin_mm=params.grid.img_origin[3],
        )
        fit = fit_endpoint(z, profile; window=ctx.ref.window, weighted=true)
        fitted_r50[index] = fit.z0
    end

    return iterations, fitted_r50
end

function run_leaf(leaf)
    cache = joinpath(
        sensitivity_out(SCENARIO, TOPOLOGY, RING),
        sensitivity_cache_name(PARAMS),
    )
    ctx = load_run_context(
        products_root=ROOT,
        scenario=SCENARIO,
        topology=TOPOLOGY,
        scanner=RING,
        crystal=CFG.crystal,
        leaf=leaf,
        sens_cache=cache,
        params=PARAMS,
    )

    attrs = shard_attrs(ctx.files[1])
    shard_generation(attrs) == "v2" ||
        error("$leaf: sigma_r_niter_v2 requires generation-2 shards")

    pool = pool_shards(ctx.files)
    pool.n_dropped == 0 || error("$leaf: pool dropped $(pool.n_dropped) LORs")
    isotope = vcat([shard_isotope(file) for file in ctx.files]...)
    length(isotope) == length(pool.coinc) ||
        error("$leaf: isotope labels are not aligned with the pooled LORs")

    survival = Float64.(attrs["washout_g"])
    all(id -> 0 <= Int(id) < length(survival), isotope) ||
        error("$leaf: isotope identifier outside the washout_g table")
    event_survival = [survival[Int(id) + 1] for id in isotope]

    n_shards = length(ctx.files)
    p_dose = THIN_DOSE / n_shards
    0 < p_dose <= 0.1 ||
        error("$leaf: keep probability must be in (0,0.1]; got $p_dose")
    event_probability = p_dose .* event_survival
    correction = finite_pool_correction(event_probability)

    event_count = length(pool.coinc)
    attenuation = lor_attenuation(ctx, pool.coinc.xstart, pool.coinc.xend)
    iterations = collect(CHECK_EVERY:CHECK_EVERY:NITER_MAX)
    endpoints_by_iteration = fill(NaN, REALIZATIONS, length(iterations))
    delay = Float64(attrs["t_del_s"])

    @printf(
        "%s: %d pooled LORs, %d shards, p_dose %.3f, N=%d\n",
        leaf,
        event_count,
        n_shards,
        p_dose,
        REALIZATIONS,
    )

    elapsed = @elapsed for realization in 1:REALIZATIONS
        rng = MersenneTwister(SEED_BASE + Int(round(delay)) * 1000 + realization)
        keep = rand(rng, event_count) .< event_probability
        xs, xe = endpoints(pool.coinc, keep)
        returned_iterations, fitted_r50 = r50_vs_iteration(
            ctx,
            xs,
            xe,
            attenuation[keep];
            device=DEV,
        )
        returned_iterations == iterations || error("iteration checkpoint mismatch")
        endpoints_by_iteration[realization, :] = fitted_r50
    end

    summaries = [sigma_R(view(endpoints_by_iteration, :, index))
                 for index in eachindex(iterations)]
    sigma_raw = [summary.sigma for summary in summaries]
    sigma = correction .* sigma_raw
    mean_r50 = [summary.mean for summary in summaries]
    failures = [summary.n_fail for summary in summaries]
    best = argmin(map(value -> isfinite(value) ? value : Inf, sigma))

    @printf(
        "  minimum sigma_R %.3f mm at %d iterations; %d s\n",
        sigma[best],
        iterations[best],
        round(Int, elapsed),
    )

    return (
        context=ctx,
        leaf=leaf,
        delay=delay,
        iterations=iterations,
        sigma_raw=sigma_raw,
        sigma=sigma,
        correction=correction,
        mean_r50=mean_r50,
        failures=failures,
    )
end

function main()
    validate_arguments()
    results = [run_leaf(leaf) for leaf in LEAVES]
    ctx = results[1].context
    out = joinpath(
        config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
        "washout_v2",
    )
    mkpath(out)
    output = joinpath(out, "sigma_r_vs_niter_v2.toml")
    open(output, "w") do io
        TOML.print(io, Dict(
            "generation" => "v2",
            "method" => "isotope-indexed window-integrated survival on labelled v2 LORs",
            "scanner" => RING,
            "crystal" => CFG.crystal,
            "dose_Gy" => THIN_DOSE,
            "realizations" => REALIZATIONS,
            "frozen_niter" => PARAMS.niter,
            "seed_base" => SEED_BASE,
            "niter" => results[1].iterations,
            "point" => [Dict(
                "leaf" => result.leaf,
                "t_del_s" => result.delay,
                "sigma_R_convention" => "finite-pool corrected",
                "finite_pool_correction" => result.correction,
                "washed_sigma_R_raw_mm" => result.sigma_raw,
                "washed_sigma_R_mm" => result.sigma,
                "washed_mean_R50_mm" => result.mean_r50,
                "n_fail" => result.failures,
            ) for result in results],
        ))
    end
    println("wrote $output")
end

main()
