# drivers/statistical_procedure.jl — reference-case validation used by the
# paper's Statistical procedure section and thinning appendix.
#
# The default case is the reference TBP BGO scanner, 1 Gy, with the 300 s
# acquisition beginning 120 s after irradiation.  It compares the ten truly
# independent 1-Gy shards with N repeated thinnings of their fixed 10-Gy pool.
#
# Run:
#   julia -t auto --project=. drivers/statistical_procedure.jl
#       [--config config/run_parameters_ring_bgo_v2.toml]
#       [--leaf del120s_ac300s_1Gy] [--realizations 200]
#       [--dose 1.0] [--seed 9000000]

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using Statistics: mean, median, std
using TOML

argvalue(flag, default, parsefn=identity) = let index = findfirst(==(flag), ARGS)
    index === nothing ? default : parsefn(ARGS[index + 1])
end

const REPO = dirname(@__DIR__)
const PRODUCTS_ROOT = joinpath(dirname(REPO), "PtCryspProds")
const CONFIG_PATH = abspath(argvalue(
    "--config",
    joinpath(REPO, "config", "run_parameters_ring_bgo_v2.toml"),
    String,
))
const PARAMS = load_run_parameters(CONFIG_PATH)
const CFG = PARAMS.config
const LEAF = argvalue("--leaf", "del120s_ac300s_1Gy", String)
const REALIZATIONS = argvalue("--realizations", 200, value -> parse(Int, value))
const DOSE_GY = argvalue("--dose", 1.0, value -> parse(Float64, value))
const SEED_BASE = argvalue("--seed", 9_000_000, value -> parse(Int, value))
const WASHED_ONLY = let index = findfirst(==("--washed-only"), ARGS)
    index !== nothing
end
const MAX_FAILURE_FRACTION = 0.05
const DEV = Metal.functional() ? MtlArray : identity

function validate_arguments()
    isfile(CONFIG_PATH) || error("configuration does not exist: $CONFIG_PATH")
    REALIZATIONS >= 20 || error("--realizations must be at least 20")
    DOSE_GY > 0 || error("--dose must be positive")
end

function fit_summary(values)
    valid = filter(isfinite, values)
    length(valid) >= 2 || error("fewer than two valid endpoint fits")
    failures = length(values) - length(valid)
    allowed = floor(Int, MAX_FAILURE_FRACTION * length(values))
    failures <= allowed ||
        error("$failures/$(length(values)) fits failed; at most $allowed are allowed")
    return (mean=mean(valid), sigma=std(valid), failures=failures)
end

function reconstruct_selection(ctx, coinc, attenuation)
    xs, xe = endpoints(coinc)
    return reconstruct_endpoint(ctx, xs, xe, attenuation; device=DEV)
end

function main()
    validate_arguments()
    cache = joinpath(
        sensitivity_out(CFG.scenario, CFG.topology, CFG.scanner),
        sensitivity_cache_name(PARAMS),
    )
    ctx = load_run_context(
        products_root=PRODUCTS_ROOT,
        scenario=CFG.scenario,
        topology=CFG.topology,
        scanner=CFG.scanner,
        crystal=CFG.crystal,
        leaf=LEAF,
        sens_cache=cache,
        params=PARAMS,
    )
    length(ctx.files) == 10 || error("reference validation requires ten shards")
    attrs = shard_attrs(ctx.files[1])
    shard_generation(attrs) == "v2" || error("reference validation requires v2 shards")

    # Direct estimate: reconstruct each statistically independent 1-Gy shard.
    shard_r50 = Float64[]
    shard_error = Float64[]
    shard_chi2 = Float64[]
    shard_counts = Int[]
    shard_seconds = @elapsed for (index, file) in enumerate(ctx.files)
        shard = read_shard(file)
        shard.n_dropped == 0 || error("shard $(index - 1) dropped LORs")
        xs, xe = endpoints(shard.coinc)
        attenuation = lor_attenuation(ctx, xs, xe)
        result = reconstruct_endpoint(ctx, xs, xe, attenuation; device=DEV)
        push!(shard_r50, result.r50_fit)
        push!(shard_error, result.z0_err)
        push!(shard_chi2, result.erfc_chi2_dof)
        push!(shard_counts, result.nev)
        @printf("shard %02d: %d events, R50 %.4f ± %.4f mm, chi2/dof %.3f\n",
                index - 1, result.nev, result.r50_fit, result.z0_err,
                result.erfc_chi2_dof)
    end
    direct = fit_summary(shard_r50)

    # Conditional estimate: independent Bernoulli draws from the same fixed pool.
    pool = pool_shards(ctx.files)
    pool.n_dropped == 0 || error("pooled sample dropped $(pool.n_dropped) LORs")
    isotope = vcat([shard_isotope(file) for file in ctx.files]...)
    length(isotope) == length(pool.coinc) || error("isotope/pool length mismatch")
    survival = Float64.(attrs["washout_g"])
    all(id -> 0 <= Int(id) < length(survival), isotope) ||
        error("isotope identifier outside washout_g table")
    event_survival = [survival[Int(id) + 1] for id in isotope]
    p_dose = DOSE_GY / length(ctx.files)
    0 < p_dose <= 0.1 || error("dose keep probability must be in (0,0.1]")
    q_nominal = p_dose
    q_washed = p_dose .* event_survival
    correction_nominal = finite_pool_correction(q_nominal)
    correction_washed = finite_pool_correction(q_washed)
    xs_all, xe_all = endpoints(pool.coinc)
    attenuation_all = lor_attenuation(ctx, xs_all, xe_all)

    nominal_r50 = Float64[]; washed_r50 = Float64[]
    nominal_error = Float64[]; washed_error = Float64[]
    nominal_chi2 = Float64[]; washed_chi2 = Float64[]
    nominal_counts = Int[]; washed_counts = Int[]
    thinning_seconds = @elapsed for realization in 1:REALIZATIONS
        rng = MersenneTwister(SEED_BASE + realization)
        uniform = rand(rng, length(pool.coinc))
        nominal_keep = uniform .< q_nominal
        washed_keep = uniform .< q_washed
        all(washed_keep .<= nominal_keep) || error("paired subset invariant failed")

        nominal = WASHED_ONLY ? nothing : reconstruct_endpoint(
            ctx, xs_all[:, nominal_keep], xe_all[:, nominal_keep],
            attenuation_all[nominal_keep]; device=DEV)
        washed = reconstruct_endpoint(
            ctx,
            xs_all[:, washed_keep],
            xe_all[:, washed_keep],
            attenuation_all[washed_keep];
            device=DEV,
        )
        push!(nominal_r50, WASHED_ONLY ? NaN : nominal.r50_fit)
        push!(washed_r50, washed.r50_fit)
        push!(nominal_error, WASHED_ONLY ? NaN : nominal.z0_err)
        push!(washed_error, washed.z0_err)
        push!(nominal_chi2, WASHED_ONLY ? NaN : nominal.erfc_chi2_dof)
        push!(washed_chi2, washed.erfc_chi2_dof)
        push!(nominal_counts, WASHED_ONLY ? 0 : nominal.nev)
        push!(washed_counts, washed.nev)
        realization % 10 == 0 && @printf("thinning %d/%d\n", realization, REALIZATIONS)
    end
    nominal = WASHED_ONLY ?
        (mean=NaN, sigma=NaN, failures=REALIZATIONS) : fit_summary(nominal_r50)
    washed = fit_summary(washed_r50)

    out = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                   "statistical_procedure")
    mkpath(out)
    output = joinpath(out, "reference_bgo_t120.toml")
    open(output, "w") do io
        TOML.print(io, Dict(
            "scanner" => ctx.ring,
            "crystal" => CFG.crystal,
            "leaf" => LEAF,
            "dose_Gy" => DOSE_GY,
            "realizations" => REALIZATIONS,
            "seed_base" => SEED_BASE,
            "washed_only" => WASHED_ONLY,
            "fit_failure_fraction_limit" => MAX_FAILURE_FRACTION,
            "finite_N_relative_sigma_uncertainty" =>
                1 / sqrt(2 * (REALIZATIONS - 1)),
            "independent_shards" => Dict(
                "count" => length(shard_r50),
                "R50_mm" => shard_r50,
                "R50_fit_error_mm" => shard_error,
                "erfc_chi2_dof" => shard_chi2,
                "events" => shard_counts,
                "mean_R50_mm" => direct.mean,
                "sigma_R_mm" => direct.sigma,
                "n_fail" => direct.failures,
                "relative_sigma_uncertainty" => 1 / sqrt(2 * (length(shard_r50) - 1)),
            ),
            "thinning" => Dict(
                "dose_probability" => p_dose,
                "nominal_R50_mm" => nominal_r50,
                "washed_R50_mm" => washed_r50,
                "nominal_R50_fit_error_mm" => nominal_error,
                "washed_R50_fit_error_mm" => washed_error,
                "nominal_erfc_chi2_dof" => nominal_chi2,
                "washed_erfc_chi2_dof" => washed_chi2,
                "nominal_events" => nominal_counts,
                "washed_events" => washed_counts,
                "nominal_mean_R50_mm" => nominal.mean,
                "washed_mean_R50_mm" => washed.mean,
                "nominal_sigma_R_raw_mm" => nominal.sigma,
                "washed_sigma_R_raw_mm" => washed.sigma,
                "nominal_finite_pool_correction" => correction_nominal,
                "washed_finite_pool_correction" => correction_washed,
                "nominal_sigma_R_corrected_mm" => correction_nominal * nominal.sigma,
                "washed_sigma_R_corrected_mm" => correction_washed * washed.sigma,
                "n_fail_nominal" => nominal.failures,
                "n_fail_washed" => washed.failures,
            ),
            "timing_s" => Dict(
                "independent_shards" => shard_seconds,
                "thinning" => thinning_seconds,
            ),
        ))
    end
    @printf("direct sigma_R %.4f mm; nominal thinning raw/corrected %.4f/%.4f mm\n",
            direct.sigma, nominal.sigma, correction_nominal * nominal.sigma)
    @printf("washed thinning raw/corrected %.4f/%.4f mm\n",
            washed.sigma, correction_washed * washed.sigma)
    println("wrote $output")
end

main()
