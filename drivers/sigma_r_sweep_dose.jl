# drivers/sigma_r_sweep_dose.jl — the range precision σ_R across a grid of
# doses, i.e. the σ_R-vs-dose curve (validation ladder rung 7, the deliverable
# per scanner). It reconstructs thinned realizations at each dose, reusing one
# pooled master, and records σ_R for every dose. (Terms — shard, realization,
# dose, R, σ_R — are defined in src/reconstruct.jl.)
#
# Options:
#
#   --realizations N     (default 50) realizations drawn at each dose. More
#                        realizations pin each σ_R tighter (spread known to
#                        ±1/√(2(N−1))).
#
#   --doses D1,D2,...    (default 1,0.5,0.2,0.1) the dose grid in Gy, highest
#                        first. Activity scales with dose, so a lower dose keeps
#                        a smaller fraction of the pooled coincidences.
#
#   --all-events         Reconstruct every kept coincidence — scatters and
#                        randoms in, uncorrected (the working protocol). The
#                        default is the trues-only selection. The thin draw is
#                        seeded identically either way, so the two selections
#                        pair realization by realization. Outputs carry the
#                        `_all` suffix (sweep_all.toml/npz).
#
# σ_R rises as the dose falls, because fewer counts make a noisier edge. The
# nominal-dose point reproduces the sigma_r_at_dose result; the fall-off maps
# how the range precision degrades with dose for this scanner. A different
# scanner has a different photon response and must be simulated separately — its
# curve cannot be produced by thinning this one.
#
# Run:  julia -t auto --project=. drivers/sigma_r_sweep_dose.jl --realizations 100
#       julia -t auto --project=. drivers/sigma_r_sweep_dose.jl --doses 1,0.5,0.2,0.1,0.05
# Writes sweep.toml + sweep.npz into the config's sigma_r/; the curve comes from
# tools/plot_sigma_r.py --sweep.

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf
using TOML

const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 50 : parse(Int, ARGS[i+1])
end
const DOSES = let i = findfirst(==("--doses"), ARGS)
    i === nothing ? [1.0, 0.5, 0.2, 0.1] : parse.(Float64, split(ARGS[i+1], ","))
end
const ALL_EVENTS = "--all-events" in ARGS


const DEV = Metal.functional() ? MtlArray : identity

function context()
    params = load_run_parameters()
    cfg = params.config        # the active arm (run_parameters.toml)
    cache = joinpath(sensitivity_out(cfg.scenario, cfg.topology, cfg.scanner),
                     sensitivity_cache_name(params))
    ctx = load_run_context(;
        products_root=joinpath(dirname(dirname(@__DIR__)), "PtCryspProds"),
        scenario=cfg.scenario, topology=cfg.topology, scanner=cfg.scanner,
        crystal=cfg.crystal, leaf=cfg.leaf, sens_cache=cache, params=params)
    write_descriptors(ctx)          # stamp geometry.toml + crystal.toml
    return ctx
end

function sweep(ctx)
    out = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                   "sigma_r")
    println("pooling $(length(ctx.files)) shards…")
    t_pool = @elapsed pool = pool_shards(ctx.files)
    M_total = length(pool.coinc)
    n_shards = length(ctx.files)
    selection = ALL_EVENTS ? "all-uncorrected" : "trues-only"
    tmask = ALL_EVENTS ? trues(M_total) : is_true(pool.coinc)
    a_all = lor_attenuation(ctx, pool.coinc.xstart, pool.coinc.xend)
    @printf("pooled %d coincidences in %.0f s; %d doses × %d realizations; %s\n",
            M_total, t_pool, length(DOSES), REALIZATIONS, selection)

    points = NamedTuple[]
    for dose in DOSES
        target = dose_to_counts(dose, 1.0, M_total, n_shards)
        fits = Float64[]
        crosses = Float64[]
        t = @elapsed for z in 1:REALIZATIONS
            keep = thin_lm(pool.coinc, target, z) .& tmask
            xs, xe = endpoints(pool.coinc, keep)
            res = reconstruct_endpoint(ctx, xs, xe, a_all[keep]; device=DEV)
            push!(fits, res.r50_fit); push!(crosses, res.r50_cross)
        end
        sf = sigma_R(fits; dose_bragg_peak=ctx.ref.dose_R80)
        sc = sigma_R(crosses; dose_bragg_peak=ctx.ref.dose_R80)
        push!(points, (dose=dose, target=target, sf=sf, sc=sc))
        @printf("dose %5.3g Gy: keep %9d | σ_R fit %.3f mm | crossing %.3f mm | %d ok/%d fail | %.0f s\n",
                dose, target, sf.sigma, sc.sigma, sf.n_ok, sf.n_fail, t)
    end

    stem = ALL_EVENTS ? "sweep_all" : "sweep"
    mkpath(out)
    npzwrite(joinpath(out, "$(stem).npz"),
             Dict("dose_Gy" => [p.dose for p in points],
                  "target_counts" => Float64.([p.target for p in points]),
                  "sigma_fit_mm" => [p.sf.sigma for p in points],
                  "sigma_crossing_mm" => [p.sc.sigma for p in points],
                  "mean_fit_mm" => [p.sf.mean for p in points],
                  "n_ok" => Float64.([p.sf.n_ok for p in points]),
                  "n_fail" => Float64.([p.sf.n_fail for p in points])))
    open(joinpath(out, "$(stem).toml"), "w") do io
        TOML.print(io, Dict(
            "realizations" => REALIZATIONS, "doses_Gy" => DOSES,
            "M_total" => M_total, "n_shards" => n_shards,
            "selection" => selection, "profile" => "whole-plane",
            "seed_base" => THINNING_SEED_BASE,
            "point" => [Dict("dose_Gy" => p.dose, "target_counts" => p.target,
                             "sigma_fit_mm" => p.sf.sigma,
                             "sigma_crossing_mm" => p.sc.sigma,
                             "mean_fit_mm" => p.sf.mean, "offset_mm" => p.sf.offset,
                             "n_ok" => p.sf.n_ok, "n_fail" => p.sf.n_fail)
                        for p in points]))
    end
    println("wrote $(joinpath(out, "$(stem).toml")) (+ $(stem).npz)")
end

sweep(context())
