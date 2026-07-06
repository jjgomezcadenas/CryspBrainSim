# drivers/sigma_r_sweep_dose.jl — the range precision σ_R across a grid of
# doses, i.e. the σ_R-vs-dose curve (validation ladder rung 7, the deliverable
# per scanner). It runs the thinned method of sigma_r_at_dose at each dose,
# reusing one pooled master, and records σ_R for every dose.
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
# σ_R rises as the dose falls, because fewer counts make a noisier edge. The
# nominal-dose point reproduces the sigma_r_at_dose result; the fall-off maps
# how the range precision degrades with dose for this scanner. A different
# scanner has a different photon response and must be simulated separately — its
# curve cannot be produced by thinning this one.
#
# Run:  julia -t auto --project=. drivers/sigma_r_sweep_dose.jl --realizations 100
#       julia -t auto --project=. drivers/sigma_r_sweep_dose.jl --doses 1,0.5,0.2,0.1,0.05
# Writes out/sigma_r/sweep.toml + sweep.npz; the curve comes from
# tools/plot_sigma_r.py --sweep.

include(joinpath(@__DIR__, "sigma_r_common.jl"))

const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 50 : parse(Int, ARGS[i+1])
end
const DOSES = let i = findfirst(==("--doses"), ARGS)
    i === nothing ? [1.0, 0.5, 0.2, 0.1] :
        parse.(Float64, split(ARGS[i+1], ","))
end

function sweep()
    s = setup()
    println("pooling $(length(s.files)) shards…")
    t_pool = @elapsed pool = pool_shards(s.files)
    M_total = length(pool.coinc)
    n_shards = length(s.files)
    tmask = is_true(pool.coinc)
    a_all = attenuation(pool.coinc.xstart, pool.coinc.xend, s.ph)
    @printf("pooled %d coincidences in %.0f s; %d doses × %d realizations\n",
            M_total, t_pool, length(DOSES), REALIZATIONS)

    points = NamedTuple[]
    for dose in DOSES
        target = dose_to_counts(dose, 1.0, M_total, n_shards)
        fits = Float64[]
        crosses = Float64[]
        t = @elapsed for z in 1:REALIZATIONS
            keep = thin_lm(pool.coinc, target, z) .& tmask
            xs, xe = endpoints(pool.coinc, keep)
            res = recon_endpoint(xs, xe, a_all[keep], s.ref, s.base, s.dev)
            push!(fits, res.r50_fit); push!(crosses, res.r50_cross)
        end
        sf = sigma_R(fits; dose_bragg_peak=s.ref.dose_R80)
        sc = sigma_R(crosses; dose_bragg_peak=s.ref.dose_R80)
        push!(points, (dose=dose, target=target, sf=sf, sc=sc))
        @printf("dose %5.3g Gy: keep %9d | σ_R fit %.3f mm | crossing %.3f mm | %d ok/%d fail | %.0f s\n",
                dose, target, sf.sigma, sc.sigma, sf.n_ok, sf.n_fail, t)
    end

    mkpath(OUT)
    npzwrite(joinpath(OUT, "sweep.npz"),
             Dict("dose_Gy" => [p.dose for p in points],
                  "target_counts" => Float64.([p.target for p in points]),
                  "sigma_fit_mm" => [p.sf.sigma for p in points],
                  "sigma_crossing_mm" => [p.sc.sigma for p in points],
                  "mean_fit_mm" => [p.sf.mean for p in points],
                  "n_ok" => Float64.([p.sf.n_ok for p in points]),
                  "n_fail" => Float64.([p.sf.n_fail for p in points])))
    open(joinpath(OUT, "sweep.toml"), "w") do io
        TOML.print(io, Dict(
            "realizations" => REALIZATIONS, "doses_Gy" => DOSES,
            "M_total" => M_total, "n_shards" => n_shards,
            "seed_base" => THINNING_SEED_BASE,
            "point" => [Dict("dose_Gy" => p.dose, "target_counts" => p.target,
                             "sigma_fit_mm" => p.sf.sigma,
                             "sigma_crossing_mm" => p.sc.sigma,
                             "mean_fit_mm" => p.sf.mean, "offset_mm" => p.sf.offset,
                             "n_ok" => p.sf.n_ok, "n_fail" => p.sf.n_fail)
                        for p in points]))
    end
    println("wrote $(joinpath(OUT, "sweep.toml")) (+ sweep.npz)")
end

sweep()
