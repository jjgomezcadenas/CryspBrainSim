# drivers/sigma_r_at_dose.jl — the range precision σ_R at a single dose
# (validation ladder rung 6). σ_R is the standard deviation of the fitted range
# endpoint R across many acquisitions; this driver measures it at one dose, two
# ways. Run both at the nominal dose and they must agree — that agreement is the
# rung-6 gate that certifies the thinned method.
#
# Options:
#
#   --from-shards        The reference method. Reconstruct each of the ten
#                        shards on its own and take the spread of the ten range
#                        values. A shard is one full, independent simulated
#                        acquisition at the nominal dose, so this σ_R carries no
#                        thinning assumption. Ten values pin σ_R to about ±24%.
#
#   --realizations N     (default 50) The production method. Pool the ten shards
#                        into one master, then draw N realizations by thinning —
#                        randomly keeping a fraction of the pooled coincidences,
#                        each realization a fresh emulated acquisition. σ_R is
#                        the spread of the N range values. Fifty pins σ_R to
#                        about ±10%.
#
#   --dose D             (default 1.0 = the nominal/top dose) The emulated dose
#                        in Gy for the thinned method. Activity scales with dose,
#                        so a lower D keeps a smaller fraction of the pool —
#                        fewer counts, a noisier edge, a larger σ_R. The
#                        --from-shards method is always at the nominal dose.
#
# At the nominal dose the thinned run is compared against the --from-shards
# reference (run it first); the two σ_R values agreeing is the gate.
#
# Run:  julia -t auto --project=. drivers/sigma_r_at_dose.jl --from-shards
#       julia -t auto --project=. drivers/sigma_r_at_dose.jl --realizations 50
#       julia -t auto --project=. drivers/sigma_r_at_dose.jl --realizations 50 --dose 0.5
# Writes into out/sigma_r/; figures come from tools/plot_sigma_r.py.

include(joinpath(@__DIR__, "sigma_r_common.jl"))

const FROM_SHARDS = "--from-shards" in ARGS
const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 50 : parse(Int, ARGS[i+1])
end
const DOSE_GY = let i = findfirst(==("--dose"), ARGS)
    i === nothing ? 1.0 : parse(Float64, ARGS[i+1])
end

# The reference: σ_R from the ten shards reconstructed independently, no thinning.
function from_shards(s)
    println("$(length(s.files)) shards; niter $(PARAMS.niter), " *
            "ROI $(PARAMS.roi.radius_mm) mm, window $(round.(s.ref.window; digits=2))")
    results = NamedTuple[]
    t_total = @elapsed for f in s.files
        r = read_shard(f)
        xs, xe = endpoints(r.coinc, is_true(r.coinc))
        t = @elapsed res = recon_endpoint(xs, xe, attenuation(xs, xe, s.ph),
                                          s.ref, s.base, s.dev)
        push!(results, (shard=Int(r.attrs["realization"]), res...))
        @printf("shard %d: R50 fit %8.3f ± %.3f mm | crossing %8.3f mm | %d ev | %.1f s\n",
                results[end].shard, res.r50_fit, res.z0_err, res.r50_cross, res.nev, t)
    end

    fits = [r.r50_fit for r in results]
    errs = [r.z0_err for r in results]
    sf = sigma_R(fits; dose_bragg_peak=s.ref.dose_R80)
    sc = sigma_R([r.r50_cross for r in results]; dose_bragg_peak=s.ref.dose_R80)
    report_sigma("from the 10 shards at 1 Gy", sf, sc)
    @printf("  each fit's own error averages %.3f mm, %.1f× the measured spread\n",
            mean(errs), mean(errs) / sf.sigma)
    @printf("  wall-clock %.0f s\n", t_total)

    mkpath(OUT)
    npzwrite(joinpath(OUT, "from_shards.npz"),
             Dict("shards" => Float64.([r.shard for r in results]),
                  "r50_fit" => fits, "z0_err" => errs,
                  "r50_cross" => [r.r50_cross for r in results],
                  "w" => [r.w for r in results],
                  "nev" => Float64.([r.nev for r in results])))
    open(joinpath(OUT, "from_shards.toml"), "w") do io
        TOML.print(io, Dict(
            "method" => "ten independent shards, no thinning",
            "n_shards" => length(results), "dose_Gy" => 1.0,
            "niter" => PARAMS.niter, "window_mm" => collect(s.ref.window),
            "sens" => Dict("cache" => SENS_CACHE,
                           "recocrysp_sha" => s.meta["recocrysp_sha"]),
            "per_shard" => Dict(
                "shard" => [r.shard for r in results],
                "n_events" => [r.nev for r in results],
                "r50_fit_mm" => fits, "z0_err_mm" => errs,
                "r50_crossing_mm" => [r.r50_cross for r in results],
                "w_mm" => [r.w for r in results]),
            "sigma_R_fit" => Dict("mean_mm" => sf.mean, "sigma_mm" => sf.sigma,
                                  "sem_mm" => sf.sem, "offset_mm" => sf.offset),
            "sigma_R_crossing" => Dict("mean_mm" => sc.mean, "sigma_mm" => sc.sigma,
                                       "sem_mm" => sc.sem, "offset_mm" => sc.offset),
            "fit_error_over_spread" => mean(errs) / sf.sigma,
            "reference" => Dict("dose_R80_mm" => s.ref.dose_R80,
                                "activity_R50_fit_mm" => s.ref.activity_R50_fit,
                                "activity_R50_crossing_mm" => s.ref.activity_R50),
            "timing_s" => t_total))
    end
    println("wrote $(joinpath(OUT, "from_shards.toml")) (+ from_shards.npz)")
end

# The production method: σ_R from N thinned realizations at the chosen dose.
function thinned(s)
    println("pooling $(length(s.files)) shards…")
    t_pool = @elapsed pool = pool_shards(s.files)
    M_total = length(pool.coinc)
    n_shards = length(s.files)
    tmask = is_true(pool.coinc)
    a_all = attenuation(pool.coinc.xstart, pool.coinc.xend, s.ph)
    target = dose_to_counts(DOSE_GY, 1.0, M_total, n_shards)
    @printf("pooled %d coincidences in %.0f s; dose %.3g Gy keeps %d per realization (fraction %.4f)\n",
            M_total, t_pool, DOSE_GY, target, target / M_total)

    results = NamedTuple[]
    t_total = @elapsed for z in 1:REALIZATIONS
        keep = thin_lm(pool.coinc, target, z) .& tmask
        xs, xe = endpoints(pool.coinc, keep)
        res = recon_endpoint(xs, xe, a_all[keep], s.ref, s.base, s.dev)
        push!(results, res)
        @printf("realization %2d: R50 fit %8.3f ± %.3f mm | crossing %8.3f mm | %d ev\n",
                z, res.r50_fit, res.z0_err, res.r50_cross, res.nev)
    end

    fits = [r.r50_fit for r in results]
    sf = sigma_R(fits; dose_bragg_peak=s.ref.dose_R80)
    sc = sigma_R([r.r50_cross for r in results]; dose_bragg_peak=s.ref.dose_R80)
    report_sigma("from $REALIZATIONS thinned realizations at $(DOSE_GY) Gy", sf, sc)

    # At the nominal dose, compare against the shard reference — the rung-6 gate.
    gate = Dict{String,Any}()
    ref_toml = joinpath(OUT, "from_shards.toml")
    if DOSE_GY == 1.0 && isfile(ref_toml)
        rt = TOML.parsefile(ref_toml)
        σ_ref = rt["sigma_R_fit"]["sigma_mm"]
        n_ref = rt["n_shards"]
        ratio = sf.sigma / σ_ref
        band = sqrt(1 / (2 * (sf.n_ok - 1)) + 1 / (2 * (n_ref - 1)))
        pass = abs(ratio - 1) < 2 * band
        @printf("  gate vs the shard reference %.3f mm (n = %d): ratio %.2f, within ±%.0f%% → %s\n",
                σ_ref, n_ref, ratio, 200band, pass ? "AGREES" : "check")
        gate = Dict("reference_sigma_mm" => σ_ref, "ratio" => ratio,
                    "band_2sigma" => 2band, "pass" => pass)
    end
    @printf("  wall-clock %.0f s\n", t_total)

    tag = dose_tag(DOSE_GY)
    mkpath(OUT)
    npzwrite(joinpath(OUT, "at_dose_$(tag).npz"),
             Dict("r50_fit" => fits, "z0_err" => [r.z0_err for r in results],
                  "r50_cross" => [r.r50_cross for r in results],
                  "nev" => Float64.([r.nev for r in results])))
    open(joinpath(OUT, "at_dose_$(tag).toml"), "w") do io
        TOML.print(io, Dict(
            "method" => "thinned realizations of the pooled master",
            "realizations" => REALIZATIONS, "dose_Gy" => DOSE_GY,
            "target_counts" => target, "M_total" => M_total, "n_shards" => n_shards,
            "seed_base" => THINNING_SEED_BASE,
            "sigma_R_fit" => Dict("mean_mm" => sf.mean, "sigma_mm" => sf.sigma,
                                  "sem_mm" => sf.sem, "offset_mm" => sf.offset),
            "sigma_R_crossing" => Dict("mean_mm" => sc.mean, "sigma_mm" => sc.sigma),
            "gate" => gate, "timing_s" => t_total))
    end
    println("wrote $(joinpath(OUT, "at_dose_$(tag).toml")) (+ at_dose_$(tag).npz)")
end

let s = setup()
    FROM_SHARDS ? from_shards(s) : thinned(s)
end
