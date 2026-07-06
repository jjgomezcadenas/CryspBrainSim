# drivers/shard_crosscheck.jl — ladder rung 6, first half: the ten stored
# shards are ten INDEPENDENT full-statistics acquisitions at the top dose
# (matched master_seed, distinct realization indices), so fitting each one
# independently and taking the std gives a bias-free σ_R at 1 Gy with no
# thinning involved. The std of n = 10 samples carries a relative error of
# 1/√(2(n−1)) ≈ 24% — this anchors σ_R's scale; the thinned σ_R at top dose
# (rung 6, second half, once thinning lands) must agree within that band.
#
# Each shard runs the frozen single-shard chain (config/run_parameters.toml):
# trues-only → ellipsoid attenuation → MLEM at the frozen iteration count →
# disc-ROI profile → windowed erfc fit + windowed crossing. Both conventions
# are collapsed by sigma_R; the ratio σ_R / mean(z0_err) calibrates the
# per-fit counting-statistics error against the true ensemble spread.
#
# Run:  julia -t auto --project=. drivers/shard_crosscheck.jl
# Writes out/shard_crosscheck/results.toml + endpoints.npz; the figure comes
# from tools/plot_crosscheck.py.

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf
using Statistics: mean
using TOML

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const SCEN = joinpath(ROOT, "uniform_headep_sobp_1e8")

const PARAMS = load_run_parameters()
const N = PARAMS.grid.n
const VS = PARAMS.grid.voxsize
const ORG = PARAMS.grid.img_origin
const SENS_CACHE = joinpath(dirname(@__DIR__), "out", "sensitivity",
    "crysp_ring_1m_grid64x64x96_1.5mm_orgm47.25_m47.25_m119.25_n$(PARAMS.n_sens)")

function reconstruct_shard(file, ref, base, ph, dev)
    r = read_shard(file)
    tmask = is_true(r.coinc)
    xs, xe = endpoints(r.coinc, tmask)
    nev = size(xs, 2)
    mult = attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes,
                                 centre=ph.centre, mu_mm_inv=ph.mu_mm_inv)
    sens = scaled_sensitivity(base, nev, PARAMS.n_sens)
    model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                 img_origin=ORG, voxsize=VS, mult=dev(mult))
    x = mlem(model, dev(Float32.(sens .> 0)); niter=PARAMS.niter)
    z, prof = depth_profile(Array(x); voxel_size_mm=VS, beam_axis=3,
                            roi_radius_mm=PARAMS.roi.radius_mm, z_origin_mm=ORG[3])
    fit = fit_endpoint(z, prof; window=ref.window, weighted=true)
    cross = windowed_crossing(z, prof, ref.window)
    return (shard=Int(r.attrs["realization"]), nev=nev,
            r50_fit=fit.z0, z0_err=fit.z0_err, w=fit.w, r50_cross=cross)
end

function main()
    ref = characterize(SCEN)
    base, meta = load_sensitivity(SENS_CACHE)
    g = meta["grid"]
    (g["n"] == collect(N) && Float32.(g["img_origin"]) == collect(ORG) &&
     Float32.(g["voxsize"]) == collect(VS)) ||
        error("shard_crosscheck: sensitivity cache grid ≠ run-parameter grid")
    ph = phantom_attenuation(SCEN)
    dev = Metal.functional() ? MtlArray : identity
    leaf = leaf_dir(ROOT; scenario="uniform_headep_sobp_1e8",
                    scanner="crysp_ring_1m", crystal="bgo", leaf="fast_1Gy")
    files = shard_files(leaf)
    println("$(length(files)) shards; frozen niter $(PARAMS.niter), " *
            "ROI $(PARAMS.roi.radius_mm) mm, window $(ref.window)")

    results = NamedTuple[]
    t_total = @elapsed for f in files
        t = @elapsed res = reconstruct_shard(f, ref, base, ph, dev)
        push!(results, res)
        @printf("shard %d: R50 fit %8.3f ± %.3f mm | crossing %8.3f mm | w %.2f mm | %d ev | %.1f s\n",
                res.shard, res.r50_fit, res.z0_err, res.r50_cross, res.w, res.nev, t)
    end

    fits = [r.r50_fit for r in results]
    crosses = [r.r50_cross for r in results]
    errs = [r.z0_err for r in results]
    sf = sigma_R(fits; dose_bragg_peak=ref.dose_R80)
    sc = sigma_R(crosses; dose_bragg_peak=ref.dose_R80)
    calib = sf.sigma / mean(errs)

    @printf("\nsigma_R at 1 Gy (n = %d shards, std rel. error ≈ %.0f%%):\n",
            sf.n_ok, 100 / sqrt(2 * (sf.n_ok - 1)))
    @printf("  fit convention:      mean %8.3f mm | σ_R %.3f mm | sem %.3f mm | offset to dose-R80 %.3f mm\n",
            sf.mean, sf.sigma, sf.sem, sf.offset)
    @printf("  crossing convention: mean %8.3f mm | σ_R %.3f mm | sem %.3f mm | offset to dose-R80 %.3f mm\n",
            sc.mean, sc.sigma, sc.sem, sc.offset)
    @printf("  z0_err calibration:  σ_R / mean(z0_err) = %.3f / %.3f = %.2f\n",
            sf.sigma, mean(errs), calib)
    @printf("  total wall-clock %.1f s\n", t_total)

    out = joinpath(dirname(@__DIR__), "out", "shard_crosscheck")
    mkpath(out)
    npzwrite(joinpath(out, "endpoints.npz"),
             Dict("shards" => Float64.([r.shard for r in results]),
                  "r50_fit" => fits, "z0_err" => errs, "r50_cross" => crosses,
                  "w" => [r.w for r in results],
                  "nev" => Float64.([r.nev for r in results])))
    open(joinpath(out, "results.toml"), "w") do io
        TOML.print(io, Dict(
            "n_shards" => length(results),
            "selection" => PARAMS.truth_selection, "niter" => PARAMS.niter,
            "window_mm" => collect(ref.window),
            "sens" => Dict("cache" => SENS_CACHE,
                           "recocrysp_sha" => meta["recocrysp_sha"]),
            "per_shard" => Dict(
                "shard" => [r.shard for r in results],
                "n_events" => [r.nev for r in results],
                "r50_fit_mm" => fits, "z0_err_mm" => errs,
                "r50_crossing_mm" => crosses, "w_mm" => [r.w for r in results]),
            "sigma_R_fit" => Dict("mean_mm" => sf.mean, "sigma_mm" => sf.sigma,
                                  "sem_mm" => sf.sem, "offset_mm" => sf.offset),
            "sigma_R_crossing" => Dict("mean_mm" => sc.mean, "sigma_mm" => sc.sigma,
                                       "sem_mm" => sc.sem, "offset_mm" => sc.offset),
            "z0_err_calibration" => calib,
            "reference" => Dict("dose_R80_mm" => ref.dose_R80,
                                "activity_R50_fit_mm" => ref.activity_R50_fit,
                                "activity_R50_crossing_mm" => ref.activity_R50),
            "timing_s" => t_total))
    end
    println("wrote $(joinpath(out, "results.toml")) (+ endpoints.npz)")
end

main()
