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
#
# Second half (`--thinned Z`, default Z = 50): pool the ten shards, draw Z
# Bernoulli-thinned realizations at the top dose (target = one acquisition's
# count, keep-probability ≈ 1/10 over the union), reconstruct each with the
# same frozen chain, and compare the thinned σ_R against the per-shard σ_R
# stored by the first half — the rung-6 agreement gate. Writes
# results_thinned.toml + endpoints_thinned.npz.

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf
using Statistics: mean
using TOML

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const SCEN = joinpath(ROOT, "uniform_headep_sobp_1e8")

const THINNED = let i = findfirst(==("--thinned"), ARGS)
    i === nothing ? 0 : (i < length(ARGS) ? parse(Int, ARGS[i+1]) : 50)
end

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

function setup()
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
    return (ref=ref, base=base, meta=meta, ph=ph, dev=dev,
            files=shard_files(leaf))
end

# Reconstruct one event list (endpoints + per-event attenuation) with the
# frozen chain and read the endpoint in both conventions.
function recon_fit(xs, xe, mult, ref, base, dev)
    nev = size(xs, 2)
    sens = scaled_sensitivity(base, nev, PARAMS.n_sens)
    model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                 img_origin=ORG, voxsize=VS, mult=dev(mult))
    x = mlem(model, dev(Float32.(sens .> 0)); niter=PARAMS.niter)
    z, prof = depth_profile(Array(x); voxel_size_mm=VS, beam_axis=3,
                            roi_radius_mm=PARAMS.roi.radius_mm, z_origin_mm=ORG[3])
    fit = fit_endpoint(z, prof; window=ref.window, weighted=true)
    return (nev=nev, r50_fit=fit.z0, z0_err=fit.z0_err, w=fit.w,
            r50_cross=windowed_crossing(z, prof, ref.window))
end

function thinned_main()
    s = setup()
    println("pooling $(length(s.files)) shards…")
    t_pool = @elapsed pool = pool_shards(s.files)
    M_total = length(pool.coinc)
    n_shards = length(s.files)
    tmask = is_true(pool.coinc)
    t_att = @elapsed a_all = attenuation_ellipsoid(pool.coinc.xstart, pool.coinc.xend;
                                                   semi_axes=s.ph.semi_axes,
                                                   centre=s.ph.centre,
                                                   mu_mm_inv=s.ph.mu_mm_inv)
    target = dose_to_counts(1.0, 1.0, M_total, n_shards)
    @printf("pooled %d LORs in %.0f s (attenuation %.1f s); top-dose target %d (p = %.4f)\n",
            M_total, t_pool, t_att, target, target / M_total)

    results = NamedTuple[]
    t_total = @elapsed for z in 1:THINNED
        keep = thin_lm(pool.coinc, target, z) .& tmask
        xs, xe = endpoints(pool.coinc, keep)
        res = recon_fit(xs, xe, a_all[keep], s.ref, s.base, s.dev)
        push!(results, res)
        @printf("realization %2d: R50 fit %8.3f ± %.3f mm | crossing %8.3f mm | %d ev\n",
                z, res.r50_fit, res.z0_err, res.r50_cross, res.nev)
    end

    fits = [r.r50_fit for r in results]
    sf = sigma_R(fits; dose_bragg_peak=s.ref.dose_R80)
    sc = sigma_R([r.r50_cross for r in results]; dose_bragg_peak=s.ref.dose_R80)
    @printf("\nthinned σ_R at 1 Gy (Z = %d, std rel. error ≈ %.0f%%):\n",
            sf.n_ok, 100 / sqrt(2 * (sf.n_ok - 1)))
    @printf("  fit convention:      mean %8.3f mm | σ_R %.3f mm | sem %.3f mm\n",
            sf.mean, sf.sigma, sf.sem)
    @printf("  crossing convention: mean %8.3f mm | σ_R %.3f mm\n", sc.mean, sc.sigma)

    out = joinpath(dirname(@__DIR__), "out", "shard_crosscheck")
    shard_toml = joinpath(out, "results.toml")
    gate = Dict{String,Any}()
    if isfile(shard_toml)
        sh = TOML.parsefile(shard_toml)
        σ_sh = sh["sigma_R_fit"]["sigma_mm"]
        n_sh = sh["n_shards"]
        ratio = sf.sigma / σ_sh
        # Agreement band: combined relative error of two std estimates.
        band = sqrt(1 / (2 * (sf.n_ok - 1)) + 1 / (2 * (n_sh - 1)))
        pass = abs(ratio - 1) < 2 * band
        @printf("  gate vs per-shard σ_R %.3f mm (n = %d): ratio %.2f, band ±%.0f%% (2σ) → %s\n",
                σ_sh, n_sh, ratio, 200band, pass ? "AGREES" : "DISAGREES")
        gate = Dict("shard_sigma_mm" => σ_sh, "ratio" => ratio,
                    "band_2sigma" => 2band, "pass" => pass)
    else
        println("  (no per-shard results.toml — run the first half for the gate)")
    end
    @printf("  total wall-clock %.0f s\n", t_total)

    mkpath(out)
    npzwrite(joinpath(out, "endpoints_thinned.npz"),
             Dict("r50_fit" => fits, "z0_err" => [r.z0_err for r in results],
                  "r50_cross" => [r.r50_cross for r in results],
                  "nev" => Float64.([r.nev for r in results])))
    open(joinpath(out, "results_thinned.toml"), "w") do io
        TOML.print(io, Dict(
            "Z" => THINNED, "dose_Gy" => 1.0, "target_counts" => target,
            "M_total" => M_total, "n_shards" => n_shards,
            "seed_base" => THINNING_SEED_BASE,
            "sigma_R_fit" => Dict("mean_mm" => sf.mean, "sigma_mm" => sf.sigma,
                                  "sem_mm" => sf.sem, "offset_mm" => sf.offset),
            "sigma_R_crossing" => Dict("mean_mm" => sc.mean, "sigma_mm" => sc.sigma),
            "gate" => gate, "timing_s" => t_total))
    end
    println("wrote $(joinpath(out, "results_thinned.toml")) (+ endpoints_thinned.npz)")
end

function main()
    s = setup()
    ref, base, meta, ph, dev, files = s.ref, s.base, s.meta, s.ph, s.dev, s.files
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

THINNED > 0 ? thinned_main() : main()
