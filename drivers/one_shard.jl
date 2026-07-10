# drivers/one_shard.jl — the single-shard chain (build step 5, ladder
# rung 5): reconstruct one stored shard's trues at full statistics on the
# corridor grid, profile over the whole transverse plane (the settled
# protocol), fit the distal erfc endpoint in the fixed window, and score R
# against the rung-1 truth reference. Sweeps R50 vs MLEM iteration
# (semi-convergence plateau), perturbs the window (stability), and records
# the wall-clock that sizes the sweep.
#
# Run:  julia -t auto --project=. drivers/one_shard.jl [shard_index] [--all-uncorr]
# Default selection is the frozen trues-only run parameter. `--all-uncorr` is the
# scatter-effect study (step a): let scatters and randoms into the event list
# with NO correction — same run parameters, same model otherwise — and tag the outputs
# `_all_uncorr` so the two runs compare side by side.
# Writes out/one_shard/results_<tag>.toml, profile + image NPZ, and the
# R50-vs-iteration table; figures come from tools/plot_one_shard.py.

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf
using TOML

const SHARD = let a = filter(!startswith("--"), ARGS)
    isempty(a) ? 0 : parse(Int, a[1])
end
const ALL_UNCORR = "--all-uncorr" in ARGS

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")

# The frozen run parameters (config/run_parameters.toml). NITER_MAX runs past
# the frozen iteration count so every run re-verifies the plateau.
const PARAMS = load_run_parameters()
const CFG = PARAMS.config      # the active arm (run_parameters.toml)
const SCENARIO, TOPOLOGY, RING = CFG.scenario, CFG.topology, CFG.scanner
const N = PARAMS.grid.n
const VS = PARAMS.grid.voxsize
const ORG = PARAMS.grid.img_origin

const ROI_MM = PARAMS.roi.radius_mm    # nothing = whole plane (the settled protocol)
const NITER_MAX = 2 * PARAMS.niter
const CHECK_EVERY = 10

fit_r50(img, window; roi=ROI_MM) = begin
    z, prof = depth_profile(img; voxel_size_mm=VS, beam_axis=3,
                            roi_radius_mm=roi, z_origin_mm=ORG[3])
    (z=z, prof=prof, fit=fit_endpoint(z, prof; window=window, weighted=true))
end

function main()
    # --- inputs via the shared context loader (reference, sensitivity cache
    # grid-verified, phantom, shard files, and the configuration identity).
    cache = joinpath(sensitivity_out(SCENARIO, TOPOLOGY, RING),
                     sensitivity_cache_name(PARAMS))
    ctx = load_run_context(; products_root=ROOT, scenario=SCENARIO,
                           topology=TOPOLOGY, scanner=RING, crystal=CFG.crystal,
                           leaf=CFG.leaf, sens_cache=cache, params=PARAMS)
    write_descriptors(ctx)          # stamp geometry.toml + crystal.toml
    ref, base, meta, ph = ctx.ref, ctx.base, ctx.meta, ctx.phantom
    n_sens = PARAMS.n_sens
    show(stdout, MIME"text/plain"(), ref); println()

    file = ctx.files[SHARD+1]
    println("shard: ", basename(file))
    r = read_shard(file)
    tmask = ALL_UNCORR ? trues(length(r.coinc)) : is_true(r.coinc)
    selection = ALL_UNCORR ? "all-uncorrected" : PARAMS.truth_selection
    xs, xe = endpoints(r.coinc, tmask)
    nev = size(xs, 2)
    println("selection: $selection ($(nev) events)")

    # Grid coverage of the true origins (edges = voxel centres ± half voxel).
    o = r.coinc.origin[:, tmask]
    lox, hix = ORG[1] - VS[1] / 2, ORG[1] + (N[1] - 1 + 0.5f0) * VS[1]
    loz, hiz = ORG[3] - VS[3] / 2, ORG[3] + (N[3] - 1 + 0.5f0) * VS[3]
    out_frac = count(i -> !(lox <= o[1, i] <= hix && lox <= o[2, i] <= hix &&
                            loz <= o[3, i] <= hiz), 1:nev) / nev
    @printf("events: %d;  origins outside grid: %.4f%%\n", nev, 100out_frac)

    # --- model: attenuation mult + scaled sensitivity, on Metal when present
    t_att = @elapsed mult = attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes,
                                                  centre=ph.centre,
                                                  mu_mm_inv=ph.mu_mm_inv)
    sens = scaled_sensitivity(base, nev, n_sens)
    dev = Metal.functional() ? MtlArray : identity
    model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                 img_origin=ORG, voxsize=VS, mult=dev(mult))
    x = dev(Float32.(sens .> 0))

    # --- MLEM with R50 checkpoints; the image at the FROZEN iteration count
    # is the headline, the run continues to 2× to re-verify the plateau.
    iters = Int[]
    r50s = Float64[]
    werrs = Float64[]
    img = zeros(Float32, N)
    t_recon = @elapsed for it in CHECK_EVERY:CHECK_EVERY:NITER_MAX
        x = mlem(model, x; niter=CHECK_EVERY)
        res = fit_r50(Array(x), ref.window)
        push!(iters, it); push!(r50s, res.fit.z0); push!(werrs, res.fit.z0_err)
        it == PARAMS.niter && (img = Array(x))
        @printf("iter %3d: R50 = %8.3f ± %.3f mm  (w = %.2f mm)%s\n",
                it, res.fit.z0, res.fit.z0_err, res.fit.w,
                it == PARAMS.niter ? "   <- frozen niter" : "")
    end

    # --- final numbers at the frozen iteration: fit, crossing, stability
    fin = fit_r50(img, ref.window)
    cross = windowed_crossing(fin.z, fin.prof, ref.window)
    stab = Dict{String,Float64}()
    for dz in (-2.0, 2.0)
        w = (ref.window[1] + dz, ref.window[2] + dz)
        stab["window_$(dz)mm"] = fit_r50(img, w).fit.z0
    end
    spread = maximum(abs.(collect(values(stab)) .- fin.fit.z0))

    @printf("\nfinal (frozen iter %d): R50 fit = %.3f ± %.3f mm | crossing = %.3f mm\n",
            PARAMS.niter, fin.fit.z0, fin.fit.z0_err, cross)
    @printf("reference:  truth fit = %.3f mm | truth crossing = %.3f mm | dose-R80 = %.3f mm\n",
            ref.activity_R50_fit, ref.activity_R50, ref.dose_R80)
    @printf("stability spread (window ±2 mm): %.3f mm\n", spread)
    @printf("timing: attenuation %.2f s | %d MLEM iters (with checkpoints) %.1f s\n",
            t_att, NITER_MAX, t_recon)

    # --- artifacts
    out = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                   "one_shard")
    mkpath(out)
    tag = @sprintf("shard%03d", SHARD) * (ALL_UNCORR ? "_all_uncorr" : "")
    npzwrite(joinpath(out, "recon_$(tag).npz"),
             Dict("image" => img, "z" => collect(fin.z), "profile" => fin.prof,
                  "iters" => Float64.(iters), "r50s" => r50s, "z0_errs" => werrs))
    open(joinpath(out, "results_$(tag).toml"), "w") do io
        TOML.print(io, Dict(
            "shard" => SHARD, "selection" => selection, "n_events" => nev,
            "out_of_grid_frac" => out_frac,
            "grid" => Dict("n" => collect(N), "img_origin" => Float64.(collect(ORG)),
                           "voxsize" => Float64.(collect(VS))),
            "profile" => "whole-plane", "niter" => PARAMS.niter,
            "niter_plateau_check" => NITER_MAX,
            "window_mm" => collect(ref.window),
            "sens" => Dict("cache" => cache, "n_sens" => n_sens,
                           "recocrysp_sha" => meta["recocrysp_sha"]),
            "r50_vs_iter" => Dict("iters" => iters, "r50_mm" => r50s,
                                  "z0_err_mm" => werrs),
            "final" => Dict("r50_fit_mm" => fin.fit.z0,
                            "z0_err_mm" => fin.fit.z0_err,
                            "w_mm" => fin.fit.w,
                            "r50_crossing_mm" => cross,
                            "stability" => stab, "stability_spread_mm" => spread),
            "reference" => Dict("activity_R50_fit_mm" => ref.activity_R50_fit,
                                "activity_R50_crossing_mm" => ref.activity_R50,
                                "dose_R80_mm" => ref.dose_R80),
            "timing_s" => Dict("attenuation" => t_att, "recon" => t_recon)))
    end
    println("wrote $(joinpath(out, "results_$(tag).toml")) (+ recon_$(tag).npz)")
end

main()
