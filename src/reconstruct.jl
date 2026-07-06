# reconstruct.jl — the reconstruction chain that turns one list of
# coincidences into a range endpoint, plus the run context that supplies its
# fixed inputs. This is the shared analysis surface the σ_R drivers build on
# (drivers/sigma_r_at_dose.jl, sigma_r_sweep_dose.jl).
#
# Vocabulary used here and in the drivers:
#
#   coincidence (LOR)  one recorded event — a pair of 511 keV annihilation
#                      photons caught in time coincidence, which defines a line
#                      of response (LOR) through the head.
#   shard              one complete simulated acquisition: every coincidence the
#                      scanner records for the proton field at the nominal (top)
#                      dose. Ten shards exist, each an independent repeat of the
#                      same source with a different random seed.
#   pooled master      the shards concatenated — the largest count sample
#                      available, and the source each realization draws from.
#   dose               the delivered radiation dose, in Gy. The β⁺ activity, and
#                      so the coincidence count, is proportional to it.
#   realization        one emulated acquisition, built by randomly keeping a
#                      fraction of the pooled coincidences (thinning). A valid
#                      stand-in for a fresh acquisition with no new simulation;
#                      the kept fraction sets the emulated dose.
#   R (R50)            the range endpoint — the depth of the distal activity
#                      edge, read from the reconstructed depth profile.
#   σ_R                the standard deviation of R across realizations: the
#                      statistical precision of the range measurement at one dose.
#
# The chain is frozen in config/run_parameters.toml and identical for every run:
# scale the sensitivity to the list's count → MLEM at the frozen iteration count
# → sum the image over the fixed disc ROI at each depth → fit the distal edge
# with an erfc. R is read two ways — the fit's half-height (the primary,
# noise-averaged value) and the profile's half-height crossing (a direct read).

# Short filename tag for a dose: 1.0 → "1Gy", 0.5 → "0p5Gy", 0.05 → "0p05Gy".
dose_tag(d::Real) =
    (d == round(d) ? string(Int(round(d))) : replace(string(d), "." => "p")) * "Gy"

"""
    load_run_context(; products_root, scenario, scanner, crystal, leaf,
                     sens_cache, params=load_run_parameters()) -> NamedTuple

Load the fixed inputs one scanner configuration reconstructs against, and
return them as `(ref, base, meta, phantom, files, params)`:

- `ref`: the truth reference (dose-R80, the true activity edge, the fixed fit
  window), from [`characterize`](@ref).
- `base`: the cached sensitivity image, checked to match the frozen grid.
- `meta`: the sensitivity provenance record.
- `phantom`: the head's attenuation parameters, from [`phantom_attenuation`](@ref).
- `files`: the shard files under the config leaf, in shard-index order.
- `params`: the frozen run parameters.

A different scanner has its own tree leaf and sensitivity cache — point the
arguments there. The compute device is the caller's to choose (pass it to
[`reconstruct_endpoint`](@ref)), so this stays free of any GPU dependency.
"""
function load_run_context(; products_root::AbstractString, scenario::AbstractString,
                          scanner::AbstractString, crystal::AbstractString,
                          leaf::AbstractString, sens_cache::AbstractString,
                          params=load_run_parameters())
    scen = joinpath(products_root, scenario)
    ref = characterize(scen)
    base, meta = load_sensitivity(sens_cache)
    g = meta["grid"]
    (g["n"] == collect(params.grid.n) &&
     Float32.(g["img_origin"]) == collect(params.grid.img_origin) &&
     Float32.(g["voxsize"]) == collect(params.grid.voxsize)) ||
        error("load_run_context: sensitivity cache grid ≠ frozen run-parameter grid")
    ph = phantom_attenuation(scen)
    files = shard_files(leaf_dir(products_root; scenario=scenario, scanner=scanner,
                                 crystal=crystal, leaf=leaf))
    return (ref=ref, base=base, meta=meta, phantom=ph, files=files, params=params)
end

"""
    lor_attenuation(ctx, xs, xe) -> Vector{Float32}

Per-LOR survival factors `exp(-μ · chord)` through the context's head
ellipsoid, for the coincidences with endpoints `xs`/`xe`.
"""
lor_attenuation(ctx, xs, xe) =
    attenuation_ellipsoid(xs, xe; semi_axes=ctx.phantom.semi_axes,
                          centre=ctx.phantom.centre, mu_mm_inv=ctx.phantom.mu_mm_inv)

"""
    reconstruct_endpoint(ctx, xs, xe, mult; device=identity) -> NamedTuple

Run the frozen reconstruction chain on one list of coincidences (`xs`/`xe`
endpoints in mm, `mult` the per-LOR attenuation) and read its range endpoint.
Scales the sensitivity to this list's count, runs MLEM, sums the image over the
disc ROI into a depth profile, and fits the distal edge. Pass `device =
MtlArray` (with Metal loaded) to run on the GPU. Returns
`(nev, r50_fit, z0_err, w, r50_cross)`: the count, R in both conventions, the
fitted edge width, and the fit's own error estimate.
"""
function reconstruct_endpoint(ctx, xs, xe, mult; device=identity)
    p = ctx.params
    nev = size(xs, 2)
    sens = scaled_sensitivity(ctx.base, nev, p.n_sens)
    model = ListmodePoissonModel(device(xs), device(xe), device(sens);
                                 img_origin=p.grid.img_origin, voxsize=p.grid.voxsize,
                                 mult=device(mult))
    x = mlem(model, device(Float32.(sens .> 0)); niter=p.niter)
    z, prof = depth_profile(Array(x); voxel_size_mm=p.grid.voxsize, beam_axis=3,
                            roi_radius_mm=p.roi.radius_mm, z_origin_mm=p.grid.img_origin[3])
    fit = fit_endpoint(z, prof; window=ctx.ref.window, weighted=true)
    return (nev=nev, r50_fit=fit.z0, z0_err=fit.z0_err, w=fit.w,
            r50_cross=windowed_crossing(z, prof, ctx.ref.window))
end
