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
# → sum the image over the whole transverse plane at each depth (the settled
# protocol; `[roi]` carries no radius) → fit the distal edge with an erfc,
# free baseline. R is read two ways — the fit's half-height (the primary,
# noise-averaged value) and the profile's half-height crossing (a direct read).

# Short filename tag for a dose: 1.0 → "1Gy", 0.5 → "0p5Gy", 0.05 → "0p05Gy".
dose_tag(d::Real) =
    (d == round(d) ? string(Int(round(d))) : replace(string(d), "." => "p")) * "Gy"

"""
    load_run_context(; products_root, scenario, scanner, crystal, leaf,
                     sens_cache, params=load_run_parameters()) -> NamedTuple

Load the fixed inputs one PET configuration reconstructs against. `crystal` is
the material key (e.g. `"BGO"`); `topology` is `"closed"` or `"open"`. Returns a
NamedTuple with the inputs and the configuration's identity (the latter names
the output directory via [`config_out`](@ref)):

- `ref`: the truth reference (dose-R80, the true activity edge, the fixed fit
  window), from [`characterize`](@ref).
- `base`, `meta`: the cached sensitivity image (checked to match the frozen
  grid) and its provenance.
- `phantom`: the head's attenuation parameters, from [`phantom_attenuation`](@ref).
- `files`: the shard files under the config leaf, in shard-index order.
- `params`: the frozen run parameters.
- `scenario`, `topology`, `ring`, `crystal`: the identity — `ring` is the
  scanner name, `crystal` the material+thickness label (e.g. `"bgo_3X0"`,
  computed from the ring wall).

The compute device is the caller's to choose (pass it to
[`reconstruct_endpoint`](@ref)), so this stays free of any GPU dependency.
"""
function load_run_context(; products_root::AbstractString, scenario::AbstractString,
                          topology::AbstractString="closed",
                          scanner::AbstractString, crystal::AbstractString,
                          leaf::AbstractString, sens_cache::AbstractString,
                          params=load_run_parameters())
    scen = joinpath(products_root, scenario)
    files = shard_files(leaf_dir(products_root; scenario=scenario, scanner=scanner,
                                 crystal=crystal, leaf=leaf))
    # v2 shards centre the phantom on the tumour; shift the truth reference by the
    # stamped offset so its window/edges land in the reconstructed image's frame.
    zoff = let a = shard_attrs(files[1])
        shard_generation(a) == "v2" ? Float64(a["source_z_offset_mm"]) : 0.0
    end
    ref = characterize(scen; z_offset_mm=zoff)
    base, meta = load_sensitivity(sens_cache)
    g = meta["grid"]
    (g["n"] == collect(params.grid.n) &&
     Float32.(g["img_origin"]) == collect(params.grid.img_origin) &&
     Float32.(g["voxsize"]) == collect(params.grid.voxsize)) ||
        error("load_run_context: sensitivity cache grid ≠ frozen run-parameter grid")
    ph = phantom_attenuation(scen)
    geo = scanner_geometry(joinpath(scen, scanner))
    return (ref=ref, base=base, meta=meta, phantom=ph, files=files, params=params,
            geo=geo, scenario=scenario, topology=topology, ring=scanner,
            crystal_material=crystal, crystal=crystal_label(crystal, geo.wall_mm))
end

"""
    write_descriptors(ctx) -> (geometry_path, crystal_path)

Stamp the self-describing scanner descriptors for `ctx`'s configuration into
the output tree: `geometry.toml` at the ring tier and `crystal.toml` at the
crystal tier. Idempotent — the drivers call it so a reader of `out/` can
answer "what scanner produced this?" from the results alone.
"""
function write_descriptors(ctx)
    g = write_ring_geometry(ctx.geo; scenario=ctx.scenario,
                            topology=ctx.topology, ring=ctx.ring)
    # Detector response from the shard attributes (constant across a config's
    # shards): energy + spatial resolution, energy cut, coincidence window.
    a = shard_attrs(ctx.files[1])
    det = (energy_resolution_fwhm=a["eres"], sigma_xyz_mm=a["sigma_xyz_mm"],
           emin_keV=a["emin_keV"], tau_ns=a["tau_ns"])
    c = write_crystal_spec(; scenario=ctx.scenario, topology=ctx.topology,
                           ring=ctx.ring, crystal=ctx.crystal,
                           material=ctx.crystal_material, wall_mm=ctx.geo.wall_mm,
                           detector=det)
    return (g, c)
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
Scales the sensitivity to this list's count, runs MLEM, sums the image over
the whole transverse plane into a depth profile, and fits the distal edge. Pass `device =
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
