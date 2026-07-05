# Using RecoCrysp from CryspBrainSim — the full recipe

**For:** the CryspBrainSim repo (data generation + reconstruction + analysis).
**Written:** 2026-07-05, against branch `tools-split` (commit `e3e3c0b`).
**Contract:** depend on `RecoCryspTools` (which re-exports the `RecoCrysp` core);
never depend on `RecoExamples` (a thin compatibility shim) or anything under
`recoExamples/<study>/` (per-study code, churns freely).

## 1. Dependency setup (verified recipe)

Two dep entries on ONE repository at ONE pinned rev — Julia resolves `[sources]`
only for the active project, so the unregistered core must be listed explicitly
alongside the subdir package:

```julia
using Pkg
SHA = "<pinned commit>"     # record this SHA in your provenance artifacts
Pkg.add([
    PackageSpec(url = "https://github.com/jjgomezcadenas/RecoCrysp.git", rev = SHA),
    PackageSpec(url = "https://github.com/jjgomezcadenas/RecoCrysp.git", rev = SHA,
                subdir = "RecoCryspTools"),
])
```

Then a single import brings everything (engine + tools):

```julia
using RecoCryspTools
```

Requirements: Julia ≥ 1.11 (this machine runs 1.12; keep both repos on the
same version). GPU: `using Metal` and pass `MtlArray`s — nothing else changes.
NOTE: as of writing, `tools-split` is not yet merged/pushed; pin the SHA of the
merge commit once it lands on `main`.

## 2. Conventions (read once, save hours)

- **Units**: mm everywhere; `Float32` everywhere (Apple GPUs have no Float64).
- **Images**: `(n0, n1, n2)` arrays; `img_origin` is the world coordinate of the
  CENTER of the first voxel; `voxsize` a 3-tuple. Center the grid with
  `org = ntuple(i -> -(n[i]-1)/2 * vs[i], 3)`.
- **LORs/events**: columns of two `(3, N)` matrices (`xstart`, `xend`).
- **Determinism**: the back projector uses atomic Float32 adds — multithreaded
  results are order-nondeterministic at rounding level. Compare with
  tolerances, never `==` (see `test/diag_backproj_determinism.jl`).
- **Reproducibility**: pass explicit `MersenneTwister(seed)` everywhere; record
  seeds in provenance.

## 3. Reading PTCRYSP listmode

```julia
c = read_coincidences("<prod>/lors_det.h5")   # MCCoincidences
c.xstart, c.xend      # (3, N) detector hits, mm (the reconstruction LORs)
c.origin              # (3, N) true annihilation points (truth; TOF synthesis)
tmask = is_true(c); smask = is_scatter(c); rmask = is_random(c)
xs, xe = endpoints(c, tmask)                  # column-select by mask
```

The reader rescales the Int16 storage via the file's `xyz_scale_mm`/`e_scale_keV`
attributes and **drops degenerate coincidences** (identical endpoints or
non-finite coordinates) with a warning — they hit undefined behavior in the
projectors. The schema contract lives in `RecoCryspTools/src/mc_listmode.jl`:
when PTCryspMC adds provenance attrs, extend it THERE (additively — this repo's
studies read the same files).

## 4. Attenuation for the ellipsoid + brain μ

Analytic (exact for a uniform body — the pattern the NEMA studies use):

```julia
mu = 0.0096f0                                  # water-ish, mm^-1 at 511 keV
chord(i) = ellipsoid_chord((xs[1,i], xs[2,i], xs[3,i]),
                           (xe[1,i], xe[2,i], xe[3,i]);
                           axes = (a, b, cz), center = (0, 0, 0))
a_fac = Float32[exp(-mu * chord(i)) for i in 1:size(xs, 2)]
```

Non-uniform μ-map alternative (exact for any voxelized μ): `a_fac =
exp.(-joseph3d_fwd(xs, xe, mumap, org, vs))`. Compute attenuation on CPU
coordinate arrays; move results to the device afterward.

## 5. The sensitivity job (the normalization deliverable)

The MLEM sensitivity is `sens = scale · Aᵀ(a)` over a large sample of GEOMETRIC
LORs, independent of the event list. **The `scale` convention is the one subtle
knob** (read the `sensitivity_image` docstring): with an independent sample of
`n_sens` LORs and `n_events` events, `scale = n_events / n_sens`; the sample's
`1/√n_sens` Monte-Carlo noise imprints on the image, so sample BIG — the locked
recipe from the MC studies is **n_sens = 5×10⁸** (below the noise floor at
2.5 mm voxels).

```julia
sc = ContinuousPET(diameter = D_mm, afov = AFOV_mm)     # scanner surface
rng = MersenneTwister(seed)
base = zeros(Float32, n)                                # accumulate in chunks
CHUNK = 20_000_000
done = 0
while done < n_sens
    nb = min(CHUNK, n_sens - done)
    gxs, gxe = sample_lors(sc, nb; rng = rng)           # uniform detector chords
    ga = ...attenuation for gxs/gxe (Section 4)...
    base .+= sensitivity_image(gxs, gxe, n, org, vs; weights = ga)
    done += nb
end
sens = base .* Float32(n_events / n_sens)
```

Cache `base` to disk (NPZ) next to a provenance record: **RecoCrysp SHA,
scanner geometry, μ/phantom parameters, n_sens, seed, grid**. On Metal the
backprojection runs at ~170 Mlors/s → 5×10⁸ LORs is minutes, dominated by LOR
generation, not projection. DOI-aware sampling variants (`surface_doi_lors`,
`emission_sens_lors`) exist for measure studies; the φ-gap dual-head sampler is
scenario code — write it in CryspBrainSim, feeding `sensitivity_image` the same
way. TOF does NOT change `sens` (the TOF adjoint summed over bins equals the
non-TOF adjoint) — one sensitivity serves TOF and non-TOF reconstructions.

## 6. Reconstruction

```julia
model = ListmodePoissonModel(xs, xe, sens;
    img_origin = org, voxsize = vs,
    mult = a_ev,               # per-event n·a (attenuation at the event LORs)
    counts = nothing,          # 1/event (raw listmode) or multiplicities
    contamination = nothing)   # additive scatter+randoms if/when modeled
x0 = Float32.(sens .> 0)       # uniform start inside the FOV
x  = mlem(model, x0; niter = 30)                        # or osem(subset_models(...))
xr = bsrem(model, x0, RelativeDifferencePrior(β, 2f0, ε); niter = 60)
```

Guidance from the studies (don't rediscover): MLEM is semi-convergent — early
stopping or the iterate-then-post-filter recipe (`gaussian_blur`) for noise
control; for edge-preserving priors use `bsrem` (convergent, stable under
attenuation), NOT `osl_mlem` (needs its clamp, diverges at small `sens`); prior
strength β scales with the magnitude of `sens` and the activity — tune per
problem, `ε` of the RDP well below the activity scale. If an additive
background is ever modeled: normalize it as an ABSOLUTE intensity over the
sampled LOR measure, `∫c dμ = N_bg/n_prompt` (`background_intensity` /
`sino_measure3` — the construction that closed the NEMA scatter residual; see
`scatter_correction_wip.md`).

## 7. TOF

PTCRYSP files carry NO timing — synthesize per-event bins from the true origin
(σ_TOF is a study knob; −1 flags events outside the bin range, filter them with
their endpoint columns):

```julia
p  = TOFParameters(num_tofbins = T, tofbin_width = w_mm, sigma_tof = σ_mm)
tb = tofbin_from_origin(c.xstart, c.xend, c.origin, p; rng = rng)
keep = tb .>= 0
model = ListmodePoissonModel(c.xstart[:, keep], c.xend[:, keep], sens;
    img_origin = org, voxsize = vs, mult = a_ev[keep],
    tofbin = tb[keep], tof_parameters = p)              # same sens as non-TOF
```

Comparison pitfall (measured, documented in the tutorial TOF example): TOF
converges FASTER, so fixed-early-iteration noise comparisons can invert the
ranking — compare at matched resolution (converge + one common post-filter).

## 8. GPU

```julia
using Metal
to_dev = MtlArray
model = ListmodePoissonModel(to_dev(xs), to_dev(xe), to_dev(sens); ...,
                             mult = to_dev(a_ev))       # all arrays one backend
```

Throughput (M5 Pro, 128³): non-TOF 170/171 Mlors/s fwd/back; TOF listmode
150/186 Mevents/s (≈ non-TOF speed); CPU back projection does NOT scale with
threads (atomic adds) — use the GPU for iteration loops.

## 9. Where to look

- Package docs: `julia --project=docs docs/make.jl` → `docs/build/` (full API +
  five executable examples incl. TOF and penalized).
- `statusmd/pp_status.md` (core), `statusmd/reco_status.md` (Tools),
  `statusmd/tof_port_scope.md` (TOF conventions),
  `statusmd/scatter_correction_wip.md` (background normalization principle).
- Working end-to-end precedents: `recoExamples/nema/nema_la_water_bgo/compare_methods.jl`
  (sens caching, contamination modes), `tutorial/examples/tof/run.jl` (TOF
  simulation → reconstruction → matched-resolution comparison).
