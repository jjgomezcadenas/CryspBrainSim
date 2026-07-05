# CryspBrainSim — development plan

The reconstruction + range-analysis repo: read the list-mode LORs a PET scanner would record
for a proton SOBP field, reconstruct the β⁺ activity, extract the distal range endpoint R, and
report its statistical precision σ_R vs dose — one curve per scanner geometry.

This document is the build plan. It points at the upstream contracts and records what is ours to
decide and build.

## Upstream contracts (the reference set)

Vendored snapshots live in `dev/reference/` (provenance and SHAs in `dev/reference/README.md`).

- **Products tree** — `dev/reference/PRODUCTS.md` (the `PtCryspProds/` layout and file contract).
- **Data-generation strategy** — `dev/reference/data_generation_strategy.md` (shards vs
  realizations, the `thin_lm` spec, the ten-shard cross-check).
- **Endpoint method** — `dev/reference/range_verification_recipe.md` (the σ_R recipe).
- **LOR schema** — `dev/reference/SCHEMA.md`.
- **RecoCrysp usage** — `dev/reference/RecoCryspUse.md` (dep setup, reader, attenuation,
  sensitivity, MLEM).

The two inputs are the whole contract: a `lors_shardNNN.h5` (the coincidence list) and the
scenario phantom (`phantom_regions.csv` + `material_*_meta.csv` → the μ map). This repo reads
PTCryspMC's files.

## Pipeline context

```
ptcrysp-scenarios → PTCryspMC.jl → PtCryspProds/ → THIS repo
 (proton+phantom    (photon MC,    (LOR shards,     (pool+thin → recon → profile
  → annihilations)   detector)      per scanner)      → erfc endpoint → σ_R)
```

Reconstruction machinery (projectors, MLEM, sensitivity, the HDF5 reader, `ellipsoid_chord`)
comes from **RecoCryspTools** as a pinned dependency. This repo owns products-tree navigation, the
μ-map builder from the scenario files, the Bernoulli thinning, the depth profile, the endpoint
estimator (ported to Julia), the sweep drivers, and — later — the φ-gap dual-head sensitivity
sampler.

---

## Dependencies

`Project.toml` depends on **RecoCryspTools** (which re-exports the RecoCrysp core), pinned by SHA.
Per `RecoCryspUse.md`, list two `PackageSpec`s on the one repo at one rev — the unregistered core
explicitly alongside the subdir package:

```julia
SHA = "<pinned merge commit on main>"        # tools-split is merged & published; pin it
Pkg.add([
    PackageSpec(url="https://github.com/jjgomezcadenas/RecoCrysp.git", rev=SHA),
    PackageSpec(url="https://github.com/jjgomezcadenas/RecoCrysp.git", rev=SHA,
                subdir="RecoCryspTools"),
])
```
`using RecoCryspTools` brings engine + tools. Depend on RecoCryspTools alone — it is the stable
surface RecoCrysp maintains for us.

Other deps: `LsqFit`, `SpecialFunctions` (the endpoint fit), `HDF5`, `TOML`, `NPZ` (sensitivity
cache), optional `Metal` (GPU iteration loops). Julia ≥ 1.11 (both repos on 1.12).

Record the RecoCrysp SHA in every cached artifact (sensitivity file, sweep outputs) so each result
carries the exact engine that produced it. Confirm the two PackageSpec URLs resolve at the pinned
rev before writing them into `Project.toml`.

---

## Repository structure

```
CryspBrainSim/
  Project.toml                deps (RecoCryspTools @SHA, LsqFit, SpecialFunctions, HDF5, TOML, NPZ, Metal)
  dev/PLAN.md                 this file
  src/
    CryspBrainSim.jl          module top; re-exports the public surface
    products.jl               PtCryspProds navigation: glob shards, pool, verify provenance attrs + truth/ bundle
    qa.jl                     shard statistics + sanity numbers (struct drivers assert on)
    characterize.jl          truth/ bundle → activity(z), dose(z), dose-R80/activity-R50 offset (the reference)
    mumap.jl                  scenario phantom → attenuation (ellipsoid_chord and/or voxel μ-map)
    sensitivity.jl            ContinuousPET sampler + chunked sens; cache `base` + provenance
    thinning.jl               thin_lm (pooled Bernoulli, own seed namespace)
    profile.jl                depth_profile + distal_window (fixed ROI on the beam axis)
    endpoint.jl               fit_endpoint + sigma_R  (Julia port of the py estimator)
    dualhead_sampler.jl       φ-gap sensitivity sampler — arrives with the dual-head data (see Deferred)
  config/                     the FROZEN knobs, common-mode across arms (see below)
  drivers/
    one_shard.jl              one shard → R          (the single-shard chain; the atomic unit)
    shard_crosscheck.jl       ten shards → bias-free σ_R at top dose  (the standing gate)
    sweep.jl                  thinned realizations × doses × arms → σ_R-vs-dose  (the deliverable)
  tools/                      Python quick-looks + plotting (reproducible, on disk)
    plot_shard.py             the 3×3 detector QA panel (adapted from PTCryspMC plot_prod.py)
    plot_truth.py             activity / depth_dose / dose_activity / sobp_plateau from truth/
    shard_summary.py          attrs + truth breakdown; run first on every new shard
    origin_profile.py         truth-origin depth histogram (detected-subset cross-check)
  test/runtests.jl            synthetic self-tests (tolerance-based) + Python cross-validation
  latex/                      update depth_profile.tex code-map table Python → Julia
```

Julia package name: `CryspBrainSim` (matches the repo). The analysis is Julia end to end; Python
serves quick-looks and figures. The existing `py/` modules serve as the frozen cross-validation
reference for the port, then retire.

---

## Code inventory — consume vs write

### Consume from RecoCryspTools
- `read_coincidences(path)` → `MCCoincidences`; `endpoints(c, mask)`, `is_true`/`is_scatter`/
  `is_random`, `c.origin` (true annihilation points). Rescales Int16 storage and drops degenerate
  LORs, reporting the count.
- `ellipsoid_chord(p1, p2; axes, center)` — analytic attenuation chord for the uniform head.
- `joseph3d_fwd` / `_back`, `sensitivity_image(...; weights, scale)`, `ListmodePoissonModel`,
  `mlem` / `osem` / `bsrem`, `ContinuousPET`, `sample_lors`, `gaussian_blur`.
- (Later) `TOFParameters`, `tofbin_from_origin`, the TOF projectors.

### Write here
1. **products.jl** — glob `lors_shard*.h5` under a config leaf, pool them, and verify the full
   provenance attr set on read (scenario, scanner, crystal, budget, dose_Gy, master_seed,
   realization, n_phi/n_z, windows, nevents, quantization scales), failing loudly if any is absent.
   Log the per-shard degenerate-drop count. Read `scanner_geometry.json`, the `phantom/` files, and
   the `truth/` bundle from the products tree — it is the self-describing interface. (The `truth/`
   bundle is delivered — contract and sanity numbers in `dev/upstream_response_truth_bundle.md`;
   the tree is the only input, with no reach into `ptcrysp-scenarios/`.)
2. **qa.jl** — compute a shard's statistics as a struct the drivers assert on: nrows, truth
   fractions, acceptance = nrows/nevents, per-hit energy and radius ranges, DOI range (r − R_inner),
   Δt/τ occupancy, source bounding box. Turns "the shard looks fine" into an automatable gate; the
   3×3 figure is `tools/plot_shard.py`. Runs on a shard at full statistics (distinct from a thinned
   realization). Run first on every new shard.
3. **characterize.jl** — read the `truth/` bundle and lock the reference before reconstruction:
   compose the true activity(z) from `activity_profile_<budget>.csv` (per-isotope + total; detector-
   independent), read dose(z) from `depth_dose.csv`, and compute **dose-R80** (distal 80%-of-max of
   the dose) and **true activity-R50** (half-height of the activity edge) and their **offset** — the
   ground-truth activity-to-range distance the reconstructed edge is scored against. Figures
   (activity, depth_dose, dose_activity, sobp_plateau) come from `tools/plot_truth.py`. The offset,
   the nominal activity edge, and the target-box depths feed the fixed distal window in `config/`.
4. **mumap.jl** — build attenuation from the scenario files: `phantom_regions.csv` supplies the
   solid, semi-axes, centre, and material; `material_*_meta.csv` supplies μ at 511 keV. For the
   uniform ellipsoid, call `ellipsoid_chord` with the scenario's axes, centre, and μ. Provide the
   voxel-μ-map route (`exp.(-joseph3d_fwd(...))`) as well — it is the path the multi-region head
   phantoms take later. Read μ and the ellipsoid centre from the files at runtime. Size the μ-map
   grid to cover the whole head (y down to −117 mm) while the activity grid covers the beam
   corridor plus margins — the two grids are independent and each is sized to its job.
5. **sensitivity.jl** — build `ContinuousPET(diameter, afov)` from the ring (r_inner 387 mm →
   D 774 mm, half_length 512 mm → afov 1024 mm) and accumulate the **unscaled** `base = Aᵀ(a_geom)`
   over chunked `sample_lors` draws (surface chords, per the recipe). Cache `base` to NPZ with a
   provenance record (RecoCrysp SHA, scanner geometry, μ/phantom params, n_sens, seed, grid). Apply
   the per-realization scale `n_events/n_sens` at reconstruction time — `base` is
   realization-independent and reused across the whole sweep, the scale rides on each realization's
   event count.
6. **thinning.jl** — `thin_lm(shard_files, target_counts, realization_index)`: pool across all
   shards, keep each event independently with probability `p = target/M_total` over the union, using
   an RNG seeded by the downstream realization index (its own namespace, operating purely on the
   produced LORs). Keeping events independently makes the realization count fluctuate as
   Binomial(M_total, p) ≈ Poisson(target) — that fluctuation is part of σ_R. Set
   `target_counts = round(p·M_total)` with `p = dose/top_dose`: the shards are the master at the top
   dose, so scaling by the dose ratio anchors the count physically. (Confirm `p = dose/top_dose`
   against the recipe's `dose_to_counts`; adopt an explicit yield model there if the recipe intends
   one.)
7. **profile.jl** — `depth_profile(image; voxsize, beam_axis, roi_radius, roi_centre, z_origin)`:
   integrate (sum) over a fixed transverse disc at each depth slice to preserve the counting
   statistics the endpoint fit rides on; centre the disc on the beam axis (0,0). `distal_window(
   z_edge_nominal, δ_prox, δ_dist)` brackets the falloff from the nominal edge (≈ z −5 mm for
   `uniform_headep`), fixed across arms and realizations.
8. **endpoint.jl** — Julia port of `depth_profile.py`'s windowed, Poisson-weighted erfc fit
   (`LsqFit` + `SpecialFunctions`) returning R50/R80/R20 + z0_err, plus `sigma_R` (mean/std/sem over
   the finite endpoints, counting the dropped fits). Match the covariance convention to scipy's
   `absolute_sigma=True` with σ=√max(P,1): weight each point by √counts and confirm LsqFit's `vcov`
   returns inv(JᵀWJ) directly, so z0_err carries the counting-statistics scale.
9. **dualhead_sampler.jl** — φ-gap partial-ring sensitivity sampler, feeding `sensitivity_image` the
   same way the closed-ring sampler does. Build it when the upstream engine produces the dual-head
   arm.

---

## Fixed knobs (`config/`) — the common-mode discipline

Identical across every scanner, crystal, dose, and realization; consumed unchanged by
`sensitivity.jl`, `profile.jl`, and `sweep.jl`. Keeping them common-mode is what makes σ_R
differences purely geometric.

| Knob | Value | Set at |
|---|---|---|
| voxel grid (activity) | TBD; z-sampling ~1.5–2 mm for the edge | single-shard stage |
| μ-map grid | covers whole head (y to −117 mm) | single-shard stage |
| ROI radius / centre | ρ ≈ n·√(σ_spot²+σ_PSF²)+R_β⁺ ≈ 12–15 mm; centre (0,0) | single-shard stage |
| distal fit window | (z_nom − δ_prox, z_nom + δ_dist), z_nom ≈ −5 mm | single-shard stage |
| MLEM iteration count | on the R50-vs-iterations plateau | single-shard stage |
| n_sens | 5×10⁸, revalidated at our grid | single-shard stage |
| truth selection | trues-only (truth==0), first pass | fixed |

Run plain `mlem` at the fixed iteration count on the headline pass, starting from
`x0 = Float32.(sens .> 0)` (uniform inside the FOV).

---

## Validation ladder

Each rung reuses the previous one; the chain grows without rebuilding.

1. **Starting-point characterization** (`characterize.jl`) — read `truth/`, produce the four truth
   figures, and lock the reference: **true activity-R50**, **dose-R80**, and their offset. Runs
   before any reconstruction; it is the "check the starting point" gate and the reference every later
   R is scored against. Also feeds the fixed distal window (nominal edge, target-box depths).
2. **Shard QA** (`qa.jl` + `tools/plot_shard.py`) — the statistics/sanity struct + 3×3 panel on each
   stored shard: acceptance, truth fractions, energy/radius/DOI ranges, source fill, Δt/τ. The gate
   that clears a shard before it enters the pipeline.
3. **Synthetic self-tests** — the Julia port reproduces the two `py` self-tests (erfc edge → R50)
   on shared input arrays (save the numpy-generated data and feed both implementations the same
   arrays). Compare with documented tolerances.
4. **Reconstruction-free reference** — fit the endpoint on the truth-origin depth profile
   (`c.origin`, no reconstruction). This isolates the estimator from the recon; being the detected
   subset, it cross-checks against — not replaces — the `truth/` activity-R50 from rung 1.
5. **Single shard** (`one_shard.jl`) — the full chain on shard 0 at full statistics. Acceptance:
   (a) the reconstructed profile overlays sensibly on `depth_dose.csv` and its activity-R50 matches
   the rung-1 truth activity-R50; (b) the erfc fit converges with sub-voxel z0_err in the fixed
   window; (c) R50 holds under small ROI/window perturbations; (d) wall-clock per reconstruction is
   measured (it sizes the sweep). This stage **freezes the knobs**.
6. **Ten-shard cross-check** (`shard_crosscheck.jl`) — fit each of the 10 shards independently and
   take the std → a bias-free σ_R at top dose (~24% precision, n=10). The thinned σ_R at top dose
   agrees within that band. Wire this gate in from the start; it runs once shards 1–9 land.
7. **Thinned sweep** (`sweep.jl`) — Z≈100–200 realizations × dose grid × arms → σ_R-vs-dose.
   Thin the same realization index across arms from matched shards to keep the source common-mode
   exact.

---

## Stage checks that set the knobs (single-shard stage)

- **Sensitivity noise at our grid.** The 5×10⁸ LOR count was certified at 2.5 mm voxels; at a finer
  grid each voxel is crossed by fewer sampled LORs, so confirm the sensitivity's Monte-Carlo mottle
  sits below the other noise at our voxel size before freezing n_sens — compute `base` with two
  seeds (or double n_sens) and confirm R50 is stable. Minutes on Metal.
- **Tolerance-based comparison.** The atomic Float32 backprojection accumulates in
  order-dependent fashion, so images and R50 reproduce to a small tolerance while the thinned event
  list reproduces bit-for-bit. Write every test and the cross-check against a documented tolerance
  and record its value.
- **Iteration count.** MLEM is semi-convergent, and the erfc edge width (hence R50) tracks the
  iteration number. Sweep R50 vs iterations, pick a value on the stable plateau, and freeze it
  common-mode.

---

## Deferred register (recorded so the next instance folds it in cheaply)

- **Randoms scale as activity².** Bernoulli thinning at probability f reproduces dose f·D_top for
  the activity-linear classes (trues, same-annihilation scatters). Randoms follow rate ≈ 2τS², so a
  realization thinned to f should carry randoms ∝ f² — thin `truth==2` events with p² while
  trues/scatters get p (or regenerate the randoms pass per dose upstream). Fold this in when the
  scatter/randoms axis is added, or when a scenario/budget pushes randoms above the sub-percent
  level. Today they sit at 0.16% at top dose and the first pass is trues-only.
- **Dual-head and mixed-crystal arms** arrive when the simulator gains `phi_gaps`/planar geometry
  and per-block crystals; `dualhead_sampler.jl` arrives with them. The closed-ring size variants
  (1 m / head / children) and the crystal axis (bgo/csi) are producible now.
- **Multi-region head phantoms** (`mird_head`, `headep`) arrive with the upstream multi-region
  scenario reader; the voxel-μ-map route in `mumap.jl` is the path that carries them.
- **TOF** is a cheap later axis: the same sensitivity image serves it, and per-event bins synthesize
  from the truth origin. Compare TOF and non-TOF at matched resolution (converge both, or match with
  a common post-filter), since TOF converges faster. It follows the headline result.

---

## First concrete target (data on disk today)

`PtCryspProds/uniform_headep_sobp_1e8/crysp_ring_1m/bgo/fast_1Gy/lors_shard000.h5`:
80.18 M decays → 17.43 M LORs (81.1% true / 18.7% scatter / 0.16% random). Head ellipsoid
(semi-axes 72/87/102 mm) at centre (0, −30, 0); beam and target on the axis (0,0); nominal distal
edge z ≈ −5 mm; brain μ(511) = 0.009913 mm⁻¹.

All ten BGO shards (000–009) and the `truth/` bundle are on disk; the `csi/` arm is a further
production run into the same tree, adding no machinery here. Ladder rungs 1–2 have run (locked
reference: dose-R80 −5.58 mm, activity-R50 −16.45 mm, offset −10.87 mm; all shards pass QA), and
rungs 4–6 are unblocked.

## Build order

1. Scaffold `Project.toml` (pin the RecoCrysp SHA, confirm the URLs), the module, `test/runtests.jl`.
2. Port the endpoint estimator (`endpoint.jl`, `profile.jl`) with the synthetic self-tests and the
   Python cross-validation (ladder rungs 3–4). Independent of the products data.
3. `products.jl`, `qa.jl`, `characterize.jl`; `tools/plot_shard.py`, `plot_truth.py`,
   `shard_summary.py`, `origin_profile.py` — the statistics + starting-point gates (rungs 1–2).
4. `mumap.jl`, `sensitivity.jl`.
5. `thinning.jl`.
6. `drivers/one_shard.jl` → freeze the knobs (ladder rung 5 and the stage checks above).
7. `drivers/shard_crosscheck.jl` and `drivers/sweep.jl` (rungs 6–7) as the remaining data lands.
8. Update `latex/depth_profile.tex` (code-map table → Julia; keep the single windowed estimator).
