# Request to PTCryspMC — export a per-scenario `truth/` bundle in `PtCryspProds/`

**From:** the CryspBrainSim analysis repo (reconstruction + σ_R).
**Goal:** let the analysis characterize the *starting point* (the source activity, the dose, their
relationship) and lock the truth reference the reconstruction is scored against — reading only the
products tree, with no reach into `ptcrysp-scenarios/`.

## Why

The analysis measures a **reconstructed activity edge** and reports it against two truth curves and
their offset:

- **true activity(z)** — the detector-independent β⁺ source falloff; the target a perfect
  reconstruction recovers (activity-R50).
- **dose(z)** — the physical Bragg/SOBP distal edge (dose-R80), the clinically meaningful range.

Their offset (activity-R50 to dose-R80) is the physics of the method and the reference every
reconstructed edge is compared to. It is the truth version of `dose_activity.png`.

Today `PtCryspProds/<scenario>/` exports the μ-map inputs (`phantom/`) and the LOR shards, but **not
the dose-side truth or the clean activity truth** — those live only in `ptcrysp-scenarios/`. Starting
the analysis without them means reaching into a second repo for the reference, which is the broken
pipe we want to close before building.

The h5 shards do carry `x0/y0/z0` (true annihilation points), but only for the **detected subset**,
which differs by scanner — a per-detector cross-check, not the common-mode truth. And the shards have
**no per-event isotope column**, so per-isotope activity(z) (the isotopes have different distal
falloffs) must come from the scenario, not the h5. Hence the request to export the truth explicitly.

## The ask

Extend `publish_prod` to write a **`truth/` bundle**, alongside the existing `phantom/`, carrying the
scenario characterization. It is detector-independent, so it lives at the `<scenario>/` level and is
shared across every scanner/crystal (exactly like `phantom/`). `publish_prod` already reads the
scenario, so most of these are direct copies.

```
PtCryspProds/<scenario>/
  phantom/                          (as now)
  truth/                            ← new
    depth_dose.csv                  # dose(z): the Bragg/SOBP distal edge → dose-R80   [copy]
    sobp_layers.csv  (+ _meta)      # the SOBP beam design                             [copy]
    run_meta.csv                    # target-box depths, Np/Gy, normalization          [copy]
    sampling_budget_<budget>.csv (+ _meta)   # per-isotope N_expected for the mix       [copy]
    activity_profile_<budget>.csv   # binned activity(z), per isotope + total          [derived]
```

### `activity_profile_<budget>.csv` (the one derived file)

A binned depth profile of the true β⁺ activity — the clean, detector-independent source curve the
analysis reads instead of dragging the full `emitters.csv` (2.1 M rows) into the tree.

- **Rows:** one per z-bin (suggest the same 0.8 mm binning / z-frame as `depth_dose.csv`, so the two
  overlay directly).
- **Columns:** `z_mm`, then one activity column per isotope (`O15, C11, N13, C10, O14`) and a
  `total`, in the same units and dose/budget scaling the shard uses (so it composes with the LORs).
- **How:** histogram `emitters.csv`'s `anh_z_mm` per `isotope_id`, weighted by the budget's
  `N_expected` and the dose scaling. Per-budget because the isotope mix depends on the acquisition
  timing; dose is a linear scale (shape is dose-invariant, so the top-dose profile is the shape at
  any dose).
- **Fallback:** if deriving it is inconvenient, shipping `emitters.csv` in `truth/` is acceptable —
  the analysis will bin it — but the derived profile keeps the tree small and emitter-free.

`depth_dose.csv`, `sobp_layers.csv`, `run_meta.csv`, and `sampling_budget_*` are budget/detector-
independent Stage-A physics (except the budget file, which is per-budget) and are verbatim copies.

## Contract update

Add the `truth/` level to `PRODUCTS.md`'s tree and axes table (detector-independent, `<scenario>/`
level, shared like `phantom/`), and have `publish_prod` stamp it. Columns reference the scenario
`SCHEMA.md`. Nothing about the LOR files or the existing layout changes.

## What the analysis does with it

A characterization step (run before reconstruction) reads `truth/` and produces the `activity`,
`depth_dose`, `dose_activity`, and `sobp_plateau` figures, and computes the **ground-truth
activity-R50-to-dose-R80 offset** — the locked reference the reconstructed edge is measured against.
