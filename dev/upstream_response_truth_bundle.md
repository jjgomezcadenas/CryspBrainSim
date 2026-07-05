# Response from PTCryspMC — the `truth/` bundle is delivered

**Re:** `upstream_request_truth_bundle.md`.
**Status:** DONE. Implemented, tested, and backfilled into the current master.
(PTCryspMC commit `4e1554a`.)

## What was delivered

A new **`truth/` level at `PtCryspProds/<scenario>/`** — detector-independent, shared across every
scanner/crystal exactly like `phantom/`. Written by `publish_prod` (once per scenario) and available
for the existing master:

```
PtCryspProds/uniform_headep_sobp_1e8/truth/
  depth_dose.csv                      # dose(z) → dose-R80                     [copy]
  sobp_layers.csv  (+ _meta)          # the SOBP beam design                   [copy]
  run_meta.csv                        # target-box depths, Np/Gy, normalization [copy]
  sampling_budget_fast.csv (+ _meta)  # per-isotope N_expected for the mix      [copy]
  activity_profile_fast.csv (+ _meta) # binned activity(z), per isotope+total   [derived]
```

Everything requested is present. The five copies are verbatim from the scenario; `activity_profile`
is the one derived file, plus a self-describing `_meta` companion.

## `activity_profile_<budget>.csv` — how to read it

- **Columns:** `z_mm, O15, C11, N13, C10, O14, total` (isotopes in id order; per-isotope + `total`).
- **z-frame:** the **exact** `z_mm` column of `depth_dose.csv` (same bins), so activity(z) and dose(z)
  overlay directly with no re-binning. Verified identical.
- **Units / scaling:** decays = **expected (mean) counts** at the leaf's dose (`fast_1Gy` → 1 Gy),
  with escaped positrons excluded — i.e. each isotope column integrates to `N_expected · f_inside`,
  the **same source scaling the LOR shards materialize from**. So the profile composes with the pooled
  LORs directly (no rescale). Dose is a linear scale (shape is dose-invariant); the top-dose profile is
  the shape at any dose. See `activity_profile_fast_meta.csv` for the stamped scenario/budget/dose.

Sanity numbers for `uniform_headep_sobp_1e8 / fast`: O15 ∫ = 5.519e7, C11 1.907e7, N13 4.488e6,
C10 1.798e5, O14 1.246e6, **total 8.018e7** (≈ one shard's ΣM). The ~0.4% deficit vs the raw
`sampling_budget` N_expected is the escaped-positron loss — the same events the shards drop.

## Notes for the characterization step

- **activity-R50** comes from `activity_profile.total` (or per-isotope, since the isotopes have
  different distal falloffs); **dose-R80** from `depth_dose.dose_core_Gy`. Both share the z-frame, so
  the **activity-R50 → dose-R80 offset** is a direct read — the locked reference for scoring.
- Per-budget: swap `fast` for `inroom` / `offline` (different isotope mix; not thin-able between
  budgets). Only `fast` is published so far.
- Column meanings trace to the scenario `SCHEMA.md`; nothing about the LOR files or the existing tree
  changed. `PtCryspProds/README.md` (the contract, from PTCryspMC `dev/PRODUCTS.md`) now documents the
  `truth/` level in the tree, axes, and file-contract tables.

You can build directly on the products tree now — no reach into `ptcrysp-scenarios/`.
