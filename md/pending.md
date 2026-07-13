# Pending (smaller items, independent of the isotope-washout study)

- **Ring-CsI del180 washed σ_R N=200 firm-up** — the v2 ring-CsI del180 washed point (0.199) dips
  ~2σ below trend at N=100 (the small bores are smooth); N=200 halves the band (±7%→±5%) and would
  either pull it into line or confirm it. Re-run `drivers/sigma_r_v2.jl --realizations 200` on the
  ring-CsI v2 arm (md/results.md → "Generation-2 σ_R study").
- **CHS + R35: remaining nine shards per arm** — obsolete for v2 (the v2 R35/40 arms already ship
  10 shards); kept only for the legacy single-shard CHS/R35 extrapolations of §1–7.
- **BGO all-events bore-radius sensitivity** — at matched dose the BGO all-events penalty grows
  monotonically as the bore shrinks (ring < R35 < CHS) while trues match everywhere, with near-equal
  event mix; understand the scatter background shape in compact geometries. Related: the ring BGO
  all-events 1/√dose extrapolation sits 1.6σ below the ten-shard 1 Gy measurement (counting-only
  thinning caveat).

- **`latex/cbs.tex`** (the living draft, separate from `endpoint_precision.tex`): revert eq:sigmaR
  from R_p back to z0/R50 (keep the R_p accuracy paragraph), cite Zapien-Campos, fold in the
  two-scanner numbers. NB: `cbs.tex` carries the user's own uncommitted edits — commit only on
  explicit instruction.
- **Composite-erfc edge model** (2–3 isotope components, offsets/widths frozen from per-isotope
  truth profiles, free amplitudes + global shift); adopt only if σ and rung stability improve
  incl. at 0.1 Gy. Directly relevant to isotope washout — see [isotope-washout.md](isotope-washout.md).

## Deferred / on-request

- **Washout spatial non-uniformity / redistribution** — the only open IW item and the one route to a
  genuine (non-calibratable) washout bias: a heterogeneous perfusion field, or cleared atoms decaying
  elsewhere in the FOV. Both are a downstream perfusion/compartment transport model on the production
  points (not Geant4, not the range estimator). Trigger: if the edge proves sensitive to the
  clearance field. See [isotope-washout.md](isotope-washout.md), washout-g4-formulation.md.
- **Scatter correction** — machinery exists (`recon_scatters.jl`, `scatter_profile.py`); not needed
  at present precision (calibration systematics dominate ≫ 0.1 mm). The trigger for revisiting: a
  window slope approaching the edge gradient on some future configuration.
- **Geometry axis** (part b continued) — ring length (sensitivity), open geometries (angular
  coverage); the analysis is configuration-blind, so each runs through the identical battery.
