# Pending (smaller items)

## Data-driven / BGO current work

- **Finite-pool validation tension.** On the dd TBP reference (d120s300) the
  nominal thinned finite-pool-corrected σ_R (0.084 mm) sits ~2.7σ above the
  direct ten-shard spread (0.049 mm) — see [results.md](results.md). Flagged to
  the paper; the likely resolution is more shards to pin the C_pool correction.
- **CsI and the rest of the grid stay native.** The dd rerun covers BGO at TBP +
  CAFOV only. Re-running CsI or the other bores on the dd source is deferred —
  the detector study is BGO-only.

## Deferred / on-request

- **Washout spatial non-uniformity / redistribution** — the only open
  isotope-washout item and the one route to a genuine, non-calibratable washout
  bias: a heterogeneous perfusion field, or cleared atoms decaying elsewhere in
  the FOV. A downstream perfusion/compartment transport model on the production
  points (not Geant4, not the range estimator). Trigger: if the edge proves
  sensitive to the clearance field. See [isotope-washout.md](isotope-washout.md),
  [washout-g4-formulation.md](washout-g4-formulation.md).
- **Scatter correction** — machinery exists (`recon_scatters.jl`,
  `scatter_profile.py`); not needed at present precision (calibration systematics
  dominate ≫ 0.1 mm). Trigger: a window slope approaching the edge gradient on
  some future configuration.
- **Composite-erfc edge model** (2–3 isotope components, offsets/widths frozen
  from per-isotope truth profiles, free amplitudes + global shift); adopt only if
  σ and rung stability improve, including at 0.1 Gy. Relevant to isotope washout —
  see [isotope-washout.md](isotope-washout.md).
- **Geometry axis** — ring length (sensitivity), open geometries (angular
  coverage); the analysis is configuration-blind, so each runs through the
  identical battery.
