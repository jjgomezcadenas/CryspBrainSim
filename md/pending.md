# Pending (smaller items, independent of the isotope-washout study)

- **`latex/cbs.tex`** (the living draft, separate from `endpoint_precision.tex`): revert eq:sigmaR
  from R_p back to z0/R50 (keep the R_p accuracy paragraph), cite Zapien-Campos, fold in the
  two-scanner numbers. NB: `cbs.tex` carries the user's own uncommitted edits — commit only on
  explicit instruction.
- **Composite-erfc edge model** (2–3 isotope components, offsets/widths frozen from per-isotope
  truth profiles, free amplitudes + global shift); adopt only if σ and rung stability improve
  incl. at 0.1 Gy. Directly relevant to isotope washout — see [isotope-washout.md](isotope-washout.md).

## Deferred / on-request

- **Scatter correction** — machinery exists (`recon_scatters.jl`, `scatter_profile.py`); not needed
  at present precision (calibration systematics dominate ≫ 0.1 mm). The trigger for revisiting: a
  window slope approaching the edge gradient on some future configuration.
- **Geometry axis** (part b continued) — ring length (sensitivity), open geometries (angular
  coverage); the analysis is configuration-blind, so each runs through the identical battery.
