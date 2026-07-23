# Pending (smaller items, independent of the isotope-washout study)

- **BIG: bounded fit supersedes the old σ_R numbers.** The bounded `fit_endpoint`
  (results.md → "Statistical-procedure study") tightens σ_R ~25–40% vs the old unbounded fit, worse
  at low counts. So **`endpoint_precision.tex` §8 (and everything quoting old `sigma_r_v2` σ_R)** is
  systematically high and should be regenerated with the bounded fit / `statistical_procedure_jobs.jl`,
  or explicitly caveated. `cbs.tex` Results is already reworked to bounded-fit numbers; the rest of the
  notes are not.
- **Merge branch `paper/statistical-procedure`** to main when ready (currently 6+ commits ahead;
  cbs.tex rework + `plot_statproc_delay_csi.py` + last-point results still uncommitted on it).
- **`cbs.tex` Fig 9 caption** (the "Statistical procedure" section, another instance's) says
  acquisition window `[120,300] s`, but the statistical-procedure data is the full `[120,420] s`
  del120 window — likely a stale number to correct.
- ~~Ring-CsI del180 N=200 firm-up~~ **DONE** — N=200 pulled it to 0.221 (was 0.199 fluctuation);
  superseded anyway by the bounded-fit statistical-procedure re-measure.
- **CHS + R35: remaining nine shards per arm** — obsolete for v2 (the v2 R35/40 arms already ship
  10 shards); kept only for the legacy single-shard CHS/R35 extrapolations of §1–7.
- **BGO all-events bore-radius sensitivity** — at matched dose the BGO all-events penalty grows
  monotonically as the bore shrinks (ring < R35 < CHS) while trues match everywhere, with near-equal
  event mix; understand the scatter background shape in compact geometries. Related: the ring BGO
  all-events 1/√dose extrapolation sits 1.6σ below the ten-shard 1 Gy measurement (counting-only
  thinning caveat).

- **`latex/cbs.tex`** (the living draft, separate from `endpoint_precision.tex`): the **Results
  section is reworked** (2026-07-23) around two bounded-fit plots + one summary table (results.md);
  the old σ_R tables are gone. Remaining smaller items if wanted: eq:sigmaR wording, cite
  Zapien-Campos. NB: `cbs.tex` carries other instances' uncommitted edits — commit only on explicit
  instruction.
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
