# CryspBrainSim — repository context

Orientation for any Claude Code session on this repo. (Imported by `CLAUDE.md`,
which holds the general working rules.)

## Purpose

Take the list-mode LORs a PET scanner would record for a proton SOBP field, reconstruct the β⁺
activity, extract the distal range endpoint R, and report its statistical precision σ_R vs dose —
one curve per scanner geometry. This is the analysis end of the chain
`ptcrysp-scenarios → PTCryspMC.jl → PtCryspProds/ → here`; the reconstruction engine comes from
`RecoCryspTools` as a pinned dependency.

## Read first

- **`md/`** — the project state, kept out of this file so it stays lean:
  [`infrastructure.md`](md/infrastructure.md) (package, drivers, tools),
  [`results.md`](md/results.md) (the science done + numbers + data on disk),
  [`sigma-r-investigation.md`](md/sigma-r-investigation.md) (the four-geometry σ_R +
  positron-range investigation, **closed** by the v2 exact per-isotope test),
  [`isotope-washout.md`](md/isotope-washout.md) (the washout loss study, done),
  [`pending.md`](md/pending.md) (smaller items). Keep these current; this file points to them.
- **`dev/PLAN.md`** — the build plan: structure, dependencies, the consume-vs-write code inventory,
  the frozen run parameters, the validation ladder, the deferred register, and the build order.
- **`dev/reference/`** — vendored snapshots of the upstream contracts (products tree, data-generation
  strategy, σ_R recipe, LOR schema, RecoCrysp usage), with provenance in `dev/reference/README.md`.

## Status

The build is complete and the range-precision study is finished on the **data-driven source with
the BGO detector**. Details live in `md/` (see "Read first") — this is the one-paragraph state.

**Where we are — data-driven BGO σ_R, two scanner sizes, five protocols (DONE, branch `BGOv2`):**
the analysis runs on scenario `uniform_headep_sobp_1e8_dd`, whose emitters are sampled from the
**nominal fitted production cross sections** (a fit to EXFOR data) rather than Geant4's internal
model (upstream ptcryspg4 `workshop/xsections_phases.md`; ×1.32 ¹¹C, ×1.49 ¹³N vs native; the dd
fitted edge sits −0.29 mm from native). Two BGO scanners are measured — **TBP** (`crysp_ring_1m`,
1 m AFOV, the reference) and **CAFOV** (`crysp_r40_35cm`, 35 cm) — across five acquisition protocols
(delay/scan s): d120s300 d180s300 d300s300 d120s120 d180s120. Washed σ_R at 1 Gy, finite-pool-corrected,
N=100:

| protocol | TBP | CAFOV |
|---|---|---|
| d120s300 | 0.128 | 0.142 |
| d180s300 | 0.142 | 0.174 |
| d300s300 | 0.198 | 0.249 |
| d120s120 | 0.162 | 0.184 |
| d180s120 | 0.196 | 0.224 |

All ≤ 0.25 mm; TBP leads CAFOV by ~10–25%. The reference point (TBP, d120s300) is σ_R = 0.128 mm
washed (0.084 nominal). Numbers + figures: [`md/results.md`](md/results.md); auto-generated tables
`latex/sigma_r_bgo_table.tex` and `latex/statistical_procedure_results.tex`, figures `latex/figs/`;
toolchain [`md/infrastructure.md`](md/infrastructure.md).

**Systematics (ptcryspg4 side, folded into the paper):** cross-section fit u_xs = ±0.13 mm;
transport / physics-list envelope u_transport = ±0.10 mm; tumour-composition +0.04 mm (soft tissue)
/ +0.13 mm (water) at d120s300. These are the calibration systematic alongside the statistical σ_R.

**Finite-pool validation (open thread):** the thinning correction C_pool = 1/√(1−q) should reproduce
a true-Poisson spread; the numerical check on the dd TBP reference (nominal N=100 thinned, corrected
0.084 mm, vs the direct ten-shard spread 0.049 mm) sits ~2.7σ high — flagged for the paper, likely
the ten-shard estimate being too coarse to pin the correction.

**How we got here (method lineage; absolute numbers superseded by the dd run above):** earlier
campaigns on the Geant4-internal source ("generation v2": tumour-centred, per-LOR isotope column,
Mizuno `washout_g`) built and validated the machinery — the whole-plane erfc endpoint (R50), the
delayed-start axis, the isotope-washout loss study (no bias; ~1.5× σ_R cost for the ~57% count
loss), the bounded-fit re-measure (the old unbounded fit ran ~25–40% high), and finite-pool-corrected
thinning. Those studies also ran a **cryogenic CsI** arm and several bore sizes; the detector study
is now **BGO-only**, so the CsI and multi-bore numbers are history. Record:
[`md/results.md`](md/results.md), [`md/isotope-washout.md`](md/isotope-washout.md),
[`md/sigma-r-investigation.md`](md/sigma-r-investigation.md).

**The paper** (`cbs.tex`) lives in a separate repository (`~/Papers/CryspBrain`); there is no local
copy in this repo.
