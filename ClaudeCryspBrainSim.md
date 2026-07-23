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

The build is complete; the physics studies are done through the acquisition-start axis, the
isotope-washout loss study, and the generation-2 exact-washout + per-isotope study across six
scanners. Details live in `md/` (see "Read first") — this is the one-paragraph state.

**Where we are:** both representative scanners (BGO 195 K, cryogenic CsI, same 2X0 ring) measure
the distal edge with **σ_R ≈ 0.11 mm per 1 Gy run** at the working protocol; precision scales as
1/√dose (→ 0.16–0.23 mm at the 0.1 Gy exploratory dose) and survives a realistic in-room
acquisition start (~0.11 mm at t_start = 180 s). Written up in `latex/endpoint_precision.tex`
(the two-scanner note, compiles clean). Two single-shard geometries are measured against the
ring — the **compact head scanner (CHS**, r 200 mm) and **R35** (r 350 mm, AFOV 512 mm): on
trues all three tie per dose despite the compact arms' 0.65× counts; on all-events the BGO
penalty grows as the bore shrinks (ring < R35 < CHS) while CsI stays flat (see results.md). The
whole-plane protocol now lives in the Julia chain too (`[roi]` carries no radius). Full numbers
+ method: [`md/results.md`](md/results.md); the toolchain:
[`md/infrastructure.md`](md/infrastructure.md).

**Isotope washout (IW) — loss study DONE (both arms):** productions are physical-decay-only; IW is
the missing simulation-to-patient physics. Washout-as-loss is a per-isotope survival scalar g_i, so
the study runs fully downstream (zero upstream, no Geant4; detection is isotope-blind, the
detected σ_R comes from per-event thinning by w(z₀,t_decay) — only columns already in the shards).
Result: **no bias** — the edge shift (+0.22 mm) calibrates away, parameter systematic ±0.02 mm
(≪ σ_R) — but a **~1.5× σ_R cost** (0.11 → ~0.16 mm at 1 Gy, roughly flat vs t_start, both arms),
the ordinary penalty for a ~57% near-uniform count loss. Measured with the thinned method; the
earlier "σ_R survives / free at t=0" was an n=10 artifact, now corrected. In the note as §7
(`latex/endpoint_precision.tex`). **Only open IW item: spatial non-uniformity** — the one genuine
bias route, needing a downstream perfusion-transport model (not pursued). Detail:
[`md/isotope-washout.md`](md/isotope-washout.md),
[`md/washout-g4-formulation.md`](md/washout-g4-formulation.md), `latex/washout_brain.tex`.

**Generation-2 (v2) σ_R study — DONE (six scanners):** upstream regenerated the products as v2
(tumour-centred, irradiation-end clock, per-LOR isotope column, stamped Mizuno `washout_g`, fixed
`del{120,180,300}` scenario leaves). Consumed downstream by `drivers/sigma_r_v2.jl` (N=100, 1 Gy):
washout as the **exact per-species g_i keep** (isotope label present) and per-isotope σ_R from
**pure** selection. Across CsI (ring/R35-50/R35-35) and BGO (ring/r40-50/r40-35, +cryostat): **no
bias** (ΔR₅₀ ±0.08 mm), washout **~1.5×** cost (tracks counts, not bore), **BGO more precise** at
every size-class (2.1× counts), and the **positron-range hypothesis definitively refuted** — ¹⁵O is
more precise per count than ¹¹C, with BGO giving the first clean ¹¹C point. In the note as §8. Prep:
`generation` guard + `shard_isotope` (products.jl), tumour-centring `z_offset` into `characterize`,
recentred grid + sensitivities, flagship configs `run_parameters_{csi_v2,ring_bgo_v2}.toml`. Full
numbers: [`md/results.md`](md/results.md); toolchain: [`md/infrastructure.md`](md/infrastructure.md).

**Statistical-procedure study + bounded fit — DONE (branch `paper/statistical-procedure`):** a new
restartable driver (`drivers/statistical_procedure_jobs.jl`, shard/ensemble/combine, `--tend` window
sub-cut, finite-pool correction) plus a **bounded** `fit_endpoint` re-measure σ_R rigorously. The
headline: the old unbounded fit **inflated σ_R ~25–40%** (worse at low counts), so the old
`sigma_r_v2` / `endpoint_precision.tex` tables are systematically high. Bounded-fit washed σ_R
(del120, 1 Gy): BGO ring **0.10**, CsI ring **0.15**; grid across TBP/LAFOV/CAFOV and a CsI
delay×duration series (all < 0.35 mm). Method cross-checked: direct 10-shard σ_R ≈ thinned
corrected; washout inflation 1.44 = counting. `cbs.tex` Results reworked to two bounded-fit plots +
one summary table. Detail: [`md/results.md`](md/results.md).

**Smaller pending:** **regenerate `endpoint_precision.tex` §8 with the bounded fit** (old numbers
~25–40% high), merge the `paper/statistical-procedure` branch, composite-erfc model —
[`md/pending.md`](md/pending.md).
