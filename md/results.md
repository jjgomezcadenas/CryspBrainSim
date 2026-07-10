# Science results

## The finding in one paragraph

Both representative scanners locate the distal edge with **σ_R ≈ 0.11 mm per 1 Gy run** at the
working protocol (all events, no corrections); precision scales as 1/√dose (0.16–0.23 mm at
0.1 Gy → the exploratory-dose case) and **survives a realistic in-room acquisition start** (still
~0.11 mm at t_start = 180 s). The activity–dose offset is a per-scanner calibration constant
measured to 0.01–0.04 mm/term; scatters bias it ≤ 0.04 mm (no correction needed at this level —
other calibration terms dominate). Written up in `latex/endpoint_precision.tex` (8 pages, compiles
clean via `tools/latex_compile.py endpoint_precision`).

## Endpoint study, part (a): distal-edge estimation — DONE (2026-07-08)

The study splits: (a) distal-edge estimation, (b) scanner comparison.

**Settled protocol:** whole-plane profiles (no ROI — the 13 mm disc clips the depth-widening halo,
shifting R_p proximally 2.4–3.0 mm), erfc edge fit with FREE baseline (forcing b = 0 shifts R50 by
0.4–0.6 mm at our statistics; the fitted b is a small negative shape-slack, −1.5…−8% of amplitude,
not a background), **R50 (= fitted z0) as THE observable**. Its ~−11 mm offset to the dose R50 is a
calibration constant fixed by the reference simulation; the measurement delivers variations against
that anchor.

- **Model cross-check** (fit lab `--model erfc|sigmoid|both`, `--no-baseline`): the logistic sigmoid
  (Zapien-Campos Eq. 3 + baseline) gives identical R50 information — constant +0.10 mm offset,
  identical spreads, rung shifts equal to 0.02 mm — and slightly better χ²; erfc stays primary,
  sigmoid is the built-in cross-check.
- **R_p** (tangent endpoint) demoted to a qualitative accuracy statement: swings ~1.2 mm across
  rungs (the shard-0 "blur-stable" reading was a fluctuation; ten-shard mean shift truth→recon is
  +0.32 mm), 4× the shard spread of R50. **R_x** (1%) is a tail diagnostic only (model-dependent by
  ~7 mm).
- **Ten-shard ladder** (`tools/ten_shards.py` → `ten_shards/results.toml`): Δ_R50 = R50(act) −
  R50(dose), erfc — truth activity −10.743; origins (acceptance only) −10.988 ± 0.010 [std 0.031];
  recon(trues) −11.216 ± 0.018 [0.057]; recon(all events) −11.194 ± 0.022 [0.071]. Budget:
  acceptance −0.25, reconstruction −0.23, scatters +0.02 mm.
- **Dose sweep** (`drivers/ten_shards_dose.jl`, `tools/ten_shards.py --dose-sweep`): Δ_R50
  dose-invariant and σ_R ∝ 1/√dose exactly: 0.057/0.072/0.117/0.177 mm at 1.0/0.5/0.2/0.1 Gy.
  **Test-dose statement: a single 0.1 Gy acquisition locates the edge to σ_R ≈ 0.18 mm, bias
  ≤ 0.04 mm** (trues, fast window, closed ring).

Literature anchor: Zapien-Campos et al., Med Phys 2025 (`papers/`, untracked) — logistic fit,
PAR = R50, whole-plane, offset stable to 0.4–0.5 mm; their range-shift transfer test is their
problem (IMPT spots), not ours (fixed field, dose axis). Note: the numbers in this block are from
the **frozen** `crysp_ring_1m/bgo` master; the two-scanner campaign supersedes them per scanner.

## Two-scanner campaign (part b, first round): DONE (2026-07-10)

Upstream published two ten-shard masters on the same 2X0 ring from identical annihilation sets
(`dev/reference/PRODUCTS.md`/`SCHEMA.md`; production note PTCryspMC.jl `latex/scanner_prods.tex`):

- **BGO 195 K** — 15.3 M LORs/shard, 72.3/27.4/0.3% true/scatter/random, eres 15%, 413 keV, τ 5 ns.
- **cryogenic CsI** — 6.1 M, 86.6/13.3/0.1%, eres 6%, 472 keV, τ 1.5 ns.

Both σ_xyz 1.486 mm (3.5 mm FWHM), carry `t_decay_s`. The old `crysp_ring_1m/bgo` master is
**frozen** (methods reference; do not pool/compare).

Results (whole-plane erfc, free b), ladder Δ_R50 (origins / recon-trues / all-events, std in
brackets):

| | BGO 195 K | CsI |
|---|---|---|
| origins | −10.938 [0.023] | −10.962 [0.039] |
| recon (trues) | −11.165 [0.065] | −11.232 [0.084] |
| recon (all events) | −11.144 [0.113] | −11.192 [0.114] |

Scatter shifts +0.021 / +0.040 mm (no bias correction needed either arm); dose sweeps flat,
σ_R ∝ 1/√dose. **Headline: at the uncorrected all-events protocol the arms tie — σ_R(1 Gy) = 0.113
(BGO) vs 0.114 mm (CsI); BGO's 2.1× statistics cancels its 2.1× scatter fraction. On trues they
differ (0.065 vs 0.084 mm) — a scatter correction is the tiebreaker.** Paired shards buy nothing
for R50 (detection sampling dominates; paired std = quadrature sum). Upstream nit: the BGO arm's
`scanner_geometry.json` says material BGO_77K while shards say BGO_195K (stale key, harmless —
constants come from shard attrs).

## Delayed acquisition-start study: DONE (2026-07-10)

`shard_t_decay` (src/products.jl, asserts zero dropped LORs for alignment),
`drivers/ten_shards_tstart.jl` (t_decay ≥ t_start cut, trues + all events),
`tools/ten_shards.py --t-start`, `tools/plot_tstart.py` (cross-scanner figure →
`out/<scenario>/closed/comparison/figures/tstart_r50.png`).

Results (all events, Δ_R50 mean ± σ): kept 77/61/48% at t_start 60/120/180 s (identical both
scanners, cut acts on shared decay times). **σ_R stays ≈ 0.11 mm through 180 s on both scanners —
at/below the counting prediction; the zero-delay precision survives the realistic in-room start.**
Calibration walks −0.44/−0.99/−1.6 mm (¹⁵O drains from the mix; its lower production threshold =
deeper edge) → start-time sensitivity 8–10 μm/s, 0.01 mm at 1 s timing. Edge-sharpening hypothesis
refuted: w grows 10.8→11.4 mm (¹³N foot gains weight; ¹¹C positron-range gain invisible under the
intrinsic width). In the note as §6 + Table 5 + Fig. 7.

**Scatter-correction stance settled:** NOT needed at this level (calibration systematics dominate
≫ 0.1 mm); CsI's point is parity with less scatter sensitivity.

## Data on disk

Three ten-shard masters under `PtCryspProds/uniform_headep_sobp_1e8/`, all **physical-decay-only**
(no isotope washout):

- `crysp_ring_1m/bgo/fast_1Gy/` — the frozen reference (174.3 M LORs pooled; **do not** pool/compare
  with the new arms).
- `crysp_ring_1m_bgo_2x0/bgo_195k/fast_1Gy/` — BGO 195 K, 153 M pooled.
- `crysp_ring_1m_csi_2x0/csi/fast_1Gy/` — CsI, 61 M pooled.

Plus the shared `truth/` bundle (`activity_profile_fast.csv` has per-isotope depth columns
`O15,C11,N13,C10,O14,total`; `depth_dose.csv`). Each new-arm shard carries `t_decay_s`. Config
`[configuration]` is parked on the BGO 195 K arm.
