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

Results (all events, Δ_R50 mean ± σ): kept 77/61/48/32% at t_start 60/120/180/300 s (identical
both scanners, cut acts on shared decay times). **σ_R stays ≈ 0.11 mm through 180 s on both
scanners — at/below the counting prediction; the zero-delay precision survives the realistic
in-room start.** The ¹⁵O drain removes a variance source with the counts (CsI transiently
over-compensates: 0.057 mm at 60 s, a 2σ dip at n = 10); at 300 s the compensation is exhausted
on CsI (0.233 mm ≈ counting) while BGO still reads 0.104 mm. Calibration walks
−0.44/−0.99/−1.6/−2.9 mm (¹⁵O drains from the mix; its lower production threshold = deeper edge)
→ start-time sensitivity 8–11 μm/s, 0.01 mm at 1 s timing. Edge-sharpening hypothesis refuted:
w grows 10.8→11.4 mm (¹³N foot gains weight; ¹¹C positron-range gain invisible under the
intrinsic width). In the note as §6 + Table 5 + Fig. 7 (story figures 4/6 centered, 5 in shared
log fractional units — the offsets live in the tables, the spreads in the figures).

**Scatter-correction stance settled:** NOT needed at this level (calibration systematics dominate
≫ 0.1 mm); CsI's point is parity with less scatter sensitivity.

## Single-shard geometries — CHS and R35, first look: DONE (2026-07-10)

Upstream published one shard (realization 0, same annihilation set as the ring shards) per arm of
two new geometries, both with the ring's 2X0 walls and detector constants:

- **CHS** (compact head scanner): r_inner 200 mm, 25 × 7 blocks, axial 358 mm. Per 1 Gy shard:
  BGO 195 K 10.14 M LORs, CsI 3.94 M — 0.66×/0.65× the ring's counts; true/scatter/random mix
  identical per crystal (72.5/27.3/0.2 and 86.6/13.3/0.1%).
- **R35** (half-metre ring): r_inner 350 mm, 43 × 10 blocks, axial 512 mm. Per 1 Gy shard:
  BGO 9.67 M, CsI 3.83 M; mix 73.9/25.9/0.2 and 87.1/12.9/0.1%.

**Protocol unification first:** the Julia chain (`reconstruct_endpoint`, `one_shard.jl`, the
`sigma_r_*` drivers) still carried the rung-5 disc ROI while the settled protocol (part a) is
whole-plane; the ROI is now retired everywhere (`[roi]` in `run_parameters.toml` has no radius).
Gate: the Julia whole-plane fit reproduces the python campaign fit on the same stored image to
1e-6 mm, both selections. `sigma_r_sweep_dose.jl` gained `--all-events` (paired thin seeds;
`sweep_all.toml`).

With one shard per arm, σ_R at 1 Gy is not measurable (thinning is degenerate at fraction 1);
the deliverables are σ_R measured at exploratory doses and the 1/√dose extrapolation, validated
on the ring arms where the ten-shard 1 Gy measurement exists. σ_R (whole-plane erfc fit, 50
thinned realizations, spread known to ±10%), per geometry at matched dose (ring / CHS / R35):

| σ_R [mm] | 0.2 Gy | 0.1 Gy | 0.05 Gy | 1 Gy (extrap) |
|---|---|---|---|---|
| BGO trues | 0.129 / 0.122 / 0.124 | 0.184 / 0.177 / 0.147 | 0.221 / 0.223 / 0.260 | 0.054 / 0.053 / 0.052 |
| BGO all | 0.158 / 0.194 / 0.186 | 0.232 / 0.303 / 0.240 | 0.284 / 0.433 / 0.361 | 0.069 / 0.093 / 0.080 |
| CsI trues | 0.192 / 0.183 / 0.188 | 0.256 / 0.241 / 0.210 | 0.332 / 0.312 / 0.305 | 0.080 / 0.075 / 0.071 |
| CsI all | 0.231 / 0.213 / 0.230 | 0.313 / 0.313 / 0.292 | 0.469 / 0.400 / 0.549 | 0.102 / 0.094 / 0.103 |

**Headline: on trues all three geometries are indistinguishable per dose (1 Gy extrapolations
0.052–0.054 mm BGO, 0.071–0.080 mm CsI) although the compact geometries collect 0.63–0.66× the
counts — per detected true they are ~1.25× more precise, cancelling the acceptance loss. On the
all-events working protocol the BGO penalty is monotonic in bore radius (ring 0.069 < R35 0.080 <
CHS 0.093 mm extrapolated; scatter fractions nearly equal, so the compact ring is more sensitive
to the background's shape, not its amount), while CsI is geometry-flat (0.094–0.103 mm).**
Extrapolation validation on the ring: trues agree with the ten-shard 1 Gy measurement (0.054 vs
0.065 ± 0.015; 0.080 vs 0.084 ± 0.020); all-events agree on CsI (0.102 vs 0.114 ± 0.027) and read
1.6σ low on BGO (0.069 vs 0.113 ± 0.027) — the thinned method emulates counting statistics only,
so quote the single-shard 1 Gy all-events extrapolations with that caveat.

Calibration anchors at 1 Gy (single shard, whole-plane; Δ_R50 vs the dose fit −3.580 mm): CHS
trues/all −11.42/−11.30 (BGO), −11.59/−11.61 (CsI); R35 −11.10/−11.03 (BGO), −11.09/−10.82
(CsI); ring ten-shard means −11.17/−11.14 and −11.23/−11.19. Per-scanner calibration constants,
spread 0.4 mm across geometries — as expected, fixed by each scanner's reference simulation.
Figures: `out/uniform_headep_sobp_1e8/closed/comparison/figures/chs_sigma_r.png` (+ `_bgo`,
`_csi` single-crystal versions; `tools/plot_chs_sigma_r.py [--crystal]`). A measured σ_R(1 Gy)
and a shard-spread on the anchors need the remaining nine shards per arm from upstream.

## Isotope washout (IW), first result: DONE at t_start = 0 (2026-07-10)

Washout modelled as loss (perfusion/metabolism clears emitters before they decay) reduces — for a
spatially-uniform brain — to a **per-isotope survival scalar** g_i (`latex/washout_brain.tex`
Eq. 7), so the whole loss study runs downstream on the frozen source, **zero upstream** (method +
G4/PTCrysp exchange: [`md/washout-g4-formulation.md`](washout-g4-formulation.md)). g_i is computed
from the Mizuno brain 3-exponential (`config/washout_brain.toml`), cross-checked closed-form vs
direct integration to 1e-6.

**Truth level** (`tools/washout.py`): g_i = O15 0.448, C11 0.376, N13 0.386, C10 0.525, O14 0.476
(¹⁵O the least suppressed of the abundant emitters — recorded decays are young). Kept fraction
f = 0.428. Central **ΔR₅₀^wo = +0.218 mm** (edge deepens: ¹⁵O owns the distal edge and survives
best) — **this calibrates away** like every other offset. The part that does *not*: the
washout-**parameter** systematic band, MC over the Mizuno uncertainties, is **±0.020 mm** — ~5×
below σ_R ≈ 0.11 mm, and the fast component's ~90% error barely propagates (spent for recorded
ages). So within this model washout is a **benign calibrated constant**, not a systematic limit.

**Detected level** (`drivers/washout_sigma_r.jl`, per-event thinning by w(z₀,t_decay) =
Σ_i P(i|z₀,t_decay) g_i, all events, ten shards). Bias transfers scanner-independently:
BGO +0.200, CsI +0.209 mm (vs truth +0.218). **Precision essentially untouched despite losing 57%
of counts** — σ_R BGO 0.113→0.126 (1.11×), CsI 0.114→0.110 (0.96×, dropping one n=10 fit outlier;
all-10 gives 0.237/2.08×), both far below the naïve 1/√f = 1.53× — the same ¹⁵O variance-drain as
the delayed-start study: washout removes the noisiest counts. Figure:
`out/uniform_headep_sobp_1e8/closed/comparison/figures/washout_sigma_r.png`
(`tools/plot_washout.py`).

**Compounded with the in-room start** (`washout_sigma_r.jl` t_start sweep, both arms, g_i
recomputed per shifted window [120+t_start, 1320]): the two effects stack on ¹⁵O, and the picture
changes from the t_start = 0 limit:

| t_start [s] | 60 | 120 | 180 | 300 |
|---|---|---|---|---|
| kept (washout × delay) | 0.31 | 0.23 | 0.18 | 0.11 |
| ΔR₅₀^wo BGO / CsI [mm] | +0.22 / +0.30 | +0.18 / +0.20 | +0.19 / +0.06 | +0.10 / −0.00 |
| σ_R washed BGO / CsI | 0.16 / 0.13 | 0.12 / 0.21 | 0.16 / 0.17 | 0.21 / 0.31 |
| σ_R nominal BGO / CsI | 0.10 / 0.06 | 0.12 / 0.09 | 0.11 / 0.11 | 0.10 / 0.23 |

Two effects: **(a)** the washout shift shrinks toward zero as the delay pre-depletes ¹⁵O (BGO
+0.22→+0.10, CsI +0.30→−0.00) — still calibrating away, just smaller; **(b)** the precision that
was free at t_start = 0 **is no longer free** — washed σ_R climbs above nominal (~1.3–2× at the
operational 180–300 s), because the ¹⁵O that supplied the variance-drain is gone by then, so
washout's ~57%-per-step count loss now bites like ordinary counting. Figure:
`out/uniform_headep_sobp_1e8/closed/comparison/figures/washout_tstart.png`. (σ_R here is n=10,
±24%, with visible wiggle — the 50-realization thinned method is needed to state the inflation
factors cleanly.)

**Caveats / open:** **model-form** uncertainty is untouched — spatial non-uniformity (the genuine
non-calibratable bias route), rabbit→human, per-species; the parameter systematic (±0.020 mm at
t_start = 0) is negligible but the σ_R firm-up is pending. **Bottom line:** with uniform brain
washout the bias always calibrates away and its parameter band is negligible, so IW does not bias
range verification; its one real cost is **precision, and only when compounded with a realistic
delayed start** — at t_start = 0 it is free, at 180–300 s it inflates σ_R ~1.3–2×. The case for
going upstream still rests on spatial non-uniformity.

## Data on disk

Under `PtCryspProds/uniform_headep_sobp_1e8/`, all **physical-decay-only** (no isotope washout):

- `crysp_ring_1m/bgo/fast_1Gy/` — the frozen reference (174.3 M LORs pooled; **do not** pool/compare
  with the new arms), ten shards.
- `crysp_ring_1m_bgo_2x0/bgo_195k/fast_1Gy/` — BGO 195 K, ten shards, 153 M pooled.
- `crysp_ring_1m_csi_2x0/csi/fast_1Gy/` — CsI, ten shards, 61 M pooled.
- `crysp_chs_bgo_2x0/bgo_195k/fast_1Gy/` — CHS BGO 195 K, **shard 000 only**, 10.1 M LORs.
- `crysp_chs_csi_2x0/csi/fast_1Gy/` — CHS CsI, **shard 000 only**, 3.9 M LORs.
- `crysp_r35_50cm_bgo_2x0/bgo_195k/fast_1Gy/` — R35 BGO 195 K, **shard 000 only**, 9.7 M LORs.
- `crysp_r35_50cm_csi_2x0/csi/fast_1Gy/` — R35 CsI, **shard 000 only**, 3.8 M LORs.

Plus the shared `truth/` bundle (`activity_profile_fast.csv` has per-isotope depth columns
`O15,C11,N13,C10,O14,total`; `depth_dose.csv`). Every 2x0 and chs shard carries `t_decay_s`.
Config `[configuration]` is parked on the ring BGO 195 K arm.
