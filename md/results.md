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

## Isotope washout (IW): DONE — no bias, ~1.5× σ_R cost (2026-07-10)

Washout modelled as loss (perfusion/metabolism clears emitters before they decay) reduces — for a
spatially-uniform brain — to a **per-isotope survival scalar** g_i (`latex/washout_brain.tex`
Eq. 7), so the whole loss study runs downstream on the frozen source, **zero upstream** (method +
G4/PTCrysp exchange: [`md/washout-g4-formulation.md`](washout-g4-formulation.md)). Written up as §7
of `latex/endpoint_precision.tex` (Table 6). g_i is computed from the Mizuno brain 3-exponential
(`config/washout_brain.toml`), cross-checked closed-form vs direct integration to 1e-6.

**Truth level** (`tools/washout.py`): g_i = O15 0.448, C11 0.376, N13 0.386, C10 0.525, O14 0.476
(¹⁵O the least suppressed of the abundant emitters — recorded decays are young). Kept fraction
f = 0.428. Central **ΔR₅₀^wo = +0.218 mm** (edge deepens: ¹⁵O owns the distal edge and survives
best) — **this calibrates away** like every other offset. The part that does *not*: the
washout-**parameter** systematic band, MC over the Mizuno uncertainties, is **±0.020 mm** — ~5×
below σ_R ≈ 0.11 mm, and the fast component's ~90% error barely propagates (spent for recorded
ages). So within this model washout is a **benign calibrated constant**, not a systematic limit.

**Detected level** (`drivers/washout_sigma_r.jl`, per-event thinning by w(z₀,t_decay) =
Σ_i P(i|z₀,t_decay) g_i, all events). **Bias** transfers scanner-independently and calibrates
away: ΔR₅₀^wo(t_start=0) ≈ +0.20–0.25 mm, matching truth +0.218; the shift shrinks toward zero as
the delay pre-depletes ¹⁵O (CsI +0.25 at 60 s → ~0 by 300 s). Washout adds **no bias** the
reference does not absorb.

**Precision — CORRECTED (2026-07-10).** Washout **inflates σ_R by ~1.5×; it is not free.** An
earlier claim here ("precision essentially untouched — BGO 1.11×, CsI 0.96×") was an **n=10
artifact**: the CsI 0.96× came only from dropping a legitimate shard as an "outlier" in the
ten-shard from_shards run. The robust **thinned** method (30 realizations, ±13%; dose-adaptive
0.2–0.5 Gy to keep the count-starved washed corner in the stable-fit regime — the 0.1 Gy corner was
verified to fail, σ_R blowing to 1.6 mm at 68 k events; dose-independence checked: CsI t=0 gives
**1.48× at 1 Gy vs 1.68× at 0.2 Gy, consistent**) measures the inflation curve on both arms
(scaled to 1 Gy):

| t_start [s] | 0 | 60 | 120 | 180 | 300 |
|---|---|---|---|---|---|
| CsI σ_R inflation (±13%) | 1.48 | 1.74 | 1.57 | 1.54 | 1.87 |
| BGO σ_R inflation (±13%) | 1.34 | 1.43 | 1.76 | 1.16 | 2.02 |

— **both arms ~1.4–1.6× and roughly flat**, with a rise to ~2× at the count-starved 300 s end,
**not** rising from 1.0, and **arm-consistent within the ±13% scatter** (both t=0 well above 1.0, so
BGO's old from_shards 1.11× was also an n=10 underestimate). The physics: washout removes ~57% of
counts **roughly uniformly** (all g_i ≈ 0.4–0.5, ~0.07 spread), so it is close to a uniform count
cut → σ_R × ~1/√0.43 ≈ 1.5× (counting), independent of t_start and arm. This is the ordinary count
penalty, **not** the ¹⁵O variance-drain that protects the *delayed start* — that drain removes early
¹⁵O; washout removes late decays roughly uniformly (the two were conflated in the first pass). At
1 Gy this takes σ_R from ~0.11 to ~0.16 mm. Figure: `.../comparison/figures/washout_inflation.png`
(`tools/plot_washout_thinned.py`); the superseded n=10 figures `washout_sigma_r.png` /
`washout_tstart.png` are retired.

**R35 (r 350 mm, 50 cm AFOV) — same cost in a smaller bore (2026-07-11).** Repeated the thinned
firm-up on both R35 arms (`crysp_r35_50cm_{bgo,csi}_2x0`, 10 shards each, 0.63× the ring's counts)
at the two mid-range starts, N=50 (±10%), dose-adaptive (BGO 0.3 Gy; CsI 0.5 Gy to clear the
fit-stability floor — at CsI 0.3 Gy the 180 s washed arm dropped to 204 k events and read 2.08,
settling to 1.92 at 0.5 Gy / 340 k). Inflation (washed/nominal, scaled to 1 Gy):

| t_start [s] | 120 | 180 |
|---|---|---|
| BGO σ_R inflation (±10%) | 1.74 | 1.90 |
| CsI σ_R inflation (±10%) | 1.74 | 1.92 |

— both arms tie (~1.74 at 120 s, ~1.9 at 180 s), consistent with the reference ring within the
point-to-point scatter: the ~1.7–1.9× cost tracks the count loss, not the bore. Nominal σ_R at 1 Gy
0.10–0.15 mm (just above the reference 0.11 mm, expected for the larger ring radius); washed
0.20–0.29 mm. In note §7 Table 6. Configs `config/run_parameters_r35_{bgo,csi}.toml`.

**Caveats / open:** model-form uncertainty untouched — spatial non-uniformity (the genuine
non-calibratable bias route), rabbit→human, per-species. **Bottom line (corrected):** uniform brain
washout adds **no bias** (calibrates away; parameter band ±0.02 mm ≪ σ_R) but costs **~1.5× in σ_R
(0.11 → ~0.16 mm at 1 Gy), roughly independent of start time** — the ordinary penalty for losing
~57% of the counts. The case for going upstream rests on spatial non-uniformity.

## Generation-2 σ_R study: exact washout + per-isotope, six scanners — DONE (2026-07-13)

Upstream regenerated the products as **generation v2** (`dev/reference/generation2_plan.md`):
tumour-centred phantom (`source_z_offset_mm` +25.58 → activity edge at world z ≈ **+9.13 mm**, not
−16.45), an **irradiation-end clock**, fixed acquisition scenarios (leaves
`del{120,180,300}s_ac300s_1Gy` = window `[t_del, t_del+300] s`, the delay axis is now the leaf
axis), an **isotope column** per LOR, and the stamped Mizuno **`washout_g`** per isotope. Guarded by
`generation="v2"` (`pool_shards` refuses to mix generations). The legacy off-centre `fast_1Gy`
products are superseded — gone from disk.

**Method** (`drivers/sigma_r_v2.jl`, N=100 thinned realizations, 1 Gy): pool each leaf and
reconstruct on the **recentred grid** (`img_origin` z −93.666, window `[−10.868, 24.132]` — the
legacy grid/window shifted +25.58; `characterize` gained a `z_offset_mm` kwarg, `load_run_context`
reads it from the shard, so the truth reference lands in the reconstructed frame — the must-fix,
else every fit targets an empty window). Washout is the **exact per-species Bernoulli keep** — keep
event with prob `p_dose·g_i[isotope(e)]` from the isotope column + stamped `washout_g` (the
recommended path of `latex/washout_brain.tex` §5, replacing the label-free marginalised w),
cross-checked vs a recompute from the stamped Mizuno params. Per-isotope σ_R from **pure**
isotope-column selection (no posterior-leakage lower bound). Zero fit failures anywhere.

**Six scanners:** CsI ring / R35-50 / R35-35, and BGO ring / r40-50 / r40-35. BGO carries a
cryostat → **+50 mm radius** at each AFOV class (ring r437, r40 bores r400) — the fair real-scanner
counterpart to CsI's R35 at each size-class.

**Washed σ_R [mm]** (all events, working protocol, ±7%):

| scanner | del120 | del180 | del300 |
|---|---|---|---|
| CsI ring    | 0.232 | 0.199* | 0.309 |
| CsI R35/50  | 0.235 | 0.266  | 0.395 |
| CsI R35/35  | 0.241 | 0.290  | 0.405 |
| BGO ring    | 0.133 | 0.168  | 0.232 |
| BGO r40/50  | 0.182 | 0.217  | 0.285 |
| BGO r40/35  | 0.218 | 0.236  | 0.278 |

(*ring-CsI del180 dips ~2σ below trend — an N=100 wobble; the small bores are smooth; N=200 would
tighten it, deferred.)

- **No bias, any scanner:** ΔR₅₀^wo within **±0.08 mm** — the near-uniform g_i (0.47/0.45/0.45/0.53/
  0.48) barely reweight the pre-depleted mix; calibrates away (far below the legacy +0.22 mm at t=0,
  since these windows already pre-deplete ¹⁵O).
- **Washout ~1.5×:** inflation 1.2–2.05 around the counting 1/√survival (survival ~0.39–0.46);
  tracks the count loss, not the geometry.
- **BGO wins on precision:** below CsI at every size-class (ring del120 washed 0.133 vs 0.232) — its
  2.1× counts outweigh the larger scatter fraction and the cryostat radius. Nominal σ_R BGO ring
  0.100/0.113/0.141 vs CsI 0.113/0.174/0.228.
- **Positron-range hypothesis definitively refuted (exact test):** ¹⁵O (longest β⁺ range) is *more*
  precise **per count** than ¹¹C on every scanner — ring ¹⁵O k = σ√N ≈ 236 vs ¹¹C ≈ 319. BGO
  delivers the first **clean ¹¹C** per-isotope point (0.4–1.1 M events, above the ~350 k fit floor on
  all three bores; count-starved 170–310 k on the CsI compact bores). Closes the Thread-B
  investigation ([`md/sigma-r-investigation.md`](sigma-r-investigation.md)).

In the note as **§8** (CsI + BGO v2 subsections + the BGO-vs-CsI comparison; Figs 9–13, Tables
7–10). Outputs `out/…/<scanner>/<label>/washout_v2/{sigma_r_washout_v2,sigma_r_per_isotope_v2}.toml`;
figures via `tools/plot_sigma_r_v2.py`, `plot_washed_v2_scanners.py`, `plot_washed_bgo_vs_csi.py`.
The two flagship reference configs are `config/run_parameters_{csi_v2,ring_bgo_v2}.toml`.

## Data on disk

Products at `PtCryspProds/uniform_headep_sobp_1e8/` are now **generation v2** (self-describing:
`generation`, `t_decay_zero=irradiation_end`, `center_on=tumour`, isotope column, `washout_g`;
schema `dev/reference/SCHEMA.md`). Per scanner+crystal, three 10-shard scenario leaves
`del{120,180,300}s_ac300s_1Gy` — for the six arms **ring / R35-50(CsI) / R35-35(CsI) / r40-50(BGO) /
r40-35(BGO)** in both flagship crystals (CsI ring, BGO ring) plus the four small bores. The legacy
off-centre `fast_1Gy` masters are **superseded and removed** (recoverable only via the committed
legacy `out/` results + git history). Shared `truth/` bundle unchanged
(`activity_profile_fast.csv` per-isotope columns `O15,C11,N13,C10,O14,total`; `depth_dose.csv`) —
still in the native frame, so the v2 chain shifts it by `source_z_offset_mm`. Config
`[configuration]` is parked on the ring BGO v2 arm.
