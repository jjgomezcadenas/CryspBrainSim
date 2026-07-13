# σ_R across geometries + the positron-range variance question (CLOSED, 2026-07-13)

Investigation opened after the washout loss study. Two linked threads: (a) a four-configuration
σ_R comparison, and (b) whether **positron range** sets a per-isotope σ_R floor. This file is the
compaction anchor; numbers + tools below.

**RESOLVED by the generation-2 exact test (see [`md/results.md`](results.md) "Generation-2 σ_R
study").** The soft-posterior per-isotope selector here gave only a lower bound (isotope leakage);
v2 carries the isotope label, so `drivers/sigma_r_v2.jl` does a **pure** single-species selection.
Across all six v2 scanners ¹⁵O (longest β⁺ range) is **more precise per count** than ¹¹C
(ring ¹⁵O k = σ√N ≈ 236 vs ¹¹C ≈ 319) — the opposite of a range penalty — with BGO giving the first
clean ¹¹C point (above the fit floor). The positron-range σ_R floor is **definitively refuted**; the
Thread-A "protection"/inversion wobbles are N-limited counting scatter, not an isotope effect.

## Thread A — four-geometry washout σ_R comparison

Extended the washout thinned firm-up (`drivers/washout_sigma_r.jl --thinned`, t_start 120/180 s,
N=50, dose-adaptive, scaled to 1 Gy) to four scanners, **same radius family**:

| config | r (mm) | AFOV (mm) | dose (BGO/CsI) |
|---|---|---|---|
| Ring 1m | ~387 | 1024 | 0.3 / 0.3 |
| R35/50 | 350 | 500 | 0.3 / **1.0** (re-run) |
| R35/35 | 350 | 350 | 0.5 / **1.0** |
| CHS | 200 | (compact) | 0.3 / 0.5 |

Figure: `out/…/comparison/figures/sigma_r_configs.png` (`tools/plot_sigma_r_configs.py`, washed
only, two start times, extensible `CONFIGS` list). R35/35 needed its 1e9 sensitivity built first
(`tools/make_sensitivity.jl`); the others had it.

**The puzzle (the user's "nothing matters"):** washed σ_R across configs is ~flat **0.22–0.29 mm**
though washed counts span **350k–885k (~2×)**. Pure counting would spread σ_R by √2 ≈ 1.4×; it
doesn't — and at 180 s it even runs backwards (R35/50 has the most counts and the highest σ_R).
Yet **within a config σ_R IS counting-limited** (dose sweeps scale as 1/√dose; R35/50 CsI washed
gives ~0.29 at 0.3/0.5/1.0 Gy). So the per-geometry constant `k = σ_R·√N` differs by config →
σ_R across geometries is set by something other than raw counts.

## What is ruled out

- **Data artifact** — R35/35 CsI's 10 shards are internally identical and statistically identical
  to R35/50 CsI (same t_decay median ~171 s, frac(t≥180)=0.483). No bad shard.
- **Reconstruction iteration count** — `drivers/washout_niter_scan.jl` (records R50 at every MLEM
  checkpoint from one pass) + `tools/plot_sigma_r_vs_niter.py` → `sigma_r_vs_niter.png`. On R35/50
  and R35/35 CsI: σ_R **rises monotonically from niter=10 and plateaus by ~50** (no interior
  minimum — the user's parabola idea is not what happens), the mean R50 (bias) has **converged by
  ~40–50**, and the config ordering is identical at every niter. So the frozen **niter=50 is
  well-chosen** and iterations do NOT cause the cross-config differences.
- **The apparent R35/50 > R35/35 "inversion" at 180 s** — was N=50 noise. At **N=200**, R35/35 CsI
  washed = 0.249 (120 s) / 0.243 (180 s) — tied within ±5%; the N=50 value 0.222 was a downward
  fluctuation.

## What survives (real, ~3σ at N=200)

R35/35 CsI, 1 Gy: nominal σ_R **rises** 0.134→0.177 (32%, ordinary count loss), washed stays
**flat** 0.249→0.243 → washout **inflation falls 1.85→1.37**. So the delayed start / washout
**protects** σ_R against the count loss. Loose end: R35/50 CsI inflation was flat 1.56/1.56 but
only at N=50 (±14%), so the protection's config-specificity is unconfirmed (needs R35/50 CsI N=200).

## Thread B — the candidate physics: positron-range variance (the current test)

**Hypothesis (user):** ¹⁵O has the highest β⁺ endpoint energy → longest positron range (~2.5 mm vs
¹¹C ~1.1 mm) → largest per-LOR position noise → a σ_R penalty *at fixed counts*. As ¹⁵O drains
(later start / washout), that variance term shrinks → σ_R protected. **Only σ_R matters; the
absolute R (production depth) is irrelevant.**

**Gating fact — unresolved, upstream:** do the productions carry positron transport (annihilation ≠
decay point)? The activity meta says "escaped positrons excluded" (⟹ transported, range IN); the LOR
schema calls z0 the "annihilation **(decay)** point" (⟹ decay=annihilation, range OUT). Cannot be
closed downstream. Coarse downstream check: truth per-isotope erfc edge widths give O-15 **9.19 mm**
≈ C-11 **9.42 mm** (O-15 not broader) → suggested the range term is small, but the erfc-over-window
is a poor proxy for the per-LOR variance (that was the wrong quantity — see below).

**Decisive test — per-isotope σ_R (`drivers/sigma_r_per_isotope.jl`):** thin the pooled master by
each isotope's posterior P(i | z0, t_decay) = r_i/Σr_j (the same per-event rates the washout weight
marginalises), reconstruct N realizations, take σ_R. **Natural statistics, NOT count-matched** (the
user's call): ¹⁵O enters with ~3× ¹¹C's counts, so counting alone predicts **σ_R(¹⁵O) ≈ 0.6·σ_R(¹¹C)**
— the count asymmetry is the built-in control. Reading:
- σ_R(¹⁵O) ≈ 0.6·σ_R(¹¹C) → range negligible;
- σ_R(¹⁵O) ≈ σ_R(¹¹C) (not 0.6×) → range is eating ¹⁵O's 3× count advantage — hypothesis supported;
- σ_R(¹⁵O) > σ_R(¹¹C) → range dominates.

Caveats: the posterior is a **soft** selector (isotope-blind data → leakage dilutes → any measured
difference is a lower bound); it only sees range that is actually in the annihilation points.

**Design / progress:**
- Scanner: **ring 1 m CsI** (both crystals have the same spatial resolution; BGO has more
  sensitivity + a bit more MS — my earlier "CsI better resolution" was wrong). CsI first because the
  effects show up more clearly there.
**DONE — ring 1 m CsI, N=100, 1 Gy** (O-15/C-11 only; N-13/O-14 fits blow up — posterior-selected
edges outside the fit window / too few counts). Figure
`…/comparison/figures/sigma_r_per_isotope.png` (`tools/plot_sigma_r_per_isotope.py`); data
`…/washout/sigma_r_per_isotope[_washed]_t{tstart}.toml`; combined = the ring-CsI washout nominal/washed.

- **Nominal σ_R vs t_start (t=0/120/180/300):** O-15 0.134/0.159/0.205/0.253 (rises as it drains —
  loses 5.5× its counts, 4.13M→749k), C-11 0.209/0.259/0.247/0.259 (flat after t=120 — barely loses
  counts), combined 0.105/0.119/0.158/0.209. O-15 < C-11 at every t_start (its counts).
- **Per-count precision k = σ_R·√N:** O-15 272/230/250/219, C-11 253/292/268/259 — no stable range
  penalty: the ~7% O-15 excess at t=0 (272 vs 253) reverses by t=120/180. Both ~260 ± ~20 (N=100 noise).
- **Washed, t=120/180:** inflation O-15 1.52/1.61, C-11 1.65/1.65, combined 1.57/1.54. Survival ~0.38
  → counting inflation ≈ 1/√0.38 ≈ 1.64; both isotopes sit there, **no isotope selectivity, no growing
  protection**; the ring's combined washed rises with the nominal (no "flat-washed" protection).

**Conclusion (Thread B): the positron-range hypothesis is NOT supported.** O-15 is never significantly
worse per count than C-11 (a ~7% t=0 hint that reverses), and washout is a plain ~1/√survival penalty
for both with no selectivity. So the R35/35 "protection" (Thread A) is **not** an isotope/range effect
— it is not reproduced on the ring, and is most likely R35/35 N=50 noise (R35/50 CsI N=200 would
confirm). The method and this negative result also motivate the per-species selector now in the
upstream washout spec (`latex/washout_brain.tex`, Sec. 7).

## State / provenance

- **Committed:** R35/35 geometry + R35/50 CsI 1 Gy; niter diagnostic + `--seed`; per-isotope driver +
  investigation doc; per-t_start filename convention; the per-isotope tomls (t=0/120/180/300 nominal,
  washed 120/180); and the upstream washout spec `latex/washout_brain.tex`.
- Drivers take `--tstart --dose --realizations [--seed] [--washed] [--isotopes]`; per-arm configs
  `config/run_parameters_{bgo,csi,r35_bgo,r35_csi,r35_35_bgo,r35_35_csi,chs_bgo,chs_csi}.toml`
  (`cp` to activate). Config stays on the arm under study (ring CsI for the per-isotope work) — the
  earlier BGO-restore habit was dropped.

## Next steps (optional / open)

1. R35/50 CsI N=200 — confirm whether the R35/35 protection is real/universal or N=50 noise (the one
   loose end left by Thread B's negative result).
2. Regenerate `sigma_r_configs.png` with the N=200 R35/35 point; fold the four-geometry comparison
   into the note.
5. If a real range effect emerges: one clean upstream question on whether the sim annihilates at the
   decay point or after positron transport.
