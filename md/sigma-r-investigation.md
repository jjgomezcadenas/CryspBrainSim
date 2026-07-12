# σ_R across geometries + the positron-range variance question (IN PROGRESS, 2026-07-12)

Live investigation, opened after the washout loss study. Two linked threads: (a) a
four-configuration σ_R comparison, and (b) whether **positron range** sets a per-isotope
σ_R floor. This file is the compaction anchor; numbers + tools below.

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
- **Phase 1 (DONE, ring 1 m CsI, t=0, nominal, N=100):** O-15 N=4.13M σ_R **0.134 ± 0.010**;
  C-11 N=1.47M σ_R **0.209 ± 0.015**. Count ratio 2.82× → counting predicts σ_R(O-15)=0.125;
  observed 0.134 (≈1σ, ~7% high). Observed σ_R(O-15)/σ_R(C-11)=0.64 vs counting-only 0.60. So
  **O-15 tracks its count advantage — the range penalty is small and not significant** (at most a
  ~7% ~1σ hint; soft-selector leakage means the true term could be a bit larger, still small). This
  matches the earlier truth-profile read (O-15 edge not broader than C-11). N-13/O-14 fits blew up
  (posterior-selected edges outside the fit window / too few counts) — only O-15 vs C-11 is clean.
  Result: `…/crysp_ring_1m_csi_2x0/…/washout/sigma_r_per_isotope.toml`.
- **Phase 2 (planned):** t=120/180, nominal AND washed — watch each isotope's σ_R and the mix as
  O-15 drains and washout bites.

## State / provenance

- **Committed:** `ed6a3f5` (R35/35 geometry + R35/50 CsI 1 Gy), `27dad1a` (niter diagnostic + a
  `--seed` flag on `washout_sigma_r.jl`).
- **Uncommitted:** N=200 R35/35 CsI toml (overwrote the N=50 on disk); `drivers/sigma_r_per_isotope.jl`;
  the 4-config figure not yet regenerated with the N=200 R35/35 point.
- Drivers take `--tstart --dose --realizations [--seed] [--washed] [--isotopes]`; per-arm configs
  `config/run_parameters_{bgo,csi,r35_bgo,r35_csi,r35_35_bgo,r35_35_csi,chs_bgo,chs_csi}.toml`
  (`cp` to activate). Config is parked on ring BGO after each run.

## Next steps

1. Read Phase 1 (per-isotope t=0) → does ¹⁵O reach the ~0.6× that counting predicts, or does range
   erode it?
2. Phase 2: t=120/180 ± washout, per isotope.
3. Optionally R35/50 CsI N=200 (confirm the protection's universality).
4. Regenerate `sigma_r_configs.png` with the N=200 R35/35 point; commit the per-isotope driver +
   results.
5. If a real range effect emerges: one clean upstream question on whether the sim annihilates at the
   decay point or after positron transport.
