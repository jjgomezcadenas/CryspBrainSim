# CLAUDE.md — CryspBrainSim

Orientation for any Claude Code session on this repo.

## Purpose

Take the list-mode LORs a PET scanner would record for a proton SOBP field, reconstruct the β⁺
activity, extract the distal range endpoint R, and report its statistical precision σ_R vs dose —
one curve per scanner geometry. This is the analysis end of the chain
`ptcrysp-scenarios → PTCryspMC.jl → PtCryspProds/ → here`; the reconstruction engine comes from
`RecoCryspTools` as a pinned dependency.

## Read first

- **`dev/PLAN.md`** — the build plan: structure, dependencies, the consume-vs-write code inventory,
  the frozen run parameters, the validation ladder, the deferred register, and the build order. Start here.
- **`dev/reference/`** — vendored snapshots of the upstream contracts (products tree, data-generation
  strategy, σ_R recipe, LOR schema, RecoCrysp usage), with provenance in `dev/reference/README.md`.

## Working agreement

- **Treat a question as a request for an answer.** When the user asks "why / what / where / can it",
  respond with the answer and stop. Take actions — edit, run, reconstruct, plot, commit — when the
  user gives an explicit instruction to do so ("do it", "run it", "build it", "commit"). When the
  intent reads as a question, answer it and name the next step in words.
- **Confirm before writing, building, committing, or pushing.** State what you are about to do and
  wait for the go-ahead.
- **Formulate affirmatively.** Say what to do, in as much detail as the reader needs to do it.
  Describe the action to take and the outcome to reach; let the obvious non-actions stay unsaid. This
  applies to plans, code comments, docstrings, and commit messages alike — `dev/PLAN.md` is written
  this way and stays that way.
- **Ask plainly, in prose.** State what you need to know directly, as a sentence.
- **Lead with measured numbers.** Give the figures first and the adjectives after; list the relevant
  options before any comparative claim.
- **State the path of any figure or artifact you produce** (e.g. `out/one_shard/figures/profile.png`)
  so it can be opened.
- **Write traceable scripts, never temporary ones.** Any script that produces an artifact — a
  figure, a reference file, a cache — lives in the repo (`tools/`, `test/`), so the artifact can be
  regenerated from a committed path. A scratchpad script that made a kept artifact is a lost
  provenance chain. When in doubt about where a script belongs, ask.
- **Always check before deleting; prefer rename over delete+add.** Never remove a file (or bundle a
  delete into a wider command) without first inspecting what is being removed and stating it. To
  rename or move a file, use `git mv` so the history reads as a rename, then apply new content on
  top — do not delete the old file and create a new one. Keep each delete as its own visible step.
- **Say "state plainly"** where you might reach for "honest / honestly" — describe the thing directly.

## Status

Follow `dev/PLAN.md` → "Build order". The build is complete and the physics studies are running;
this Status is the current-state orientation. Infrastructure and results below, then the next
problem (**biological washout**).

### Infrastructure (build steps 1–6 + output layout): DONE

The Julia package (`CryspBrainSim`, `using RecoCryspTools`, Metal GPU, tests green) provides the
whole chain: endpoint estimator (`profile.jl`, `endpoint.jl` — erfc fit, covariance = scipy
`absolute_sigma`), products navigation + provenance (`products.jl`, `shard_stats.jl`,
`characterize.jl`), attenuation + sensitivity (`mumap.jl`, `sensitivity.jl` — `sensitivity_base`
over ContinuousPET draws, NPZ+TOML cache with RecoCrysp SHA), Bernoulli thinning (`thinning.jl`),
the shared reconstruction context (`reconstruct.jl` — `load_run_context`, `write_descriptors`), the
`[configuration]`-driven run parameters (`config.jl`), and the output-layout path helpers
(`output.jl` mirrored by `tools/crysp_paths.py` — nothing hard-codes `out/<...>`; sensitivity at
the ring tier, crystal label folds material+thickness in X0). Drivers: `one_shard.jl`,
`ten_shards_dose.jl`, `ten_shards_tstart.jl`, `sigma_r_at_dose.jl`, `sigma_r_sweep_dose.jl`. Python
tools: `fit_activity_profile.py` (the fit lab — `--model`, `--roi`, `--no-baseline`, `--no-pulls`),
`ten_shards.py` (`--dose-sweep`/`--fano`/`--t-start`), `scatter_profile.py`,
`recon_scatters.jl`, `plot_recon_projections.py`, `plot_tstart.py`, `collect_note_figures.sh`,
`latex_compile.py`. **Switching scanner arm = one edit of `config/run_parameters.toml`
`[configuration]` (scenario/topology/scanner/crystal/leaf) + rerun (~25 min/arm).**

<details><summary>Original step-by-step build history (superseded by the summary above)</summary>

- **Step 2 — endpoint estimator port: DONE.** `src/profile.jl` (`depth_profile`, `distal_window`) +
  `src/endpoint.jl` (`fit_endpoint`, `sigma_R`), ported from `py/`. Covariance matches scipy
  `absolute_sigma=True` via `PrecisionWeights` (unweighted path matches the MSE-scaled default).
  Cross-validated on shared synthetic arrays against the frozen py reference
  (`test/data/endpoint_reference.npz`, regenerated by `test/make_reference.py`): profiles bit-exact,
  fit parameters ≤ 2e-5 relative, tolerances documented in `test/runtests.jl`; 36 tests green.
  Cross-check figure: `tools/plot_endpoint_port.py` → `out/endpoint_port/figures/fit_crosscheck.png`.
- **Step 3 — products, QA, characterization: DONE.** `src/products.jl` (leaf/shard navigation,
  provenance verification, pooling, geometry/phantom/truth readers), `src/qa.jl` (`ShardQA`),
  `src/characterize.jl` (`TruthReference`, `distal_crossing`); tools
  `plot_shard.py`/`plot_truth.py`/`shard_summary.py`/`origin_profile.py`; 85 tests green on
  synthetic fixtures. Rungs 1–2 ran on the real tree: reference locked at **dose-R80 −5.58 mm,
  activity-R50 −16.45 mm, offset −10.87 mm** (`out/characterize/truth_reference.toml`); all ten BGO
  shards pass QA (acceptance 21.73–21.75%, trues 81.1%, 0 degenerate). Rung-4 quick-look: erfc fit
  of the detected-origin profile −14.76 ± 0.03 mm vs −14.32 mm on the truth profile — the 0.4 mm
  tilt is the attenuation gradient (detected subset ≠ true activity; the recon corrects it).
- **Step 4 — mumap + sensitivity: DONE.** `src/mumap.jl` (analytic `attenuation_ellipsoid` +
  voxel `build_mumap`/`attenuation_mumap` routes, `centered_grid`), `src/sensitivity.jl`
  (`sensitivity_base` chunked over `ContinuousPET` draws, NPZ cache + TOML provenance with the
  RecoCrysp SHA, `scaled_sensitivity` applied per realization); 101 tests green.
  **Scale assessed** (`tools/bench_sensitivity.jl` → `out/sensitivity_scope/bench.toml`, M-series,
  6 threads, Metal): sample 34 M LORs/s, chords 934 M/s, backprojection 242 M/s (Metal) / 61 M/s
  (CPU) — **n_sens 10⁹ costs ~37 s**, 5×10⁸ ~20 s; memory is a non-issue (chunk 0.5 GiB, pooled
  174 M-event master ≈ 4.6 GiB for recon arrays, 48 GiB RAM). Two-seed MC mottle at the provisional
  64×64×96 @1.5 mm grid (`tools/make_sensitivity.jl --check`): 3.92% per image at 10⁸, 1.76% at
  5×10⁸, 1.24% at 10⁹ (exact 1/√n). **Run parameter set: n_sens = 10⁹** (PLAN.md table; base cached under
  `out/sensitivity/`), with the R50-vs-seed stage check at rung 5 confirming it at the frozen grid.
- **Step 5 — single-shard chain: DONE (rung 5; run parameters frozen).** `drivers/one_shard.jl` +
  `tools/plot_one_shard.py`. Shard 0's 14.14 M trues reconstructed on the corridor grid
  (64×64×96 @ 1.5 mm, origin (−47.25, −47.25, −119.25) — fixed from the first, clipping, centered
  grid; 0.06% of origins outside): **R50 fit −15.593 ± 0.147 mm / crossing −16.101 mm at the
  frozen 50 iterations** vs truth fit −14.32 / crossing −16.45 mm — both conventions sub-voxel,
  bracketing the truth (the residual is the disc-ROI vs full-plane profile convention, not
  attenuation). Plateau <0.01 mm per 10 iters past 40 (re-verified to 100 each run); stability
  spread 0.31 mm (ROI 12/15, window ±2 mm); 100 iters in 26 s on Metal → ~13 s per full-stat
  reconstruction at the frozen count. **Run parameters frozen into `config/run_parameters.toml`** (loaded by
  `load_run_parameters()`; drivers and tools consume it): grid, ROI 13 mm, window (−36.45, −1.45),
  niter 50, n_sens 10⁹ (corridor base: mottle 1.28%, 36.7 s, cached under `out/sensitivity/`).
  Results: `out/one_shard/results_shard000.toml`, figures under `out/one_shard/figures/`.
- **Step 6 — thinning + rung 6: DONE.** The σ_R drivers are `drivers/sigma_r_at_dose.jl`
  (σ_R at one dose: `--from-shards` reference or `--realizations N [--dose D]` thinned),
  `drivers/sigma_r_sweep_dose.jl` (σ_R across a dose grid), and shared `sigma_r_common.jl`;
  figures from `tools/plot_sigma_r.py` into `out/sigma_r/`. `src/thinning.jl` (`thin_mask`/`thin_lm`
  seeded Bernoulli over the pooled master, own seed namespace `THINNING_SEED_BASE = 1_000_000`;
  `dose_to_counts`, anchor confirmed — yield model collapses to `p = (dose/top_dose)/n_shards`
  because each shard is one top-dose acquisition).
  **σ_R(1 Gy) = 0.065 mm** from the 10 independent shards (`--from-shards`; fit; sem 0.021, ±24%),
  vs 0.239 mm crossing (3.7× noisier — the case for fitting). Mean R50 −15.553 mm, offset to
  dose-R80 −9.97 mm; each fit's own error 0.147 mm is 2.2× the true spread (MLEM correlations).
  **Gate PASSES**: thinned σ_R(1 Gy) = 0.078 mm (50 realizations) vs 0.065 mm — ratio 1.20, inside
  the ±51% 2σ band. Scatter effect (all events, uncorrected, shard 0): +65 μm fit / −6 μm crossing —
  correction deferred. 122 tests green.
- **Output layout: DONE (branch `output-layout`).** `out/` mirrors the products axes
  (`out/<scenario>/[truth,mumap]`, `.../<topology>/<ring>/sensitivity`,
  `.../<ring>/<crystal>/{shard_stats,one_shard,sigma_r,origin_profile}`, `out/validation/`) via path
  helpers in `src/output.jl` + `tools/crysp_paths.py` — nothing hard-codes `out/<...>`. Sensitivity
  is at the ring tier (geometric, shared across crystals); crystal folds material + thickness in X0
  (`bgo_3X0` = 3.7 cm BGO ≈ 3.3 X0). Renames: `qa`→`shard_stats` (`ShardStats`, `shard_stats`);
  `sensitivity_cache_name` drops the ring prefix (path carries it). Existing `out/` data relocated
  (not regenerated); all drivers + tools verified against the new tree; 138 tests green.

</details>

### Science results

- **Endpoint study, part (a): DONE (2026-07-08).** The study splits: (a) distal-edge estimation,
  (b) scanner comparison. **Settled protocol: whole-plane profiles (no ROI — the 13 mm disc clips
  the depth-widening halo, shifting R_p proximally 2.4–3.0 mm), erfc edge fit with FREE baseline
  (forcing b = 0 shifts R50 by 0.4–0.6 mm at our statistics; the fitted b is a small negative
  shape-slack, −1.5…−8% of amplitude, not a background), R50 (= fitted z0) as THE observable.**
  Its ~−11 mm offset to the dose R50 is a calibration constant fixed by the reference simulation;
  the measurement delivers variations against that anchor. Model cross-check (fit lab `--model
  erfc|sigmoid|both`, `--no-baseline`): the logistic sigmoid (Zapien-Campos Eq. 3 + baseline) gives
  identical R50 information — constant +0.10 mm offset, identical spreads, rung shifts equal to
  0.02 mm — and slightly better χ²; erfc stays primary, sigmoid is the built-in cross-check.
  R_p (tangent endpoint) demoted to a qualitative accuracy statement: it swings ~1.2 mm across
  rungs (the shard-0 "blur-stable" reading was a fluctuation; ten-shard mean shift truth→recon is
  +0.32 mm), 4× the shard spread of R50. R_x (1%) is a tail diagnostic only (model-dependent by
  ~7 mm). **Ten-shard ladder** (`tools/ten_shards.py` → `ten_shards/results.toml`, figures): Δ_R50
  = R50(act) − R50(dose), erfc: truth activity −10.743; origins (acceptance only) −10.988 ± 0.010
  [std 0.031]; recon(trues) −11.216 ± 0.018 [0.057]; recon(all events, uncorrected) −11.194 ± 0.022
  [0.071]. Budget: acceptance −0.25 mm, reconstruction −0.23 mm, scatters +0.02 mm (correction
  stays deferred). **Dose sweep** (`drivers/ten_shards_dose.jl` thins each shard, p = dose/1 Gy,
  seed idx = shard·1000 + 100·dose + seed; `tools/ten_shards.py --dose-sweep` →
  `ten_shards/dose_sweep.toml`): Δ_R50 dose-invariant (means −11.216/−11.207/−11.197/−11.248 at
  1.0/0.5/0.2/0.1 Gy, each within 1σ of the anchor) and σ_R follows 1/√dose exactly:
  0.057/0.072/0.117/0.177 mm (pred 0.081/0.127/0.180). **Test-dose statement: a single 0.1 Gy
  acquisition locates the distal edge to σ_R ≈ 0.18 mm with the calibration bias known to
  ≤ 0.04 mm** (trues, fast window, closed ring). Literature anchor: Zapien-Campos et al., Med Phys
  2025 (papers/, untracked) — logistic fit, PAR = R50, whole-plane, offset stable to 0.4–0.5 mm;
  their range-shift transfer test is their problem (IMPT spots), not ours (fixed field, dose axis).
  Upstream request written: `dev/upstream_request_lor_decay_time.md` (t_decay_s Float32 only;
  isotope id dropped as non-actionable).
- **Summary note: `latex/endpoint_precision.tex` (compiles clean).** Self-contained note on the
  estimator + precision: definitions (R50, R80, ΔR), setup, calibration budget (Table 1), scatter
  check, low-dose operation (1 Gy vs 0.1 Gy, exploratory-dose conclusion), outlook (delayed start,
  composite model, scanner comparison). Figures collected from `out/` into `latex/figs/` by
  `tools/collect_note_figures.sh` (fit_recon_activity, recon_projections, scatters_profile,
  ladder_delta_r50, dose_sweep_r50). New tools: `plot_recon_projections.py` (orthogonal projections
  + phantom outline + R50 plane), `--no-pulls` in the fit lab (publication figures; standing set on
  disk is pull-free, whole-plane, erfc).
- **Two-scanner campaign (part b, first round): DONE (2026-07-10).** Upstream published two
  ten-shard masters on the same 2X0 ring from identical annihilation sets (see
  `dev/reference/PRODUCTS.md`/`SCHEMA.md`, refreshed; production note: PTCryspMC.jl
  latex/scanner_prods.tex): **BGO 195 K** (15.3 M LORs/shard, 72.3/27.4/0.3% true/scatter/random,
  eres 15%, 413 keV, τ 5 ns) and **cryogenic CsI** (6.1 M, 86.6/13.3/0.1%, eres 6%, 472 keV,
  1.5 ns); both σ_xyz 1.486 mm (3.5 mm FWHM), with the requested `t_decay_s` column (delivered —
  the delayed-start study is unblocked). The old `crysp_ring_1m/bgo` master is **frozen** (methods
  reference only; the endpoint_precision.tex numbers refer to it). Plumbing: `[configuration]`
  block in `config/run_parameters.toml` names the active arm; Julia (`params.config`) and Python
  (`crysp_paths.active_config()`) read it — switching arms is one edit + the full chain
  (~25 min/arm). Labels from the arm name: `bgo_195k_2X0`, `csi_2X0`. Results (whole-plane erfc,
  free b, frozen chain): ladder Δ_R50 BGO −10.938 [0.023] / −11.165 [0.065] / −11.144 [0.113] and
  CsI −10.962 [0.039] / −11.232 [0.084] / −11.192 [0.114] (origins/recon/all-events, std in
  brackets); scatter shifts +0.021 / +0.040 mm (no correction needed for bias on either arm);
  dose sweeps flat with σ_R ∝ 1/√dose. **Headline: at the uncorrected all-events protocol the
  arms tie — k = σ_R(1 Gy) = 0.113 (BGO) vs 0.114 mm (CsI); BGO's 2.1× statistics cancels
  against its 2.1× scatter fraction. On trues they differ (0.065 vs 0.084 mm) — a scatter
  correction is the tiebreaker and the next technology question.** Also measured: paired shards
  buy nothing for R50 (detection sampling dominates; paired std = quadrature sum). Upstream nit:
  the BGO arm's scanner_geometry.json says material BGO_77K while shards say BGO_195K (stale key,
  harmless — constants come from shard attrs).
- **Delayed-start study: DONE (2026-07-10).** `shard_t_decay` (src/products.jl, asserts zero
  dropped LORs for alignment), `drivers/ten_shards_tstart.jl` (t_decay ≥ t_start cut, trues + all
  events, 20 recons/scanner/point), `tools/ten_shards.py --t-start`, `tools/plot_tstart.py`
  (cross-scanner figure → out/<scenario>/closed/comparison/figures/tstart_r50.png). Results
  (all events, Δ_R50 mean ± σ): kept 77/61/48% at t_start 60/120/180 s (identical both scanners);
  **σ_R stays ≈ 0.11 mm through 180 s on both scanners — at/below the counting prediction; the
  zero-delay precision survives the realistic in-room start.** Calibration walks −0.44/−0.99/−1.6 mm
  (¹⁵O drains from the mix; its lower threshold = deeper edge) → start-time sensitivity 8–10 μm/s,
  0.01 mm at 1 s timing. Edge-sharpening hypothesis refuted: w grows 10.8→11.4 mm (¹³N foot gains
  weight; ¹¹C positron-range gain invisible under the intrinsic width). In the note as §6 +
  Table 5 + Fig. 7. Scatter-correction stance settled: NOT needed at this level (calibration
  systematics dominate ≫ 0.1 mm); CsI's point is parity with less scatter sensitivity.
### The finding so far (endpoint_precision.tex — the two-scanner note, compiles clean)

Both representative scanners locate the distal edge with **σ_R ≈ 0.11 mm per 1 Gy run** at the
working protocol (all events, no corrections); precision scales as 1/√dose (0.16–0.23 mm at
0.1 Gy → the exploratory-dose case) and **survives a realistic in-room acquisition start** (still
~0.11 mm at t_start = 180 s). The activity–dose offset is a per-scanner calibration constant
measured to 0.01–0.04 mm/term; scatters bias it ≤ 0.04 mm (no correction needed at this level —
other calibration terms dominate). The note has 8 pages: setup + two-scanner table, estimator (R50
defined graphically), calibration budgets, the 1 Gy comparison, low-dose, §6 acquisition-start.

### Next problem: BIOLOGICAL WASHOUT (the current productions have NONE)

Every production to date models **physical decay only** — the β⁺ emitters sit where they were
created until they decay. In tissue they do not: dissolved gases and metabolites (¹⁵O, ¹¹C, ¹³N)
are cleared by perfusion/diffusion before decaying (literature: sum-of-exponentials fast/medium/
slow biological components per isotope/chemical form; folds with the physical half-life into an
effective one). Washout is thus the missing physics between our simulation and a real patient. Two
questions to answer, mirroring the endpoint study: **(a) does it degrade precision** — it removes
counts, worst at late times, compounding the delayed-start loss; **(b) does it bias the edge** —
only if clearance is depth-/tissue-dependent (in the uniform-brain phantom it is roughly uniform,
so first-order the edge should hold; the test is whether Δ_R50 moves under washout).

Open scoping question for the next session — **where washout enters**: (i) **upstream** — a new
production with washout kinetics applied (cleanest, but needs a PTCryspMC feature + run); or
(ii) **downstream reweighting** on the stored runs — an event decaying at `t_decay_s` must have
survived clearance from its creation to its decay, so it could be Bernoulli-thinned by a survival
probability. That needs the per-event **creation time and isotope** (the isotope id we explicitly
*dropped* from the decay-time request; creation time is not stored either), so downstream is not
currently possible without more columns. The **truth-level** washout study is computable now:
`truth/activity_profile_fast.csv` has per-isotope depth profiles (columns `O15,C11,N13,C10,O14,
total`), and **spatially-uniform washout is a per-isotope multiplicative reweighting** of those
columns (each isotope's P(z) scaled by its window-integrated survival fraction under its clearance
model), then refit for Δ_R50 — no new detected-level data needed. This is the natural first step;
it directly tests question (a) at truth level and, since the reweighting shifts the isotope mix,
whether the edge moves. Spatially-dependent clearance would need more than the depth columns.
Decide the axis, then whether an upstream request (or per-event isotope+creation-time columns) is
warranted.

### Pending (smaller, independent of washout)

- **latex/cbs.tex** (the living draft, separate from endpoint_precision.tex): revert eq:sigmaR from
  R_p back to z0/R50 (keep the R_p accuracy paragraph), cite Zapien-Campos, fold in the two-scanner
  numbers.
- **Composite-erfc edge model** (2–3 isotope components, offsets/widths frozen from per-isotope
  truth profiles, free amplitudes + global shift); adopt only if σ and rung stability improve
  incl. at 0.1 Gy. Directly relevant to washout (the per-isotope decomposition is the same object).

### Data on disk

Three ten-shard masters under `PtCryspProds/uniform_headep_sobp_1e8/`, all physical-decay-only:
the frozen reference `crysp_ring_1m/bgo/fast_1Gy/` (174.3 M LORs pooled; **do not** pool/compare
with the new arms) and the two 2X0 arms `crysp_ring_1m_bgo_2x0/bgo_195k/fast_1Gy/` (153 M) and
`crysp_ring_1m_csi_2x0/csi/fast_1Gy/` (61 M), plus the shared `truth/` bundle. Each new-arm shard
carries `t_decay_s` (absolute decay time; enables the pure-cut start-time study). Config
`[configuration]` is parked on the BGO 195 K arm.
