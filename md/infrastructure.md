# Infrastructure (build steps 1–6 + output layout): DONE

The Julia package (`CryspBrainSim`, `using RecoCryspTools`, Metal GPU, tests green) provides the
whole chain:

- **Endpoint estimator** — `src/profile.jl` (`depth_profile`, `distal_window`), `src/endpoint.jl`
  (`fit_endpoint`, `sigma_R` — erfc fit, covariance = scipy `absolute_sigma` via `PrecisionWeights`;
  cross-validated against the frozen py reference `test/data/endpoint_reference.npz`). The chain
  reads **whole-plane** profiles (the settled protocol; `[roi]` in `run_parameters.toml` carries
  a centre but no radius since 2026-07-10) and matches the python fit lab to 1e-6 mm on the same
  image.
- **Products navigation + provenance** — `src/products.jl` (leaf/shard navigation,
  `REQUIRED_SHARD_ATTRS` verification, pooling, geometry/phantom/truth readers, `shard_t_decay`),
  `src/shard_stats.jl` (`ShardStats`), `src/characterize.jl` (`TruthReference`, `distal_crossing`).
- **Attenuation + sensitivity** — `src/mumap.jl` (`attenuation_ellipsoid`, voxel routes,
  `centered_grid`), `src/sensitivity.jl` (`sensitivity_base` over ContinuousPET draws, NPZ+TOML
  cache stamped with the RecoCrysp SHA, `scaled_sensitivity`). n_sens = 10⁹ (mottle 1.28%, ~37 s).
- **Thinning** — `src/thinning.jl` (`thin_mask`/`thin_lm` seeded Bernoulli over the pooled master,
  own seed namespace `THINNING_SEED_BASE = 1_000_000`; `dose_to_counts`, anchor confirmed
  `p = (dose/top_dose)/n_shards`).
- **Reconstruction context** — `src/reconstruct.jl` (`load_run_context`, `reconstruct_endpoint`,
  `write_descriptors`, `lor_attenuation`).
- **Run parameters** — `src/config.jl` (`load_run_parameters`) reads `config/run_parameters.toml`,
  including the `[configuration]` block (scenario/topology/scanner/crystal/leaf) that names the
  active arm.
- **Output layout** — `src/output.jl` mirrored by `tools/crysp_paths.py`; nothing hard-codes
  `out/<...>`. Tree mirrors the products axes: `out/<scenario>/[truth,mumap]`,
  `.../<topology>/<ring>/sensitivity` (ring tier — geometric, shared across crystals),
  `.../<ring>/<crystal>/{shard_stats,one_shard,ten_shards,sigma_r,origin_profile}`,
  `.../<topology>/comparison/` (cross-scanner), `out/validation/`. Crystal label folds
  material + thickness in X0 from the arm name (`bgo_3X0`, `bgo_195k_2X0`, `csi_2X0`).

**Drivers** (`drivers/`): `one_shard.jl`, `ten_shards_dose.jl`, `ten_shards_tstart.jl`,
`sigma_r_at_dose.jl`, `sigma_r_sweep_dose.jl` (`--all-events` reconstructs the uncorrected
working-protocol selection → `sweep_all.toml`; the thin seeds pair with the trues run),
`sigma_r_v2.jl` (authoritative nominal, exact isotope-labelled washout and per-isotope σ_R), and
`sigma_r_niter_v2.jl` (MLEM-iteration stability using the same v2 washout model).

**Python tools** (`tools/`): `fit_activity_profile.py` (the fit lab — `--model erfc|sigmoid|both`,
`--roi`, `--no-baseline`, `--no-pulls`; whole-plane erfc + free baseline is the default),
`ten_shards.py` (`--dose-sweep` / `--fano` / `--t-start`), `scatter_profile.py`,
`recon_scatters.jl`, `plot_recon_projections.py`, `plot_tstart.py`, `plot_one_shard.py`,
`plot_shard.py`, `plot_truth.py`, `plot_sigma_r.py`, `plot_chs_sigma_r.py` (single-shard
geometries CHS + R35 vs ring σ_R — six-arm dose sweeps + the ring ten-shard 1 Gy anchors;
`--crystal bgo|csi` for one-crystal figures), `collect_note_figures.sh`, `latex_compile.py`.

**Switching scanner arm** = activate the ready-made per-arm config,
`cp config/run_parameters_{bgo,csi}.toml config/run_parameters.toml`, then rerun the chain
(~25 min/arm). The two variants differ only in `[configuration]`; the frozen blocks are kept
identical (verified by md5). Everything resolves the active arm from `run_parameters.toml`.

## Generation-2 (v2) additions — DONE (2026-07-13)

To consume the v2 products (`dev/reference/generation2_plan.md`: tumour-centred phantom,
irradiation-end clock, per-LOR isotope column, stamped `washout_g`, `del…` scenario leaves):

- **Loader** — `src/products.jl`: `shard_isotope(file)` and `shard_generation(attrs)` readers;
  `pool_shards` now refuses to mix generations (a v2/legacy pool errors). Legacy loading is
  unchanged (`generation` absent ⇒ `"legacy"`).
- **Frame** — `src/characterize.jl`: `characterize(…; z_offset_mm=0.0)` rigidly shifts the truth
  z-frame (dose-R80, activity-R50, the fit window) into the reconstructed image frame.
  `src/reconstruct.jl`: `load_run_context` reads `source_z_offset_mm` from the shard (v2 only) and
  passes it — so the fit window lands on the tumour-centred edge (world z ≈ +9.13 mm) instead of the
  legacy −16.45. **Must-fix**: without it every v2 fit targets an empty window.
- **Grid/sensitivity** — the v2 configs carry the recentred grid (`img_origin` z −93.666, window
  `−10.868/24.132`, edge 9.132 = legacy + `source_z_offset_mm`); sensitivities rebuilt per scanner
  at that origin (cache `…orgm47.25_m47.25_m93.67…`, ~35 s each; `tools/make_sensitivity.jl`).
- **Driver** — `drivers/sigma_r_v2.jl`: loops the three `del…` leaves, pools each once, and produces
  in one pass the washout σ_R (nominal vs the **exact per-species g_i keep** — `p_dose·g_i[isotope]`
  from the isotope column + stamped `washout_g`, cross-checked vs a recompute from the stamped Mizuno
  params) and the **pure per-isotope** σ_R (isotope-column selection). Flags:
  `--realizations --dose --isotopes --leaves`.
- **Tools** — `tools/plot_sigma_r_v2.py` (the combined two-panel figure: nominal/washed σ_R + pure
  per-isotope, points-only; title reads the active arm), `tools/plot_washed_v2_scanners.py [csi|bgo]`
  (three-bore washed σ_R per crystal), `tools/plot_washed_bgo_vs_csi.py` (BGO-vs-CsI washed at the
  three size-classes). Wired into `collect_note_figures.sh` (16 figures).
- **Configs** — the two flagship reference arms are `config/run_parameters_csi_v2.toml`
  (ring CsI) and `config/run_parameters_ring_bgo_v2.toml` (ring BGO); small bores
  `run_parameters_{r35_35_csi,r35_50_csi,r40_35_bgo,r40_50_bgo}_v2.toml`. The 8 legacy off-centre
  configs are kept as archival provenance of the §1–7 results (their `fast_1Gy` products are gone;
  not re-runnable).
