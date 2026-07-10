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
`washout_sigma_r.jl` (isotope-washout σ_R: per-event thinning by the washout survival w(z₀,t_decay);
`--thinned` pools for N-realization σ_R at `--dose`/`--tstart`, `--washed-only` and the σ_R-blowup +
`min_washed_events` guards for the low-count corner).

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
