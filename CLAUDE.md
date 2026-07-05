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
  the frozen knobs, the validation ladder, the deferred register, and the build order. Start here.
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
- **Say "state plainly"** where you might reach for "honest / honestly" — describe the thing directly.

## Status

Follow `dev/PLAN.md` → "Build order". Done and next:

- **Step 1 — package scaffold: DONE.** `Project.toml` (package `CryspBrainSim`), `src/CryspBrainSim.jl`
  (`using RecoCryspTools`), `test/runtests.jl`. Deps resolved incl. Metal (GPU functional);
  `Pkg.test()` green. RecoCrysp is wired via local `[sources]` paths to the sibling checkout
  (`../RecoCrysp` @ `a6d900e`); pin that SHA into cached artifacts for provenance.
- **Step 2 — NEXT: port the endpoint estimator.** Write `src/endpoint.jl` + `src/profile.jl` as the
  Julia port of `py/depth_profile.py` (windowed, Poisson-weighted erfc fit) + `sigma_R` from
  `py/range_endpoint.py`, with the synthetic self-tests on shared arrays and the scipy-vs-LsqFit
  covariance check (ladder rungs 3–4). Independent of the products data.
- **Then** steps 3–8 in `dev/PLAN.md`.

Open threads: the `truth/` bundle request to PTCryspMC (`dev/upstream_request_truth_bundle.md`) is
pending — until it lands, `characterize.jl` reads from `ptcrysp-scenarios/`. Confirm the thinning
anchor `p = dose/top_dose` against the recipe's `dose_to_counts`.

First data target: `PtCryspProds/uniform_headep_sobp_1e8/crysp_ring_1m/bgo/fast_1Gy/lors_shard000.h5`
(shard 0, BGO).
