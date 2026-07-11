# CLAUDE.md — CryspBrainSim

Orientation for any Claude Code session on this repo.

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
  [`isotope-washout.md`](md/isotope-washout.md) (the **current problem**),
  [`pending.md`](md/pending.md) (smaller items). Keep these current; this file points to them.
- **`dev/PLAN.md`** — the build plan: structure, dependencies, the consume-vs-write code inventory,
  the frozen run parameters, the validation ladder, the deferred register, and the build order.
- **`dev/reference/`** — vendored snapshots of the upstream contracts (products tree, data-generation
  strategy, σ_R recipe, LOR schema, RecoCrysp usage), with provenance in `dev/reference/README.md`.

## Working agreement

- **Treat a question as a request for an answer.** When the user asks "why / what / where / can it",
  respond with the answer and stop. Take actions — edit, run, reconstruct, plot, commit — when the
  user gives an explicit instruction to do so ("do it", "run it", "build it", "commit"). When the
  intent reads as a question, answer it and name the next step in words.
- **YU = "your understanding".** When the user closes a message with "YU", state your understanding
  of what was said back — verified against the repo/data where possible — and stop there.
- **Confirm before writing, building, committing, or pushing.** State what you are about to do and
  wait for the go-ahead.
- **Check when in doubt.** If carrying out an instruction literally would move, delete, rename, or
  overwrite something that already exists — or the literal reading might not be the intended one —
  ask first. Noticing the ambiguity and proceeding anyway is the error; a one-line question is cheap.
  Don't manufacture trivial doubts where an obvious default exists.
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
- **Every plot comes from a committed tool — never one on the fly.** No figure in a note (or anywhere
  kept) may exist that a committed `tools/` script does not produce and reproduce from its inputs.
  NEVER generate a plot ad hoc and drop the image into `latex/` or a note — that orphans it from its
  provenance. NEVER make a plot on the fly at all: always write (or extend) a tool, and **always
  check with me before making the plot.** Wire the figure into `tools/collect_note_figures.sh` so a
  fresh clone regenerates it.
- **Always check before deleting; prefer rename over delete+add.** Never remove a file (or bundle a
  delete into a wider command) without first inspecting what is being removed and stating it. To
  rename or move a file, use `git mv` so the history reads as a rename, then apply new content on
  top — do not delete the old file and create a new one. Keep each delete as its own visible step.
- **Say "state plainly"** where you might reach for "honest / honestly" — describe the thing directly.
- **Output in plain text or bold only.** Use **bold** for the emphasis you would otherwise carry
  with colour; no coloured text.

## Status

The build is complete; the physics studies are done through the acquisition-start axis and the
isotope-washout loss study. Details live in `md/` (see "Read first") — this is the one-paragraph
state.

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

**Smaller pending:** cbs.tex fold-in, composite-erfc model — [`md/pending.md`](md/pending.md).
