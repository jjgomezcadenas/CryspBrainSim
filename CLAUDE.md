# CLAUDE.md — general working rules

These are the general working rules for any session in this repo. Everything specific to
this repository — its purpose, what to read first, and the current project state — is in
the imported file below; read it before acting.

@ClaudeCryspBrainSim.md

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
- **Ask plainly, in prose. NEVER show me a menu.** State what you need to know directly, as a
  sentence — never as a multiple-choice picker or options menu. If there are choices, name them in a
  sentence and let me answer in my own words.
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
- **No LaTeX in the terminal; write readable formulas.** When talking in the terminal, never use
  LaTeX markup (no `\`, `$`, `_{}`, `^{}`, `\oplus`, `\sqrt`, etc.). Write formulas in plain readable
  form: `sigma_R`, `R80 - R50`, `1/sqrt(dose)`, `g(15O)/g(11C)`, `e^(-lambda t)`. LaTeX belongs only
  inside `.tex` files.
