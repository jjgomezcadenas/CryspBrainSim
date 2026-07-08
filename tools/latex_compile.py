#!/usr/bin/env python3
"""
latex_compile.py — compile a LaTeX note and clean up after it, leaving only
the source (.tex), the output (.pdf) and the figures (figs/).

Runs the standard bibliography sequence — pdflatex, bibtex, pdflatex twice —
so citations and cross-references resolve, then deletes the auxiliary files
LaTeX scatters (.aux, .log, .bbl, …). The .pdf is a build product and is not
tracked; regenerate it here.

Run:  python3 tools/latex_compile.py [name ...]
With no argument, compiles every .tex in latex/. A name may be a bare stem
(`cbs`), a file (`cbs.tex`) or a path.
"""
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LATEX = os.path.join(REPO, "latex")

AUX_EXT = (".aux", ".log", ".out", ".toc", ".lof", ".lot", ".bbl", ".blg",
           ".fls", ".fdb_latexmk", ".synctex.gz", ".nav", ".snm", ".vrb",
           ".idx", ".ilg", ".ind")


def run(cmd, cwd):
    """Run a compile step; return True on success. LaTeX/BibTeX return
    non-zero on warnings we tolerate, so failures are judged on the .pdf."""
    subprocess.run(cmd, cwd=cwd, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)


def compile_one(stem):
    tex = f"{stem}.tex"
    if not os.path.exists(os.path.join(LATEX, tex)):
        print(f"skip: {tex} not found in {LATEX}")
        return False
    run(["pdflatex", "-interaction=nonstopmode", tex], LATEX)
    run(["bibtex", stem], LATEX)
    run(["pdflatex", "-interaction=nonstopmode", tex], LATEX)
    run(["pdflatex", "-interaction=nonstopmode", tex], LATEX)
    for ext in AUX_EXT:
        f = os.path.join(LATEX, stem + ext)
        if os.path.exists(f):
            os.remove(f)
    pdf = os.path.join(LATEX, f"{stem}.pdf")
    ok = os.path.exists(pdf)
    print(f"{'wrote' if ok else 'FAILED'} {os.path.relpath(pdf, REPO)}")
    return ok


def main():
    args = sys.argv[1:]
    if not args:
        stems = sorted(f[:-4] for f in os.listdir(LATEX) if f.endswith(".tex"))
    else:
        stems = [os.path.splitext(os.path.basename(a))[0] for a in args]
    if not all(compile_one(s) for s in stems):
        sys.exit(1)


if __name__ == "__main__":
    main()
