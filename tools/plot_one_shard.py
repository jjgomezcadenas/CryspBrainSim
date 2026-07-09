#!/usr/bin/env python3
"""
plot_one_shard.py — the reconstruction diagnostic from a
drivers/one_shard.jl run:

  r50_iter_shardNNN.png  R50 vs MLEM iteration (the semi-convergence plateau
                         that justifies the frozen iteration count)

The depth profile and its endpoint fit are the fit lab's job
(tools/fit_activity_profile.py, whole-plane erfc) — its figures land in the
same fits/figures/ directory, so every figure for one shard lives in one
place.

Run:  python3 tools/plot_one_shard.py [shard_index] [--all-uncorr]
Reads one_shard/recon_<tag>.npz (tag = shardNNN[_all_uncorr]); writes into
one_shard/fits/figures/.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from crysp_paths import active_config  # noqa: E402

ACTIVE = active_config()

BLUE = "#2a78d6"
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
SURFACE = "#fcfcfb"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def main():
    pos = [a for a in sys.argv[1:] if not a.startswith("--")]
    shard = int(pos[0]) if pos else 0
    tag = f"shard{shard:03d}" + ("_all_uncorr" if "--all-uncorr" in sys.argv else "")
    out = os.path.join(ACTIVE.cfg_dir, "one_shard")
    d = np.load(os.path.join(out, f"recon_{tag}.npz"))
    figdir = os.path.join(out, "fits", "figures")
    os.makedirs(figdir, exist_ok=True)

    # --- R50 vs iteration
    fig, ax = plt.subplots(figsize=(7.5, 4.2), facecolor=SURFACE)
    style(ax)
    ax.errorbar(d["iters"], d["r50s"], yerr=d["z0_errs"], color=BLUE, lw=1.8,
                marker="o", ms=4, capsize=2)
    ax.set_title(f"{tag}: R50 vs MLEM iteration (semi-convergence plateau)",
                 color=INK, fontsize=10.5, loc="left")
    ax.set_xlabel("MLEM iteration", color=INK)
    ax.set_ylabel("fitted R50 [mm]", color=INK)
    fig.tight_layout()
    p2 = os.path.join(figdir, f"r50_iter_{tag}.png")
    fig.savefig(p2, dpi=160, facecolor=SURFACE)
    print(f"wrote {p2}")


if __name__ == "__main__":
    main()
