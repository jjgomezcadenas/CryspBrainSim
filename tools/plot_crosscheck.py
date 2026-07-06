#!/usr/bin/env python3
"""
plot_crosscheck.py — the rung-6 figure from a drivers/shard_crosscheck.jl
run: per-shard R50 in both conventions with the ensemble mean and ±σ_R band.

Run:  python3 tools/plot_crosscheck.py
Reads out/shard_crosscheck/endpoints.npz + results.toml; writes
out/shard_crosscheck/figures/r50_by_shard.png.
"""
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BLUE, AQUA = "#2a78d6", "#1baf7a"
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


def panel(ax, shards, vals, errs, stats, color, title):
    mean, sig = stats["mean_mm"], stats["sigma_mm"]
    ax.axhspan(mean - sig, mean + sig, color=GRIDC, alpha=0.6, lw=0)
    ax.axhline(mean, color=MUTED, ls="--", lw=1)
    if errs is not None:
        ax.errorbar(shards, vals, yerr=errs, fmt="o", ms=5, color=color,
                    capsize=2, lw=1.4)
    else:
        ax.plot(shards, vals, "o", ms=5, color=color)
    ax.set_title(f"{title}: mean {mean:.3f} mm, σ_R {sig:.3f} mm",
                 color=INK, fontsize=10, loc="left")
    ax.set_xlabel("shard", color=INK)
    ax.set_xticks(shards)


def main():
    out = os.path.join(REPO, "out", "shard_crosscheck")
    d = np.load(os.path.join(out, "endpoints.npz"))
    with open(os.path.join(out, "results.toml"), "rb") as f:
        res = tomllib.load(f)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)
    for a in (a1, a2):
        style(a)
    panel(a1, d["shards"], d["r50_fit"], d["z0_err"], res["sigma_R_fit"],
          BLUE, "erfc fit (bars: per-fit z0_err)")
    a1.set_ylabel("R50 [mm]", color=INK)
    panel(a2, d["shards"], d["r50_cross"], None, res["sigma_R_crossing"],
          AQUA, "windowed crossing")

    fig.suptitle(f"Ten-shard cross-check at 1 Gy — bias-free σ_R "
                 f"(n = {res['n_shards']}, std rel. error ≈ 24%); grey band = ±σ_R",
                 color=INK, fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)
    path = os.path.join(figdir, "r50_by_shard.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
