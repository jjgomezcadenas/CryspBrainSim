#!/usr/bin/env python3
"""
plot_one_shard.py — the two rung-5 figures from a drivers/one_shard.jl run:

  profile_shardNNN.png   reconstructed depth profile (disc ROI) with the
                         fitted erfc, the fixed window, and the truth
                         references (activity overlay, R50 markers, dose-R80)
  r50_iter_shardNNN.png  R50 vs MLEM iteration (the semi-convergence plateau)

Run:  python3 tools/plot_one_shard.py [shard_index]
Reads out/one_shard/recon_shardNNN.npz + results_shardNNN.toml and the
scenario truth/; writes into out/one_shard/figures/.
"""
import os
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.special import erfc  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCEN = os.path.join(os.path.dirname(REPO), "PtCryspProds", "uniform_headep_sobp_1e8")

BLUE, AQUA, RED = "#2a78d6", "#1baf7a", "#e34948"
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
    shard = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    tag = f"shard{shard:03d}"
    out = os.path.join(REPO, "out", "one_shard")
    d = np.load(os.path.join(out, f"recon_{tag}.npz"))
    with open(os.path.join(out, f"results_{tag}.toml"), "rb") as f:
        res = tomllib.load(f)
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)

    win = res["window_mm"]
    fit = res["final"]
    ref = res["reference"]

    # --- profile + fit + references
    z, prof = d["z"], d["profile"]
    fig, ax = plt.subplots(figsize=(9.5, 5), facecolor=SURFACE)
    style(ax)
    ax.axvspan(*win, color=GRIDC, alpha=0.5, lw=0)
    ax.plot(z, prof, drawstyle="steps-mid", color=BLUE, lw=1.6,
            label=f"reconstructed profile (ROI {res['roi_mm']:g} mm)")

    ta = np.genfromtxt(os.path.join(SCEN, "truth", "activity_profile_fast.csv"),
                       delimiter=",", names=True)
    sel = (ta["z_mm"] >= z[0]) & (ta["z_mm"] <= z[-1])
    scale = np.percentile(prof, 95) / ta["total"][sel].max()
    ax.plot(ta["z_mm"][sel], ta["total"][sel] * scale, ls="--", color=AQUA, lw=1.8,
            label="truth activity (scaled)")

    zf = np.linspace(win[0] - 10, win[1] + 8, 600)
    base_amp = prof[(z >= win[0]) & (z <= win[1])]
    # redraw the fitted erfc from (z0, w) with plateau/base from the window
    k = max(2, base_amp.size // 5)
    b0, p0 = np.median(base_amp[-k:]), np.median(base_amp[:k])
    ax.plot(zf, b0 + (p0 - b0) * 0.5 * erfc((zf - fit["r50_fit_mm"]) /
            (np.sqrt(2) * fit["w_mm"])), color=INK, lw=1.4, ls=":",
            label=f"erfc fit  R50 = {fit['r50_fit_mm']:.2f} ± {fit['z0_err_mm']:.2f} mm")

    for x, c, lab in ((fit["r50_fit_mm"], BLUE, None),
                      (ref["activity_R50_fit_mm"], AQUA, None),
                      (ref["dose_R80_mm"], RED, f"dose-R80 {ref['dose_R80_mm']:.2f}")):
        ax.axvline(x, color=c, ls=":", lw=1.2, label=lab)
    ax.set_xlim(-80, 20)
    ax.set_title(f"{tag}: reconstructed depth profile — R50 fit "
                 f"{fit['r50_fit_mm']:.2f} (truth fit {ref['activity_R50_fit_mm']:.2f}) / "
                 f"crossing {fit['r50_crossing_mm']:.2f} (truth {ref['activity_R50_crossing_mm']:.2f}) mm",
                 color=INK, fontsize=10, loc="left")
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("activity in ROI / slice", color=INK)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="lower left")
    fig.tight_layout()
    p1 = os.path.join(figdir, f"profile_{tag}.png")
    fig.savefig(p1, dpi=160, facecolor=SURFACE)
    print(f"wrote {p1}")

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
