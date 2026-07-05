#!/usr/bin/env python3
"""
origin_profile.py — the truth-origin depth histogram of a shard's TRUE
coincidences (validation ladder rung 4's quick-look): histogram the
annihilation z0 of the detected trues, overlay the scenario's true
activity(z) (the detector-independent profile the detected subset is drawn
from), and print the interpolated distal R50 of both. The Julia rung-4 fit
runs on the same quantity through fit_endpoint; this is its picture.

Run:  python3 tools/origin_profile.py --shard <.../lors_shard000.h5>
Writes out/origin_profile/figures/<scenario>_<crystal>_shardNNN_origin.png
(override with --out). --roi-mm R restricts to sqrt(x0²+y0²) ≤ R.
"""
import argparse
import csv
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BLUE, AQUA = "#2a78d6", "#1baf7a"  # detected-origin histogram, true activity
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
SURFACE = "#fcfcfb"


def dec(v):
    return v.decode() if isinstance(v, bytes) else v


def distal_crossing(z, y, level, reference=None):
    thr = level * (reference if reference is not None else y.max())
    for i in range(len(y) - 2, -1, -1):
        if y[i] >= thr > y[i + 1]:
            return z[i] + (thr - y[i]) * (z[i + 1] - z[i]) / (y[i + 1] - y[i])
    return float("nan")


def windowed_r50(z, y, window):
    """Half-height crossing read INSIDE the fixed distal window, against the
    window's own plateau and tail medians. The detected-origin histogram is
    tilted by the attenuation gradient along the head, so a global-max
    crossing misreads its edge; the local reading matches what the erfc fit
    (src/endpoint.jl) measures."""
    lo, hi = window
    sel = (z >= lo) & (z <= hi)
    zf, yf = z[sel], y[sel]
    if zf.size < 4:
        return float("nan")
    k = max(2, yf.size // 5)
    plateau = float(np.median(yf[:k]))
    base = float(np.median(yf[-k:]))
    return distal_crossing(zf, yf, 1.0, reference=base + 0.5 * (plateau - base))


def scenario_dir_of(shard_path):
    """<scenario>/ is three levels above a homogeneous config leaf; find it
    by walking up to the directory that holds truth/."""
    d = os.path.dirname(os.path.abspath(shard_path))
    for _ in range(5):
        if os.path.isdir(os.path.join(d, "truth")):
            return d
        d = os.path.dirname(d)
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shard", required=True)
    p.add_argument("--budget", default="fast")
    p.add_argument("--roi-mm", type=float, default=None,
                   help="transverse ROI radius on the origin (default: none)")
    p.add_argument("--bin-mm", type=float, default=1.0)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    with h5py.File(args.shard, "r") as f:
        a = {k: dec(v) for k, v in f.attrs.items()}
        xs = float(a["xyz_scale_mm"])
        truth = f["truth"][:]
        x0 = f["x0_mm"][:].astype(np.float64) * xs
        y0 = f["y0_mm"][:].astype(np.float64) * xs
        z0 = f["z0_mm"][:].astype(np.float64) * xs

    keep = truth == 0
    if args.roi_mm is not None:
        keep &= np.hypot(x0, y0) <= args.roi_mm
    z0 = z0[keep]
    scenario, crystal = a.get("scenario", "?"), a.get("crystal", "?")
    shard_ix = int(a.get("realization", -1))

    lo = np.floor(z0.min())
    edges = np.arange(lo, np.ceil(z0.max()) + args.bin_mm, args.bin_mm)
    hist, _ = np.histogram(z0, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(9, 4.8), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.plot(centers, hist, drawstyle="steps-mid", color=BLUE, lw=1.6,
            label=f"detected trues origin ({z0.size:,})")

    # Overlay the detector-independent truth activity, scaled to the
    # histogram plateau, when the scenario's truth/ bundle is reachable. The
    # fixed distal window comes from the truth edge (same margins as
    # distal_window in src/profile.jl).
    r50_true = float("nan")
    window = None
    sdir = scenario_dir_of(args.shard)
    if sdir:
        path = os.path.join(sdir, "truth", f"activity_profile_{args.budget}.csv")
        with open(path) as f:
            rows = list(csv.reader(f))
        za = np.array([float(r[0]) for r in rows[1:]])
        ta = np.array([float(r[-1]) for r in rows[1:]])
        r50_true = distal_crossing(za, ta, 0.5)
        window = (r50_true - 20.0, r50_true + 15.0)
        sel = (za >= centers[0]) & (za <= centers[-1])
        scale = hist.max() / ta.max()
        ax.plot(za[sel], ta[sel] * scale, color=AQUA, lw=1.8, ls="--",
                label="true activity (scaled)")

    # The detected-origin edge, read inside the fixed window (the attenuation
    # gradient tilts the histogram, so a global-max crossing misreads it).
    if window is None:
        window = (centers[0], centers[-1])
        ax.text(0.02, 0.02, "no truth/ bundle found — window = full range",
                transform=ax.transAxes, color=MUTED, fontsize=8)
    r50_det = windowed_r50(centers, hist.astype(float), window)
    ax.axvspan(*window, color=GRIDC, alpha=0.5, lw=0)

    for r, c in ((r50_det, BLUE), (r50_true, AQUA)):
        if np.isfinite(r):
            ax.axvline(r, color=c, ls=":", lw=1.2)
    ax.set_title(f"{scenario} / {crystal} — shard {shard_ix}: truth-origin depth "
                 f"(windowed R50: detected {r50_det:.2f} / true {r50_true:.2f} mm)",
                 color=INK, fontsize=10.5, loc="left")
    ax.set_xlabel("z0 [mm]", color=INK)
    ax.set_ylabel(f"trues / {args.bin_mm:g} mm", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)

    if not args.out:
        name = f"{scenario}_{crystal.lower()}_shard{shard_ix:03d}_origin.png"
        args.out = os.path.join(REPO, "out", "origin_profile", "figures", name)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=160, facecolor=SURFACE)
    print(f"wrote {args.out}")
    print(f"windowed R50 detected-origin = {r50_det:.3f} mm | "
          f"R50 true activity = {r50_true:.3f} mm")


if __name__ == "__main__":
    main()
