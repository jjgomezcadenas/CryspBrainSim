#!/usr/bin/env python3
"""
plot_recon_projections.py — the three orthogonal projections of a
reconstructed image (sum along each axis), with the phantom ellipsoid's
outline and the fitted R50 plane overlaid: the "what the scanner
reconstructs" figure. The activity corridor left by the beam sits inside
the head outline and ends at the distal edge the endpoint fit measures.

Run:  python3 tools/plot_recon_projections.py [--shard N] [--all-uncorr]
Reads one_shard/recon_shardNNN[_all_uncorr].npz, the frozen grid, the
phantom region (products phantom/phantom_regions.csv) and the R50 of the
matching fit TOML (if present); writes fits/figures/recon_projections.png.
"""
import argparse
import csv
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

from crysp_paths import REPO
from fit_activity_profile import ACTIVE, CFG, INK, MUTED, RED, SURFACE

PRODUCTS = ACTIVE.scenario_dir


def phantom_region(region=0):
    with open(os.path.join(PRODUCTS, "phantom", "phantom_regions.csv")) as f:
        row = list(csv.DictReader(f))[region]
    return {"semi_axes": (float(row["a_mm"]), float(row["b_mm"]),
                          float(row["c_mm"])),
            "centre": (float(row["cx_mm"]), float(row["cy_mm"]),
                       float(row["cz_mm"]))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--all-uncorr", action="store_true")
    args = p.parse_args()

    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        grid = tomllib.load(f)["grid"]
    org, vs, n = grid["img_origin_mm"], grid["voxsize_mm"], grid["n"]
    lo = [org[i] - vs[i] / 2 for i in range(3)]
    hi = [org[i] + (n[i] - 0.5) * vs[i] for i in range(3)]

    tag = f"shard{args.shard:03d}" + ("_all_uncorr" if args.all_uncorr else "")
    out = os.path.join(CFG, "one_shard")
    img = np.load(os.path.join(out, f"recon_{tag}.npz"))["image"]

    # The fitted R50 of the matching whole-plane fit, when available.
    fit_toml = os.path.join(out, "fits", f"fit_{tag}.toml")
    r50 = None
    if os.path.exists(fit_toml):
        with open(fit_toml, "rb") as f:
            t = tomllib.load(f)
        act = ("recon_all_events_activity" if args.all_uncorr
               else "recon_activity")
        r50 = t.get(act, {}).get("erfc", {}).get("R50_mm")

    ph = phantom_region()
    ax_mm, cen = ph["semi_axes"], ph["centre"]

    # (title, summed axis, horizontal axis idx, vertical axis idx)
    panels = [("transaxial (x-y)", 2, 0, 1),
              ("coronal (x-z)", 1, 2, 0),
              ("sagittal (y-z)", 0, 2, 1)]
    widths = [hi[h] - lo[h] for _, _, h, _ in panels]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.4), facecolor=SURFACE,
                             width_ratios=widths)
    for ax, (title, s, h, v) in zip(axes, panels):
        # Summing over axis s leaves the remaining axes in index order
        # (rows, cols); imshow puts cols on the horizontal axis.
        proj = img.sum(axis=s)
        rem = [i for i in range(3) if i != s]
        if rem[1] != h:
            proj = proj.T
        ax.imshow(proj, origin="lower", cmap="Blues", aspect="equal",
                  extent=[lo[h], hi[h], lo[v], hi[v]])
        ax.add_patch(Ellipse((cen[h], cen[v]), 2 * ax_mm[h], 2 * ax_mm[v],
                             fill=False, ls="--", lw=1.2, ec=INK, alpha=0.8))
        if r50 is not None and h == 2:
            ax.axvline(r50, color=RED, ls=":", lw=1.2)
        ax.set_title(title, color=INK, fontsize=12, loc="left")
        ax.set_xlabel("xyz"[h] + " [mm]", color=INK, fontsize=12)
        ax.set_ylabel("xyz"[v] + " [mm]", color=INK, fontsize=12)
        ax.tick_params(colors=MUTED, labelsize=10)
        for sp in ax.spines.values():
            sp.set_color(MUTED)

    fig.suptitle(f"{ACTIVE.scanner} / {ACTIVE.label} — run "
                 f"{tag.replace('shard', '').replace('_all_uncorr', ' (all events)')}",
                 color=INK, fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    figdir = os.path.join(out, "fits", "figures")
    os.makedirs(figdir, exist_ok=True)
    name = "recon_projections" + ("_all_uncorr" if args.all_uncorr else "")
    path = os.path.join(figdir, f"{name}.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
