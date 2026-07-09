#!/usr/bin/env python3
"""
scatter_profile.py — what shape do the scatters put under the distal edge?
Profiles the scatters-only reconstruction (tools/recon_scatters.jl) on the
whole plane and quantifies it against the trues reconstruction of the same
shard:

Each MLEM image is normalized to its own event count, so the scatters-only
profile is rescaled by n_scatters/n_all to express it in the units of the
joint all-events image — the contribution the all-events fit actually sits
on (the rescaled distal level reproduces the flat tail visible in the
all-events figure). Reported against the all-events profile:

  - the pedestal: the scaled scatter level distal of the edge, as a
    fraction of the all-events plateau;
  - the slope through the fit window: a straight line fitted to the scaled
    scatter profile inside the frozen window. A flat pedestal is exactly
    absorbed by the erfc's free baseline b and cannot move R50; only the
    sloping part can. The printed equivalent-shift is the WORST-CASE bound
    (window-scale variation divided by the edge gradient, as if the fit
    absorbed none of it) — compare it to the measured trues → all-events
    shift to see how much the free-baseline fit actually absorbs.

Run:  python3 tools/scatter_profile.py [--shard N] [--show]
Reads one_shard/recon_shardNNN_scatters.npz, recon_shardNNN_all_uncorr.npz
and results_shardNNN_all_uncorr.toml (for n_all); writes
one_shard/fits/scatters_shardNNN.toml and fits/figures/scatters_activity.png.
"""
import argparse
import os
import tomllib

import matplotlib
import numpy as np

from crysp_paths import REPO
from fit_activity_profile import (
    CFG, GRIDC, INK, RED, SURFACE, profile_from_image, style, toml_dump)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--show", action="store_true")
    args = p.parse_args()
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        params = tomllib.load(f)
    window = (params["window"]["z_lo_mm"], params["window"]["z_hi_mm"])
    centre = params["roi"]["centre_mm"]

    out = os.path.join(CFG, "one_shard")
    tag = f"shard{args.shard:03d}"
    ds = np.load(os.path.join(out, f"recon_{tag}_scatters.npz"))
    da = np.load(os.path.join(out, f"recon_{tag}_all_uncorr.npz"))
    with open(os.path.join(out, f"results_{tag}_all_uncorr.toml"), "rb") as f:
        n_all = tomllib.load(f)["n_events"]
    n_sc = int(ds["n_events"])
    z, ps = profile_from_image(ds["image"], params["grid"], centre, None)
    _, pa = profile_from_image(da["image"], params["grid"], centre, None)
    ps = ps * (n_sc / n_all)     # joint-image units (see docstring)

    # The all-events plateau (proximal fifth of the window) and edge
    # gradient: max |dP/dz| inside the window.
    inw = (z >= window[0]) & (z <= window[1])
    k = max(2, inw.sum() // 5)
    plateau = float(np.median(pa[inw][:k]))
    grad = float(np.max(np.abs(np.gradient(pa[inw], z[inw]))))

    # The scaled scatter profile inside the window: level + slope.
    (slope, level), cov = np.polyfit(z[inw] - z[inw].mean(), ps[inw], 1,
                                     cov=True)
    slope_err = float(np.sqrt(cov[0, 0]))
    # Distal pedestal: scaled scatter level beyond the window.
    distal = z > window[1] + 2
    pedestal = float(np.median(ps[distal]))

    # Worst-case bound: the window-scale variation of the scatter background
    # divided by the edge gradient, as if the fit absorbed none of it.
    span = window[1] - window[0]
    shift = float(slope) * span / grad

    res = {"shard": args.shard,
           "window_lo_mm": window[0], "window_hi_mm": window[1],
           "scatter": {"n_events": n_sc, "n_all": int(n_all),
                       "scale": n_sc / n_all,
                       "window_level": float(level),
                       "window_slope_per_mm": float(slope),
                       "window_slope_err": slope_err,
                       "distal_pedestal": pedestal},
           "all_events": {"plateau": plateau, "edge_gradient_per_mm": grad},
           "ratios": {"pedestal_over_plateau": pedestal / plateau,
                      "window_level_over_plateau": float(level) / plateau,
                      "worst_case_R50_shift_mm": shift}}
    print(f"{tag}: {n_sc} scatters of {n_all} events "
          f"(scale {n_sc / n_all:.3f})")
    print(f"scaled scatter pedestal {pedestal:.0f} = "
          f"{100 * pedestal / plateau:.1f}% of all-events plateau {plateau:.0f}")
    print(f"window: scatter level {level:.0f} "
          f"({100 * level / plateau:.1f}% of plateau), slope "
          f"{slope:+.2f} ± {slope_err:.2f} per mm over {span:.0f} mm")
    # The measured trues → all-events shift, from this arm's own ladder.
    ladder = os.path.join(CFG, "ten_shards", "results.toml")
    measured = ""
    if os.path.exists(ladder):
        with open(ladder, "rb") as f:
            t = tomllib.load(f)
        d = (t["all_ev"]["erfc"]["delta_R50_mean_mm"]
             - t["recon"]["erfc"]["delta_R50_mean_mm"])
        res["ratios"]["measured_shift_mm"] = d
        measured = f" (measured trues→all-events: {d:+.3f} mm)"
    print(f"all-events edge gradient {grad:.0f} per mm → worst-case R50 "
          f"shift bound ≈ {shift:+.3f} mm{measured}")

    fits = os.path.join(out, "fits")
    figdir = os.path.join(fits, "figures")
    os.makedirs(figdir, exist_ok=True)
    toml_dump(os.path.join(fits, f"scatters_{tag}.toml"), res)
    print(f"wrote {os.path.join(fits, f'scatters_{tag}.toml')}")

    fig, ax = plt.subplots(figsize=(9.5, 4.8), facecolor=SURFACE)
    style(ax)
    ax.axvspan(*window, color=GRIDC, alpha=0.5, lw=0)
    ax.errorbar(z, ps, yerr=np.sqrt(np.clip(ps, 1.0, None)), fmt="o", ms=3,
                color=INK, mec="none", elinewidth=0.8, capsize=0,
                label="scatters-only reconstruction")
    zf = z[inw]
    ax.plot(zf, level + slope * (zf - zf.mean()), color=RED, lw=1.0,
            label=f"window line: slope {slope:+.1f}/mm "
                  f"(≈ {shift:+.3f} mm on R50)")
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("P(z)", color=INK)
    ax.set_title(f"{tag}: whole-plane profile of the scatters-only "
                 f"reconstruction", color=INK, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="best")
    fig.tight_layout()
    fpath = os.path.join(figdir, "scatters_activity.png")
    fig.savefig(fpath, dpi=160, facecolor=SURFACE)
    print(f"wrote {fpath}")
    if args.show:
        plt.show(block=False)
        plt.pause(0.15)
        input("[return to continue] ")
    plt.close(fig)


if __name__ == "__main__":
    main()
