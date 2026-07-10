#!/usr/bin/env python3
"""
plot_chs_sigma_r.py — the single-shard geometries (the compact head scanner
CHS, r 200 mm, and the R35 half-metre ring, r 350 mm / AFOV 512 mm) against
the 1 m ring: σ_R vs dose for the six arms (BGO 195 K and CsI on each
geometry), one panel per selection (trues-only and the all-events working
protocol).

Reads each arm's sigma_r/sweep.npz (trues) and sweep_all.npz (all events)
— drivers/sigma_r_sweep_dose.jl, whole-plane erfc with free baseline — and,
for the ring arms, the ten-shard measurements at 1 Gy (ten_shards/
dose_sweep.toml for trues, ten_shards/results.toml for all events) as the
independent anchor the 1/√dose extrapolation is checked against. The CHS
and R35 arms hold a single shard (realization 0), so their σ_R points come
from thinned realizations at ≤ 0.2 Gy and the 1 Gy value is the
extrapolation line alone.

Writes the cross-scanner figure into the topology-level comparison/ leaf
(chs_sigma_r.png; with --crystal bgo|csi, the one-crystal figure
chs_sigma_r_<crystal>.png).

Run:  python3 tools/plot_chs_sigma_r.py [--crystal bgo|csi|both]
"""
import argparse
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, crystal_label, scenario_out
from fit_activity_profile import BLUE, INK, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
# (label, scanner dir, crystal arm, wall mm, colour, marker, line style,
#  ten-shard 1 Gy anchor?): the six arms. Open squares/dashed = the 1 m ring
# (the ten-shard reference), filled circles/solid = CHS, open
# triangles/dotted = R35.
ARMS = [
    ("ring BGO 195 K", "crysp_ring_1m_bgo_2x0", "bgo_195k", 22.36, RED,
     "s", "--", True),
    ("CHS BGO 195 K", "crysp_chs_bgo_2x0", "bgo_195k", 22.36, RED,
     "o", "-", False),
    ("R35 BGO 195 K", "crysp_r35_50cm_bgo_2x0", "bgo_195k", 22.36, RED,
     "^", ":", False),
    ("ring CsI", "crysp_ring_1m_csi_2x0", "csi", 37.20, BLUE,
     "s", "--", True),
    ("CHS CsI", "crysp_chs_csi_2x0", "csi", 37.20, BLUE,
     "o", "-", False),
    ("R35 CsI", "crysp_r35_50cm_csi_2x0", "csi", 37.20, BLUE,
     "^", ":", False),
]
PANELS = [("trues-only", "sweep"), ("all events, uncorrected", "sweep_all")]


def arm_cfg(scanner, crystal, wall_mm):
    return config_out(SCENARIO, TOPOLOGY, scanner, crystal_label(crystal, wall_mm))


def sweep_points(cfg, stem):
    d = np.load(os.path.join(cfg, "sigma_r", stem + ".npz"))
    dose, sig, n = d["dose_Gy"], d["sigma_fit_mm"], d["n_ok"]
    return dose, sig, sig / np.sqrt(2 * (n - 1))


def ring_1gy(cfg, stem):
    """The ten-shard 1 Gy measurement (n = 10 independent shards)."""
    if stem == "sweep":
        with open(os.path.join(cfg, "ten_shards", "dose_sweep.toml"), "rb") as f:
            s = tomllib.load(f)["d1000mGy"]["erfc"]["sigma_R_mm"]
    else:
        with open(os.path.join(cfg, "ten_shards", "results.toml"), "rb") as f:
            s = tomllib.load(f)["all_ev"]["erfc"]["delta_R50_std_mm"]
    return s, s / np.sqrt(2 * (10 - 1))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crystal", choices=("bgo", "csi", "both"), default="both")
    args = p.parse_args()
    arms = [a for a in ARMS if args.crystal == "both"
            or a[2].startswith(args.crystal)]

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.6), facecolor=SURFACE,
                             sharex=True)
    dd = np.geomspace(0.04, 1.2, 100)
    for ax, (title, stem) in zip(axes, PANELS):
        style(ax)
        for label, scanner, crystal, wall, colour, mk, ls, tenshard in arms:
            cfg = arm_cfg(scanner, crystal, wall)
            dose, sig, err = sweep_points(cfg, stem)
            mfc = colour if mk == "o" else "none"
            ax.errorbar(dose, sig, yerr=err, fmt=mk, ms=6, color=colour,
                        mfc=mfc, mew=1.2, elinewidth=1.2, capsize=3,
                        label=label)
            # One thin 1/√dose line per arm (weighted anchor a = σ·√d).
            a = np.average(sig * np.sqrt(dose), weights=1 / err**2)
            ax.plot(dd, a / np.sqrt(dd), color=colour, lw=0.9, ls=ls,
                    alpha=0.6)
            if tenshard:
                s1, e1 = ring_1gy(cfg, stem)
                ax.errorbar([1.0], [s1], yerr=[e1], fmt="*", ms=11,
                            color=colour, mfc="none", mew=1.2, elinewidth=1.2,
                            capsize=3,
                            label=None if scanner != arms[0][1]
                            else "10-shard measurement (1 Gy)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0],
                      ["0.05", "0.1", "0.2", "0.5", "1"])
        ax.set_ylabel("$\\sigma_R$ [mm]", color=INK, fontsize=13)
        ax.set_title(title, color=INK, fontsize=11, loc="left")
        ax.legend(frameon=False, fontsize=10, labelcolor=INK, loc="lower left",
                  ncol=2)
    axes[1].set_xlabel("dose [Gy]", color=INK, fontsize=13)
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    tag = "" if args.crystal == "both" else f"_{args.crystal}"
    path = os.path.join(out, f"chs_sigma_r{tag}.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
