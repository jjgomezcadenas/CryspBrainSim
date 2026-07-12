#!/usr/bin/env python3
"""
plot_sigma_r_configs.py — how the washed σ_R varies across scanner
configurations, at the two acquisition-start times.

For each config (Ring 1m, R35/50, CHS) and crystal arm (BGO, CsI) it reads the
washout thinned firm-up (drivers/washout_sigma_r.jl --thinned) and takes the
washed σ_R at t_start = 120 and 180 s, scaled to 1 Gy (× √dose; the highest
available dose is used per point). Two panels (BGO | CsI, shared y), the three
configs on x, one series per start time (colour = start time) with the
±1/√(2(N−1)) thinned band.

Writes out/<scenario>/<topology>/comparison/figures/sigma_r_configs.png.
Run:  python3 tools/plot_sigma_r_configs.py
"""
import glob
import math
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, scenario_out
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
TIMES = [120.0, 180.0]
TCOLOUR = {120.0: BLUE, 180.0: RED}

# config label -> (BGO scanner/leaf, CsI scanner/leaf), ordered large bore -> compact
CONFIGS = [
    ("Ring 1m", ("crysp_ring_1m_bgo_2x0", "bgo_195k_2X0"),
                ("crysp_ring_1m_csi_2x0", "csi_2X0")),
    ("R35/50",  ("crysp_r35_50cm_bgo_2x0", "bgo_195k_2X0"),
                ("crysp_r35_50cm_csi_2x0", "csi_2X0")),
    ("R35/35",  ("crysp_r35_35cm_bgo_2x0", "bgo_195k_2X0"),
                ("crysp_r35_35cm_csi_2x0", "csi_2X0")),
    ("CHS",     ("crysp_chs_bgo_2x0", "bgo_195k_2X0"),
                ("crysp_chs_csi_2x0", "csi_2X0")),
]
ARMS = [("BGO 195 K", 1), ("CsI", 2)]   # (label, tuple index into a CONFIGS row)


def points(scanner, leaf):
    """t_start -> (nominal_1Gy, washed_1Gy, N) at 120/180 s, highest dose kept."""
    wdir = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, leaf), "washout")
    by_t = {}
    for f in glob.glob(os.path.join(wdir, "*.toml")):
        d = tomllib.load(open(f, "rb"))
        dose, nreal = d.get("dose_Gy"), d.get("realizations")
        if dose is None or "point" not in d:
            continue
        s = math.sqrt(dose)
        for p in d["point"]:
            t = p["t_start_s"]
            if t not in TIMES:
                continue
            if t not in by_t or dose > by_t[t][3]:
                by_t[t] = (p["nominal_sigma_R_mm"] * s,
                           p["washed_sigma_R_mm"] * s, nreal, dose)
    return by_t


def band(sig, n):
    return sig / math.sqrt(2 * (n - 1))


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2), facecolor=SURFACE,
                             sharey=True)
    x = np.arange(len(CONFIGS))
    labels = [c[0] for c in CONFIGS]

    for ax, (arm_name, idx) in zip(axes, ARMS):
        style(ax)
        for k, t in enumerate(TIMES):
            col = TCOLOUR[t]
            wsh, wsh_e = [], []
            for row in CONFIGS:
                scanner, leaf = row[idx]
                bt = points(scanner, leaf)[t]
                wsh.append(bt[1]);  wsh_e.append(band(bt[1], bt[2]))
            dx = (k - 0.5) * 0.08
            ax.errorbar(x + dx, wsh, yerr=wsh_e, fmt="o", ms=7, color=col,
                        ls="none", capsize=4,
                        label=f"$t_{{start}}$={t:.0f} s")
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_xlim(-0.5, len(CONFIGS) - 0.5)
        ax.set_title(arm_name, color=INK, fontsize=12, loc="left")
        ax.set_xlabel("scanner configuration", color=INK, fontsize=12)

    axes[0].set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    axes[0].set_ylim(0.0, None)
    axes[1].legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    fig.suptitle("Washed endpoint precision $\\sigma_R$ vs configuration "
                 "(two start times)", color=INK, fontsize=12)
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "sigma_r_configs.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
