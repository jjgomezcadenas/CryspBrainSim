#!/usr/bin/env python3
"""
plot_grogg_leaves_v2.py — the linear-endpoint (Grogg et al.) comparison across
the full acquisition-delay axis, on the two ambient/cold operating points that
carry all three leaves at their full 300 s windows: CsI(Tl) on R35/35 and
BGO 77 K on r40/35. Companion to plot_grogg_v2.py (which is the single-scan
mechanism figure); this one shows that the estimator ordering holds at every
delay and on both crystals.

Two panels (CsI(Tl), BGO 77 K), each plotting nominal σ_R vs t_del (log axis)
for the three estimators: the erfc R50 of this note, the raw linear intercept,
and the linear intercept with the published 7 mm smoothing. Points carry the
±1/√(2(N−1)) band. Reads each arm's washout_v2/sigma_r_grogg_v2.toml (three
points, one per leaf) from drivers/sigma_r_v2.jl.

Writes grogg_leaves_v2.png into <scenario>/<topology>/comparison/figures/. Run:
  python3 tools/plot_grogg_leaves_v2.py
"""
import math
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, scenario_out
from fit_activity_profile import AQUA, BLUE, INK, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"

# panel -> (title, scanner, crystal-label)
ARMS = [
    ("CsI(Tl), R35/35 (ambient)", "crysp_r35_35cm_csi_2x0", "csi_tl_2X0"),
    ("BGO 77 K, r40/35", "crysp_r40_35cm_bgo_2x0", "bgo_77k_2X0"),
]
# estimator -> (label, colour, marker, sigma-field)
ESTS = [
    ("$R_{50}$ (erfc)", BLUE, "o", "erfc_nominal_sigma_R_mm"),
    ("linear intercept, raw", RED, "^", "grogg_nominal_sigma_R_mm"),
    ("linear intercept, +7 mm", AQUA, "s", "grogg_sm7_nominal_sigma_R_mm"),
]


def load(scanner, crystal):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, crystal),
                        "washout_v2", "sigma_r_grogg_v2.toml")
    with open(path, "rb") as f:
        d = tomllib.load(f)
    pts = sorted(d["point"], key=lambda p: p["t_del_s"])
    band = 1.0 / math.sqrt(2 * (d["realizations"] - 1))
    t = np.array([p["t_del_s"] for p in pts])
    return t, pts, band


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), facecolor=SURFACE,
                             sharey=True)
    for ax, (title, scanner, crystal) in zip(axes, ARMS):
        style(ax)
        ax.set_yscale("log")
        t, pts, band = load(scanner, crystal)
        for lab, col, mk, field in ESTS:
            s = np.array([p[field] for p in pts])
            ax.errorbar(t, s, yerr=s * band, fmt=mk + "-", ms=7, lw=1.0,
                        color=col, capsize=3, label=lab)
        ax.set_title(title, color=INK, fontsize=12, loc="left")
        ax.set_xlabel("$t_{\\rm del}$ [s]", color=INK, fontsize=12)
        ax.set_xticks([120, 180, 300])
        ax.set_xlim(90, 330)
    axes[0].set_ylabel("nominal $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    axes[0].legend(frameon=False, fontsize=10.5, labelcolor=INK,
                   loc="center left")
    fig.suptitle("Range precision by estimator vs acquisition delay "
                 "(full 300 s windows, v2, 1 Gy)",
                 color=INK, fontsize=12.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))

    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "grogg_leaves_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
