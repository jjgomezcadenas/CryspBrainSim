#!/usr/bin/env python3
"""
plot_washed_bgo_vs_csi.py — the crystal comparison for the note: washed σ_R vs
acquisition delay, BGO vs CsI, at each of the three scanner size-classes, from
drivers/sigma_r_v2.jl. Three panels (ring / 50 cm AFOV / 35 cm AFOV), each
overlaying the two crystals (points, ±1/√(2(N−1)), no connecting lines).

Size-class pairing is by AFOV (ring / 50 cm / 35 cm). BGO carries a cryostat, so
at every size-class it sits ~50 mm larger in radius than the CsI counterpart
(BGO ring r437 vs CsI r387; BGO r40 vs CsI R35 on the compact bores) — each is
the real scanner that crystal would build, so the crystal comparison is fair.

Reads each arm's washout_v2/sigma_r_washout_v2.toml and writes
washed_bgo_vs_csi_v2.png into <scenario>/<topology>/comparison/figures/. Run:
  python3 tools/plot_washed_bgo_vs_csi.py
"""
import math
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, scenario_out
from fit_activity_profile import BLUE, INK, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"

# size-class -> (title, (bgo_scanner, bgo_leaf), (csi_scanner, csi_leaf))
CLASSES = [
    ("ring 1 m", ("crysp_ring_1m_bgo_2x0", "bgo_195k_2X0"),
     ("crysp_ring_1m_csi_2x0", "csi_2X0")),
    ("50 cm AFOV", ("crysp_r40_50cm_bgo_2x0", "bgo_195k_2X0"),
     ("crysp_r35_50cm_csi_2x0", "csi_2X0")),
    ("35 cm AFOV", ("crysp_r40_35cm_bgo_2x0", "bgo_195k_2X0"),
     ("crysp_r35_35cm_csi_2x0", "csi_2X0")),
]


def washed(scanner, leaf):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, leaf),
                        "washout_v2", "sigma_r_washout_v2.toml")
    with open(path, "rb") as f:
        d = tomllib.load(f)
    pts = sorted(d["point"], key=lambda p: p["t_del_s"])
    return (np.array([p["t_del_s"] for p in pts]),
            np.array([p["washed_sigma_R_mm"] for p in pts]),
            1.0 / math.sqrt(2 * (d["realizations"] - 1)))


def main():
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.4), facecolor=SURFACE,
                             sharey=True)
    ymax = 0.0
    for ax, (title, bgo, csi) in zip(axes, CLASSES):
        style(ax)
        for (sc, leaf), col, mk, lab in ((bgo, INK, "o", "BGO 195 K"),
                                         (csi, RED, "s", "CsI")):
            t, s, band = washed(sc, leaf)
            dx = -4 if col is INK else 4
            ax.errorbar(t + dx, s, yerr=s * band, fmt=mk, ms=7, color=col,
                        capsize=3, label=lab)
            ymax = max(ymax, float(np.max(s + s * band)))
        ax.set_title(title, color=INK, fontsize=12, loc="left")
        ax.set_xlabel("$t_{\\rm del}$ [s]", color=INK, fontsize=12)
        ax.set_xticks([120, 180, 300])
        ax.set_xlim(90, 330)
    axes[0].set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    for ax in axes:
        ax.set_ylim(0.0, 1.08 * ymax)
    axes[0].legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    fig.suptitle("Washed $\\sigma_R$: BGO vs CsI at three scanner size-classes "
                 "(v2, 1 Gy)", color=INK, fontsize=12.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "washed_bgo_vs_csi_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
