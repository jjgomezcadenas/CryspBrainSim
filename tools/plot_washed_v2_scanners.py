#!/usr/bin/env python3
"""
plot_washed_v2_scanners.py — the generation-2 geometry comparison for the note:
washed σ_R vs acquisition delay, for the three CsI bores (ring 1 m, R35/50,
R35/35), from drivers/sigma_r_v2.jl. One panel, washed σ_R only (points,
±1/√(2(N−1)), no connecting lines) — the intrinsic/nominal σ_R, the inflation,
and the per-isotope split live in the per-arm tables/figure, not here.

Reads each arm's washout_v2/sigma_r_washout_v2.toml and writes
washed_sigma_r_scanners_v2.png into <scenario>/<topology>/comparison/figures/.
Run:  python3 tools/plot_washed_v2_scanners.py
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
LEAF = "csi_2X0"
# (scanner dir, display label, colour, marker) — the three CsI bores, v2.
ARMS = [
    ("crysp_ring_1m_csi_2x0", "ring 1 m ($r$ 387 mm)", INK, "o"),
    ("crysp_r35_50cm_csi_2x0", "R35/50 ($r$ 350, AFOV 500)", BLUE, "s"),
    ("crysp_r35_35cm_csi_2x0", "R35/35 ($r$ 350, AFOV 350)", RED, "^"),
]


def washed(scanner):
    """(t_del array, washed σ_R array, band) for one arm."""
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, LEAF),
                        "washout_v2", "sigma_r_washout_v2.toml")
    with open(path, "rb") as f:
        d = tomllib.load(f)
    pts = sorted(d["point"], key=lambda p: p["t_del_s"])
    t = np.array([p["t_del_s"] for p in pts])
    s = np.array([p["washed_sigma_R_mm"] for p in pts])
    return t, s, 1.0 / math.sqrt(2 * (d["realizations"] - 1))


def main():
    fig, ax = plt.subplots(figsize=(7.4, 5.2), facecolor=SURFACE)
    style(ax)
    for scanner, label, col, mk in ARMS:
        t, s, band = washed(scanner)
        ax.errorbar(t + (0 if col is INK else (6 if col is BLUE else -6)), s,
                    yerr=s * band, fmt=mk, ms=8, color=col, capsize=3, label=label)
    ax.set_xlabel("acquisition delay $t_{\\rm del}$ [s]", color=INK, fontsize=13)
    ax.set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    ax.set_ylim(0.0, None)
    ax.set_xticks([120, 180, 300])
    ax.set_xlim(90, 330)
    ax.legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    ax.set_title("Washed $\\sigma_R$ vs acquisition delay — three CsI bores (v2, 1 Gy)",
                 color=INK, fontsize=12, loc="left")
    fig.tight_layout()

    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "washed_sigma_r_scanners_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
