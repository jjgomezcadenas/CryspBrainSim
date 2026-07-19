#!/usr/bin/env python3
"""
plot_washed_shortscan_afov.py — washed range precision at the short in-room scan
([120,300] s window, 1 Gy) vs scanner AFOV, for the two crystals. One point per
scanner size (TBP / LAFOV / CAFOV, the paper's naming), CsI and BGO as two
series with the ±1/√(2(N−1)) band. Reads each arm's
washout_v2/sigma_r_washout_v2_t120_300.toml (field washed_sigma_R_mm).

Writes washed_shortscan_afov.png into <scenario>/<topology>/comparison/figures/.
Run:  python3 tools/plot_washed_shortscan_afov.py
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

# scanner size -> (x position, tick label); ordered largest -> smallest AFOV
SIZES = [
    (0, "TBP\n(100 cm)"),
    (1, "LAFOV\n(50 cm)"),
    (2, "CAFOV\n(35 cm)"),
]
# crystal -> (label, colour, marker, [(scanner, crystal-label) per size])
ARMS = [
    ("CsI", RED, "s", [
        ("crysp_ring_1m_csi_2x0", "csi_2X0"),
        ("crysp_r35_50cm_csi_2x0", "csi_2X0"),
        ("crysp_r35_35cm_csi_2x0", "csi_2X0"),
    ]),
    ("BGO", BLUE, "o", [
        ("crysp_ring_1m_bgo_2x0", "bgo_195k_2X0"),
        ("crysp_r40_50cm_bgo_2x0", "bgo_195k_2X0"),
        ("crysp_r40_35cm_bgo_2x0", "bgo_195k_2X0"),
    ]),
]


def washed_short(scanner, crystal):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, crystal),
                        "washout_v2", "sigma_r_washout_v2_t120_300.toml")
    with open(path, "rb") as f:
        d = tomllib.load(f)
    p = d["point"][0]
    band = 1.0 / math.sqrt(2 * (d["realizations"] - 1))
    return p["washed_sigma_R_mm"], band


def main():
    fig, ax = plt.subplots(figsize=(6.4, 4.6), facecolor=SURFACE)
    style(ax)
    xs = [x for x, _ in SIZES]
    for lab, col, mk, arms in ARMS:
        s, bands = [], []
        for scanner, crystal in arms:
            v, b = washed_short(scanner, crystal)
            s.append(v)
            bands.append(b)
        s = np.array(s)
        dx = -0.03 if col is RED else 0.03
        ax.errorbar(np.array(xs) + dx, s, yerr=s * np.array(bands), fmt=mk + "-",
                    ms=8, lw=1.2, color=col, capsize=3, label=lab)
    ax.set_xticks(xs)
    ax.set_xticklabels([lab for _, lab in SIZES], color=INK, fontsize=11)
    ax.set_xlim(-0.4, 2.4)
    ax.set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    ax.set_title("Short in-room scan ([120,300] s): $\\sigma_R$ vs scanner AFOV",
                 color=INK, fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=12, labelcolor=INK, loc="center left")
    fig.tight_layout()

    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "washed_shortscan_afov.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
