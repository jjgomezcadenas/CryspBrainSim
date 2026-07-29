#!/usr/bin/env python3
"""
plot_statproc_delay_bgo.py — reference BGO-ring washed range precision sigma_R vs
acquisition start delay, for two scan durations, on the data-driven source
(nominal fitted cross sections; bounded erfc fit, finite-pool-corrected, N=100,
1 Gy). Two series:
  red  — full 300 s scan  ([t, t+300]): del120/del180/del300 leaves
  blue — short 120 s scan ([t, t+120]): the --tend sub-cuts _t120_240, _t180_300
at start delays t = 120, 180, 300 s (the 120 s scan has t = 120, 180 only).
Points with the +/-1/sqrt(2(N-1)) band stored in each combine.

Reads each case's combined/washed_N100.toml (fields corrected_sigma_R_mm,
relative_sigma_uncertainty) under crysp_ring_1m_bgo_2x0/bgo_195k_2X0/
statistical_procedure/.

Writes latex/figs/statproc_delay_bgo.png.
Run:  python3 tools/plot_statproc_delay_bgo.py
"""
import glob
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out
from fit_activity_profile import BLUE, INK, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8_dd", "closed"
SCANNER, CRYSTAL = "crysp_ring_1m_bgo_2x0", "bgo_195k_2X0"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# series -> (label, colour, marker, [(start, case-dir) ...])
SERIES = [
    ("300 s scan", RED, "o", [
        (120, "del120s_ac300s_1Gy_D1p0Gy"),
        (180, "del180s_ac300s_1Gy_D1p0Gy"),
        (300, "del300s_ac300s_1Gy_D1p0Gy"),
    ]),
    ("120 s scan", BLUE, "s", [
        (120, "del120s_ac300s_1Gy_D1p0Gy_t120_240"),
        (180, "del180s_ac300s_1Gy_D1p0Gy_t180_300"),
    ]),
]


def washed(case):
    base = os.path.join(config_out(SCENARIO, TOPOLOGY, SCANNER, CRYSTAL),
                        "statistical_procedure", case, "combined")
    f = glob.glob(os.path.join(base, "washed_N*.toml"))
    if not f:
        raise FileNotFoundError(f"no washed combine in {base}")
    with open(f[0], "rb") as fh:
        d = tomllib.load(fh)
    return d["corrected_sigma_R_mm"], d["relative_sigma_uncertainty"]


def main():
    fig, ax = plt.subplots(figsize=(6.6, 4.6), facecolor=SURFACE)
    style(ax)
    for lab, col, mk, cases in SERIES:
        t = np.array([s for s, _ in cases])
        vals = [washed(c) for _, c in cases]
        s = np.array([v for v, _ in vals])
        band = np.array([b for _, b in vals])
        ax.errorbar(t, s, yerr=s * band, fmt=mk, ls="none", ms=8, color=col,
                    capsize=3, label=lab)
    ax.set_xticks([120, 180, 300])
    ax.set_xlim(90, 330)
    ax.set_ylim(0.0, 0.38)
    ax.set_xlabel("acquisition start delay [s]", color=INK, fontsize=12)
    ax.set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    ax.set_title("BGO ring: range precision vs start delay and scan duration",
                 color=INK, fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=12, labelcolor=INK, loc="upper left")
    fig.tight_layout()

    d = os.path.join(REPO, "latex", "figs")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "statproc_delay_bgo.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
