#!/usr/bin/env python3
"""
plot_statproc_washed_grid.py — washed range precision sigma_R vs scanner AFOV
for the two crystals, from the statistical-procedure jobs (bounded erfc fit,
finite-pool-corrected, N=100). One point per scanner size (TBP / LAFOV / CAFOV),
BGO and CsI as two series with the +/-1/sqrt(2(N-1)) band. del120, 1 Gy.

Reads each arm's
  <scanner>/<crystal>/statistical_procedure/del120s_ac300s_1Gy_D1p0Gy/
      combined/washed_N100.toml   (field corrected_sigma_R_mm)

Writes statproc_washed_grid.png into <scenario>/<topology>/comparison/figures/.
Run:  python3 tools/plot_statproc_washed_grid.py
"""
import glob
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, scenario_out
from fit_activity_profile import BLUE, INK, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
CASE = "del120s_ac300s_1Gy_D1p0Gy"

# size -> (x, tick label); largest -> smallest AFOV
SIZES = [(0, "TBP\n(100 cm)"), (1, "LAFOV\n(50 cm)"), (2, "CAFOV\n(35 cm)")]
# crystal -> (label, colour, marker, [(scanner, crystal-label) per size])
ARMS = [
    ("BGO", BLUE, "o", [
        ("crysp_ring_1m_bgo_2x0", "bgo_195k_2X0"),
        ("crysp_r40_50cm_bgo_2x0", "bgo_195k_2X0"),
        ("crysp_r40_35cm_bgo_2x0", "bgo_195k_2X0"),
    ]),
    ("CsI", RED, "s", [
        ("crysp_ring_1m_csi_2x0", "csi_2X0"),
        ("crysp_r35_50cm_csi_2x0", "csi_2X0"),
        ("crysp_r35_35cm_csi_2x0", "csi_2X0"),
    ]),
]


def washed(scanner, crystal):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, crystal),
                        "statistical_procedure", CASE, "combined")
    f = glob.glob(os.path.join(path, "washed_N*.toml"))
    if not f:
        raise FileNotFoundError(f"no washed combine in {path}")
    with open(f[0], "rb") as fh:
        d = tomllib.load(fh)
    return d["corrected_sigma_R_mm"], d["relative_sigma_uncertainty"]


def main():
    fig, ax = plt.subplots(figsize=(6.6, 4.6), facecolor=SURFACE)
    style(ax)
    xs = [x for x, _ in SIZES]
    for lab, col, mk, arms in ARMS:
        s, band = zip(*(washed(sc, cr) for sc, cr in arms))
        s = np.array(s)
        ax.errorbar(np.array(xs), s, yerr=s * np.array(band), fmt=mk, ls="none",
                    ms=8, color=col, capsize=3, label=lab)
    ax.set_xticks(xs)
    ax.set_xticklabels([lab for _, lab in SIZES], color=INK, fontsize=11)
    ax.set_xlim(-0.4, 2.4)
    ax.set_ylim(0.0, 0.18)
    ax.set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    ax.set_title("Range precision vs scanner AFOV (del120, washout, N=100)",
                 color=INK, fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=12, labelcolor=INK, loc="upper left")
    fig.tight_layout()

    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "statproc_washed_grid.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
