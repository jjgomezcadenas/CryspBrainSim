#!/usr/bin/env python3
"""
plot_grogg_v2.py — the linear-endpoint (Grogg et al.) comparison for the note.
At a short in-room scan ([120,300] s), on the two 35 cm-AFOV arms (CsI R35/35
and BGO r40/35), it contrasts three range estimators computed on identical
events: the erfc R50 of this note, the raw linear x-intercept, and the linear
intercept with the published 7 mm smoothing. From drivers/sigma_r_v2.jl
(--tend 300), reading each arm's washout_v2/sigma_r_grogg_v2_t120_300.toml.

Left panel — nominal σ_R at 1 Gy by estimator (log axis; the raw intercept is
30–60× worse than R50, the smoothing recovers ~5–8× of that but stays 5–7×
above R50), both arms overlaid, points with the ±1/√(2(N−1)) band.

Right panel — the mechanism, on the CsI arm: the distribution of the linear
fit-region start (the "last distal maximum") over the 100 realisations. Raw it
scatters ~5 mm across many bins (the start-selection lottery that dominates the
raw σ_R); the 7 mm smoothing collapses it to a single bin, which is why the
smoothed σ_R drops — without adding any counts.

Writes grogg_v2.png into <scenario>/<topology>/comparison/figures/. Run:
  python3 tools/plot_grogg_v2.py
"""
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
TAG = "t120_300"

# arm -> (label, colour, marker, scanner, leaf)
ARMS = [
    ("CsI R35/35", RED, "s", "crysp_r35_35cm_csi_2x0", "csi_2X0"),
    ("BGO r40/35", BLUE, "o", "crysp_r40_35cm_bgo_2x0", "bgo_195k_2X0"),
]
# estimator -> (x position, label, sigma-field prefix)
ESTS = [
    (0, "$R_{50}$\n(erfc)", "erfc_nominal"),
    (1, "linear\nintercept", "grogg_nominal"),
    (2, "linear\n+7 mm", "grogg_sm7_nominal"),
]


def load(scanner, leaf):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, leaf),
                        "washout_v2", f"sigma_r_grogg_v2_{TAG}.toml")
    with open(path, "rb") as f:
        d = tomllib.load(f)
    return d["point"][0], 1.0 / math.sqrt(2 * (d["realizations"] - 1))


def main():
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.5), facecolor=SURFACE,
                                   gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ---- left: nominal σ_R by estimator, both arms, log axis ----
    style(axL)
    axL.set_yscale("log")
    for lab, col, mk, scanner, leaf in ARMS:
        p, band = load(scanner, leaf)
        xs = [x for x, _, _ in ESTS]
        ys = np.array([p[f"{pre}_sigma_R_mm"] for _, _, pre in ESTS])
        dx = -0.07 if col is RED else 0.07
        axL.errorbar(np.array(xs) + dx, ys, yerr=ys * band, fmt=mk, ms=8,
                     color=col, capsize=3, label=lab)
    axL.set_xticks([x for x, _, _ in ESTS])
    axL.set_xticklabels([lab for _, lab, _ in ESTS], color=INK, fontsize=11)
    axL.set_xlim(-0.5, 2.5)
    axL.set_ylabel("nominal $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    axL.legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    axL.set_title("Range precision by estimator (short scan, [120,300] s)",
                  color=INK, fontsize=11.5, loc="left")

    # ---- right: the start-selection lottery on the CsI arm ----
    style(axR)
    p, _ = load(*ARMS[0][3:])
    zf_raw = np.array(p["realizations_grogg_nominal_z_first_mm"])
    zf_sm = np.array(p["realizations_grogg_sm7_nominal_z_first_mm"])
    lo, hi = math.floor(zf_raw.min()) - 1, math.ceil(zf_raw.max()) + 2
    bins = np.arange(lo, hi, 1.5)
    axR.hist(zf_raw, bins=bins, color=RED, alpha=0.55,
             label=f"raw (std {zf_raw.std():.1f} mm)")
    axR.hist(zf_sm, bins=bins, color=INK, alpha=0.85,
             label=f"+7 mm smoothed (std {zf_sm.std():.1f} mm)")
    axR.set_xlabel("linear fit-region start $z$ [mm]", color=INK, fontsize=12)
    axR.set_ylabel("realisations", color=INK, fontsize=12)
    axR.legend(frameon=False, fontsize=10.5, labelcolor=INK, loc="upper right")
    axR.set_title("Start-selection jitter, CsI (smoothing collapses it)",
                  color=INK, fontsize=11.5, loc="left")

    fig.tight_layout()
    d = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "grogg_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
