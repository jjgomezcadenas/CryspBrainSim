#!/usr/bin/env python3
"""
plot_sigma_r_vs_niter.py — washed σ_R as a function of MLEM iterations, per
scanner configuration and start time, from drivers/washout_niter_scan.jl.

Tests whether the fixed 50-iteration operating point is what makes σ_R differ
across geometries: if the σ_R(niter) curves have their minima at different niter
(or cross), the frozen niter is imposing a config-dependent σ_R rather than the
counts. Two rows — top: washed σ_R(niter) (precision), bottom: mean washed R50
(niter) (the bias/convergence, so a low-σ_R low-niter point is not mistaken for
optimal when the edge is still under-converged). Columns are the start times;
each curve carries the ±1/√(2(N−1)) band, the frozen niter is marked, and each
σ_R curve's minimum is flagged.

Reads each config's washout/sigma_r_vs_niter.toml. Add configs to CONFIGS.
Writes out/<scenario>/<topology>/comparison/figures/sigma_r_vs_niter.png.
Run:  python3 tools/plot_sigma_r_vs_niter.py
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
TIMES = [120.0, 180.0]

# (label, scanner, leaf, colour) — extend as more configs are scanned
CONFIGS = [
    ("R35/50 CsI", "crysp_r35_50cm_csi_2x0", "csi_2X0", BLUE),
    ("R35/35 CsI", "crysp_r35_35cm_csi_2x0", "csi_2X0", RED),
]


def load(scanner, leaf):
    path = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner, leaf),
                        "washout", "sigma_r_vs_niter.toml")
    if not os.path.exists(path):
        return None
    d = tomllib.load(open(path, "rb"))
    niter = np.array(d["niter"], float)
    n = d["realizations"]
    frozen = d["frozen_niter"]
    sig = {p["t_start_s"]: np.array(p["washed_sigma_R_mm"], float) for p in d["point"]}
    r50 = {p["t_start_s"]: np.array(p["washed_mean_R50_mm"], float) for p in d["point"]}
    return niter, n, frozen, sig, r50


def main():
    fig, axes = plt.subplots(2, len(TIMES), figsize=(11.5, 8.0),
                             facecolor=SURFACE, sharex=True)
    frozen_ref = None
    for j, t in enumerate(TIMES):
        a_sig, a_r50 = axes[0, j], axes[1, j]
        style(a_sig); style(a_r50)
        for label, scanner, leaf, colour in CONFIGS:
            got = load(scanner, leaf)
            if got is None:
                continue
            niter, n, frozen, sig, r50 = got
            frozen_ref = frozen
            if t not in sig:
                continue
            band = sig[t] / math.sqrt(2 * (n - 1))
            a_sig.errorbar(niter, sig[t], yerr=band, fmt="-o", ms=5, lw=1.3,
                           color=colour, capsize=3, label=label)
            k = int(np.nanargmin(sig[t]))
            a_sig.plot(niter[k], sig[t][k], "*", ms=15, color=colour, mec=INK,
                       mew=0.6, zorder=5)
            a_r50.plot(niter, r50[t], "-o", ms=5, lw=1.3, color=colour, label=label)
        for ax in (a_sig, a_r50):
            if frozen_ref is not None:
                ax.axvline(frozen_ref, color=MUTED, lw=1.0, ls=":")
        a_sig.text(frozen_ref, a_sig.get_ylim()[1], f" frozen={frozen_ref}",
                   color=MUTED, fontsize=9, va="top", ha="left")
        a_sig.set_title(f"$t_{{start}}$ = {t:.0f} s", color=INK, fontsize=12, loc="left")
        a_r50.set_xlabel("MLEM iterations", color=INK, fontsize=12)
    axes[0, 0].set_ylabel("washed $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    axes[0, 0].set_ylim(0.0, None)
    axes[1, 0].set_ylabel("mean washed $R_{50}$ [mm]", color=INK, fontsize=13)
    axes[0, -1].legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    fig.suptitle("Washed $\\sigma_R$ (top, ★=min) and mean $R_{50}$ (bottom) "
                 "vs MLEM iterations", color=INK, fontsize=12)
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "sigma_r_vs_niter.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
