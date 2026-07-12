#!/usr/bin/env python3
"""
plot_sigma_r_v2.py — the generation-2 σ_R figure for the note, from
drivers/sigma_r_v2.jl, over the three acquisition-scenario leaves
(del120/180/300, the delay axis on the irradiation-end clock). One figure,
two panels, points only (no connecting lines):

  top:    washout σ_R vs delay — nominal vs exact-washed (points, ±band).
  bottom: pure per-isotope σ_R vs delay — O15, C11, and the combined mix.

The inflation and per-count precision k live in the note's tables, not the
figure. Reads the v2 tomls under <cfg>/washout_v2/ (active arm) and writes
sigma_r_v2.png into <cfg>/washout_v2/figures/. Run:
  cp config/run_parameters_csi_v2.toml config/run_parameters.toml
  python3 tools/plot_sigma_r_v2.py
"""
import math
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import active_config
from fit_activity_profile import BLUE, INK, RED, SURFACE, style


def load(name):
    path = os.path.join(active_config().cfg_dir, "washout_v2", name)
    with open(path, "rb") as f:
        return tomllib.load(f)


def main():
    dw = load("sigma_r_washout_v2.toml")
    di = load("sigma_r_per_isotope_v2.toml")
    band = 1.0 / math.sqrt(2 * (dw["realizations"] - 1))

    wp = sorted(dw["point"], key=lambda p: p["t_del_s"])
    t = np.array([p["t_del_s"] for p in wp])
    nom = np.array([p["nominal_sigma_R_mm"] for p in wp])
    wsh = np.array([p["washed_sigma_R_mm"] for p in wp])

    def iso_sigma(name):
        pts = sorted((p for p in di["point"] if p["isotope"] == name),
                     key=lambda p: p["t_del_s"])
        return (np.array([p["t_del_s"] for p in pts]),
                np.array([p["sigma_R_mm"] for p in pts]))

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 7.4), facecolor=SURFACE,
                                 sharex=True)
    for a in (a1, a2):
        style(a)

    # top: washout σ_R (points only)
    a1.errorbar(t, nom, yerr=nom * band, fmt="o", ms=7, color=BLUE, capsize=3,
                label="nominal (no washout)")
    a1.errorbar(t, wsh, yerr=wsh * band, fmt="s", ms=7, color=RED, mfc="none",
                capsize=3, label="washed (exact $g_i$ keep)")
    a1.set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a1.set_ylim(0.0, None)
    a1.legend(frameon=False, fontsize=10, labelcolor=INK, loc="upper left")
    a1.set_title("Range precision vs acquisition delay (ring 1 m CsI, v2, 1 Gy)",
                 color=INK, fontsize=11, loc="left")

    # bottom: per-isotope σ_R (points only)
    for name, label, col in (("O15", "$^{15}$O", BLUE), ("C11", "$^{11}$C", RED)):
        ti, si = iso_sigma(name)
        a2.errorbar(ti, si, yerr=si * band, fmt="o", ms=7, color=col, capsize=3,
                    label=label)
    a2.errorbar(t, nom, yerr=nom * band, fmt="D", ms=6, color=INK, mfc="none",
                capsize=3, label="combined (all isotopes)")
    a2.set_ylabel("per-isotope $\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a2.set_ylim(0.0, None)
    a2.set_xlabel("acquisition delay $t_{\\rm del}$ [s]", color=INK, fontsize=13)
    a2.set_xticks(t)
    a2.set_xlim(t[0] - 30, t[-1] + 30)
    a2.legend(frameon=False, fontsize=10, labelcolor=INK, loc="upper left", ncol=3)
    fig.tight_layout()

    d = os.path.join(active_config().cfg_dir, "washout_v2", "figures")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, "sigma_r_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
