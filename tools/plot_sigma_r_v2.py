#!/usr/bin/env python3
"""
plot_sigma_r_v2.py — the generation-2 σ_R figures for the note, from
drivers/sigma_r_v2.jl, over the three acquisition-scenario leaves
(del120/180/300, the delay axis on the irradiation-end clock). Two figures:

  washout_v2.png     — top: nominal vs exact-washed σ_R vs delay (points, ±band);
                       bottom: washout inflation vs delay, against the counting
                       expectation 1/√survival (dotted).
  per_isotope_v2.png — top: pure per-isotope σ_R (O15, C11) and the combined mix
                       vs delay; bottom: per-count precision k = σ_R·√N — the
                       positron-range test (O15, longest range, is flat-lowest).

Reads the v2 tomls under <cfg>/washout_v2/ (active arm = ring 1 m CsI v2) and
writes both PNGs into <cfg>/washout_v2/figures/. Run:
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
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style


def load(name):
    ac = active_config()
    path = os.path.join(ac.cfg_dir, "washout_v2", name)
    with open(path, "rb") as f:
        return tomllib.load(f)


def figdir():
    d = os.path.join(active_config().cfg_dir, "washout_v2", "figures")
    os.makedirs(d, exist_ok=True)
    return d


def plot_washout():
    d = load("sigma_r_washout_v2.toml")
    band = 1.0 / math.sqrt(2 * (d["realizations"] - 1))
    pts = sorted(d["point"], key=lambda p: p["t_del_s"])
    t = np.array([p["t_del_s"] for p in pts])
    nom = np.array([p["nominal_sigma_R_mm"] for p in pts])
    wsh = np.array([p["washed_sigma_R_mm"] for p in pts])
    infl = np.array([p["inflation"] for p in pts])
    surv = np.array([p["washout_survival"] for p in pts])

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 6.6), facecolor=SURFACE,
                                 sharex=True, height_ratios=[1.5, 1])
    for a in (a1, a2):
        style(a)
    a1.errorbar(t, nom, yerr=nom * band, fmt="-o", ms=7, lw=1.4, color=BLUE,
                capsize=3, label="nominal (no washout)")
    a1.errorbar(t, wsh, yerr=wsh * band, fmt="--s", ms=7, lw=1.2, color=RED,
                mfc="none", capsize=3, label="washed (exact $g_i$ keep)")
    a1.set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a1.set_ylim(0.0, None)
    a1.legend(frameon=False, fontsize=10, labelcolor=INK, loc="upper left")
    a1.set_title("Washout $\\sigma_R$ cost vs acquisition delay (ring 1 m CsI, v2, 1 Gy)",
                 color=INK, fontsize=11, loc="left")

    a2.errorbar(t, infl, yerr=infl * band * math.sqrt(2), fmt="-o", ms=7, lw=1.4,
                color=INK, capsize=3, label="washed / nominal")
    a2.plot(t, 1.0 / np.sqrt(surv), ":", lw=1.2, color=MUTED,
            label="counting $1/\\sqrt{\\mathrm{survival}}$")
    a2.axhline(1.0, color=MUTED, lw=0.8, ls="-", alpha=0.4)
    a2.set_ylabel("inflation", color=INK, fontsize=13)
    a2.set_xlabel("acquisition delay $t_{\\rm del}$ [s]", color=INK, fontsize=13)
    a2.set_xticks(t)
    a2.legend(frameon=False, fontsize=9.5, labelcolor=INK, loc="best")
    fig.tight_layout()
    path = os.path.join(figdir(), "washout_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def plot_per_isotope():
    di = load("sigma_r_per_isotope_v2.toml")
    dw = load("sigma_r_washout_v2.toml")
    band = 1.0 / math.sqrt(2 * (di["realizations"] - 1))

    def iso_series(name):
        pts = sorted((p for p in di["point"] if p["isotope"] == name),
                     key=lambda p: p["t_del_s"])
        return (np.array([p["t_del_s"] for p in pts]),
                np.array([p["sigma_R_mm"] for p in pts]),
                np.array([p["mean_events"] for p in pts]))

    wpts = sorted(dw["point"], key=lambda p: p["t_del_s"])
    tc = np.array([p["t_del_s"] for p in wpts])
    comb = np.array([p["nominal_sigma_R_mm"] for p in wpts])

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 6.6), facecolor=SURFACE,
                                 sharex=True, height_ratios=[1.3, 1])
    for a in (a1, a2):
        style(a)
    for name, label, col in (("O15", "$^{15}$O", BLUE), ("C11", "$^{11}$C", RED)):
        t, s, nev = iso_series(name)
        a1.errorbar(t, s, yerr=s * band, fmt="-o", ms=7, lw=1.4, color=col,
                    capsize=3, label=label)
        a2.errorbar(t, s * np.sqrt(nev), yerr=s * np.sqrt(nev) * band, fmt="-o",
                    ms=7, lw=1.4, color=col, capsize=3, label=label)
    a1.errorbar(tc, comb, yerr=comb * band, fmt="--D", ms=6, lw=1.1, color=INK,
                mfc="none", capsize=3, label="combined (all isotopes)")
    a1.set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a1.set_ylim(0.0, None)
    a1.legend(frameon=False, fontsize=10, labelcolor=INK, loc="upper left", ncol=3)
    a1.set_title("Per-isotope range precision vs delay "
                 "(ring 1 m CsI, v2, exact selection, 1 Gy)",
                 color=INK, fontsize=11, loc="left")
    a2.set_ylabel("per-count $k=\\sigma_R\\sqrt{N}$", color=INK, fontsize=13)
    a2.set_xlabel("acquisition delay $t_{\\rm del}$ [s]", color=INK, fontsize=13)
    a2.set_xticks(tc)
    a2.legend(frameon=False, fontsize=10, labelcolor=INK, loc="best")
    fig.tight_layout()
    path = os.path.join(figdir(), "per_isotope_v2.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    plot_washout()
    plot_per_isotope()


if __name__ == "__main__":
    main()
