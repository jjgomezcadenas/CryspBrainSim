#!/usr/bin/env python3
"""
plot_washout_thinned.py — the isotope-washout σ_R inflation vs acquisition
start, from the thinned firm-up (drivers/washout_sigma_r.jl --thinned), for
the two 1 m arms. This supersedes the n=10 washout_sigma_r.png / washout_tstart.png,
which underestimated the effect ("free at t=0" was an n=10 outlier artifact).

The thinned runs were done step-by-step at dose-adaptive doses (0.2–1 Gy, chosen
to keep the count-starved washed corner in the stable-fit regime), each preserved
as washout/thinned*.toml. This globs them per arm, dedups by t_start preferring
the highest dose (the 1 Gy anchor where present — inflation is dose-independent,
so this only picks the least-scaled estimate), and plots:

  top: σ_R inflation (washed/nominal) vs t_start, both arms, with the
       ±1/√(2(N−1)) band — the dimensionless, dose-independent result.
  bottom: σ_R washed and nominal, each scaled to 1 Gy (× √dose), both arms.

Writes the figure + a consolidated washout_thinned_curve.toml into the
topology comparison/ leaf.

Run:  python3 tools/plot_washout_thinned.py
"""
import glob
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, crystal_label, scenario_out
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style, toml_dump

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
ARMS = [("BGO 195 K", "crysp_ring_1m_bgo_2x0", "bgo_195k", 22.36, RED),
        ("CsI", "crysp_ring_1m_csi_2x0", "csi", 37.20, BLUE)]


def arm_curve(scanner, crystal, wall):
    """Consolidate all thinned tomls for one arm into a t_start-keyed curve,
    preferring the highest-dose estimate at each t_start."""
    wdir = os.path.join(config_out(SCENARIO, TOPOLOGY, scanner,
                                   crystal_label(crystal, wall)), "washout")
    files = (glob.glob(os.path.join(wdir, "thinned*.toml"))
             + glob.glob(os.path.join(wdir, "sigma_r_washout_thinned.toml")))
    by_t = {}
    for f in files:
        d = tomllib.load(open(f, "rb"))
        dose, nreal = d["dose_Gy"], d["realizations"]
        for p in d["point"]:
            t = p["t_start_s"]
            rec = dict(dose=dose, nreal=nreal, infl=p["inflation"],
                       sig_nom=p["nominal_sigma_R_mm"],
                       sig_wo=p["washed_sigma_R_mm"])
            if t not in by_t or dose > by_t[t]["dose"]:   # prefer higher dose
                by_t[t] = rec
    ts = sorted(by_t)
    return (np.array(ts), by_t)


def main():
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 7.0), facecolor=SURFACE,
                                 sharex=True, height_ratios=[1.3, 1])
    for a in (a1, a2):
        style(a)
    tree = {"meta": {"method": "thinned firm-up, dose-adaptive, scaled to 1 Gy",
                     "note": "supersedes the n=10 washout figures"}}
    for name, scanner, crystal, wall, colour in ARMS:
        ts, by_t = arm_curve(scanner, crystal, wall)
        infl = np.array([by_t[t]["infl"] for t in ts])
        band = np.array([infl[i] / np.sqrt(2 * (by_t[t]["nreal"] - 1))
                         for i, t in enumerate(ts)])
        # scale σ_R to 1 Gy: σ ∝ 1/√dose
        s1 = np.array([np.sqrt(by_t[t]["dose"]) for t in ts])
        sig_wo = np.array([by_t[t]["sig_wo"] for t in ts]) * s1
        sig_nom = np.array([by_t[t]["sig_nom"] for t in ts]) * s1

        a1.errorbar(ts, infl, yerr=band, fmt="-o", color=colour, ms=6, lw=1.3,
                    capsize=4, label=name)
        a2.plot(ts, sig_wo, "-o", color=colour, ms=6, lw=1.4,
                label=f"{name} washed")
        a2.plot(ts, sig_nom, "--s", color=colour, ms=5, lw=1.0, mfc="none",
                label=f"{name} nominal")
        tree[crystal] = {"t_start_s": [float(t) for t in ts],
                         "inflation": [float(x) for x in infl],
                         "sigma_R_washed_1Gy_mm": [float(x) for x in sig_wo],
                         "sigma_R_nominal_1Gy_mm": [float(x) for x in sig_nom],
                         "dose_Gy": [by_t[t]["dose"] for t in ts]}

    a1.axhline(1.0, color=MUTED, lw=0.9, ls=":")
    a1.set_ylabel("σ$_R$ inflation (washed / nominal)", color=INK, fontsize=12)
    a1.set_ylim(0.9, 2.3)
    a1.set_title("Washout inflates $\\sigma_R$ ~1.5× — not free, roughly flat "
                 "vs start", color=INK, fontsize=11, loc="left")
    a1.legend(frameon=False, fontsize=11, labelcolor=INK, loc="upper left")
    a2.set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a2.set_xlabel("acquisition start $t_{start}$ [s]", color=INK, fontsize=13)
    a2.legend(frameon=False, fontsize=10, labelcolor=INK, loc="upper left", ncol=2)
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "washout_inflation.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")
    toml_dump(os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison",
                           "washout_thinned_curve.toml"), tree)


if __name__ == "__main__":
    main()
