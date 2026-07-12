#!/usr/bin/env python3
"""
plot_sigma_r_per_isotope.py — per-isotope range precision vs acquisition start,
from drivers/sigma_r_per_isotope.jl, testing the positron-range hypothesis.

Globs each start time's sigma_r_per_isotope_t{tstart}.toml (nominal) and
sigma_r_per_isotope_washed_t{tstart}.toml (washed) for O-15 and C-11, and reads
the combined (all-isotope) σ_R from the washout nominal/washed tomls in the same
directory (scaled to 1 Gy). Two panels:

  top:    σ_R vs t_start — nominal (solid) and washed (dashed) for O-15, C-11,
          and the combined mix, with ±1/√(2(N−1)) bands.
  bottom: washout inflation (washed/nominal) vs t_start for the three.

Reading: O-15 stays better than C-11 (its counts); per-count precision shows no
stable positron-range penalty on O-15; washout inflates both ≈ 1/√survival with
no isotope selectivity. Ring 1 m CsI.

Writes out/<scenario>/<topology>/comparison/figures/sigma_r_per_isotope.png.
Run:  python3 tools/plot_sigma_r_per_isotope.py
"""
import glob
import math
import os
import re
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, scenario_out
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
SCANNER, LEAF = "crysp_ring_1m_csi_2x0", "csi_2X0"
SPECIES = [("O15", "$^{15}$O", BLUE), ("C11", "$^{11}$C", RED),
           ("combined", "combined", INK)]
NREAL = 100        # per-isotope realizations (for the band)


def wdir():
    return os.path.join(config_out(SCENARIO, TOPOLOGY, SCANNER, LEAF), "washout")


def per_isotope(washed):
    """t_start -> {iso: sigma_R} from the sigma_r_per_isotope[_washed]_t*.toml."""
    stem = "sigma_r_per_isotope_washed_t" if washed else "sigma_r_per_isotope_t"
    out = {}
    for f in glob.glob(os.path.join(wdir(), stem + "*.toml")):
        m = re.search(r"_t(\d+)\.toml$", f)
        if not m:
            continue
        t = float(m.group(1))
        d = tomllib.load(open(f, "rb"))
        out[t] = {p["isotope"]: p["sigma_R_mm"] for p in d["point"]}
    return out


def combined():
    """t_start -> (nominal_1Gy, washed_1Gy), from the washout tomls, highest dose."""
    best = {}
    for f in glob.glob(os.path.join(wdir(), "*.toml")):
        d = tomllib.load(open(f, "rb"))
        dose = d.get("dose_Gy")
        if not dose or "point" not in d:
            continue
        for p in d["point"]:
            if "nominal_sigma_R_mm" not in p:
                continue
            t = p["t_start_s"]
            s = math.sqrt(dose)
            rec = (p["nominal_sigma_R_mm"] * s,
                   p.get("washed_sigma_R_mm", float("nan")) * s, dose)
            if t not in best or dose > best[t][2]:
                best[t] = rec
    return best


def main():
    nom_i, wsh_i, comb = per_isotope(False), per_isotope(True), combined()

    def series(key, washed):
        """(t array, sigma array) for a species (or 'combined')."""
        if key == "combined":
            ts = sorted(comb)
            idx = 1 if washed else 0
            vals = [comb[t][idx] for t in ts]
        else:
            src = wsh_i if washed else nom_i
            ts = sorted(t for t in src if key in src[t])
            vals = [src[t][key] for t in ts]
        ts = [t for t, v in zip(ts, vals) if v == v]        # drop NaN
        vals = [v for v in vals if v == v]
        return np.array(ts), np.array(vals)

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 8.0), facecolor=SURFACE,
                                 sharex=True, height_ratios=[1.4, 1])
    for a in (a1, a2):
        style(a)

    for key, label, colour in SPECIES:
        tn, sn = series(key, False)
        a1.errorbar(tn, sn, yerr=sn / math.sqrt(2 * (NREAL - 1)), fmt="-o", ms=6,
                    lw=1.4, color=colour, capsize=3, label=f"{label} nominal")
        tw, sw = series(key, True)
        if len(tw):
            a1.errorbar(tw, sw, yerr=sw / math.sqrt(2 * (NREAL - 1)), fmt="--s",
                        ms=6, lw=1.2, color=colour, mfc="none", capsize=3,
                        label=f"{label} washed")
            # inflation on the shared t_start points
            sn_at = {t: v for t, v in zip(tn, sn)}
            ti = [t for t in tw if t in sn_at]
            infl = [sw[list(tw).index(t)] / sn_at[t] for t in ti]
            a2.plot(ti, infl, "-o", ms=6, lw=1.4, color=colour, label=label)

    a1.set_ylabel("$\\sigma_R$ at 1 Gy [mm]", color=INK, fontsize=13)
    a1.set_ylim(0.0, None)
    a1.legend(frameon=False, fontsize=9.5, labelcolor=INK, loc="upper left", ncol=3)
    a1.set_title("Per-isotope range precision vs acquisition start "
                 "(ring 1 m CsI, 1 Gy)", color=INK, fontsize=11, loc="left")
    a2.axhline(1.0, color=MUTED, lw=0.9, ls=":")
    a2.set_ylabel("washout inflation", color=INK, fontsize=13)
    a2.set_xlabel("acquisition start $t_{start}$ [s]", color=INK, fontsize=13)
    a2.legend(frameon=False, fontsize=10, labelcolor=INK, loc="best")
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "sigma_r_per_isotope.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
