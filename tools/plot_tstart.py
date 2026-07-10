#!/usr/bin/env python3
"""
plot_tstart.py — the acquisition-start trade-off across the two scanners:
Δ_R50 and σ_R as functions of the start time t_start, at the working
protocol (all events, no corrections).

Reads each scanner's ten_shards/results.toml (the t_start = 0 baseline) and
ten_shards/tstart_<T>.toml (drivers/ten_shards_tstart.jl +
ten_shards.py --t-start). Writes the comparison figure into the
topology-level comparison/ leaf (a cross-scanner artifact sits above the
per-scanner tiers). The two scanners' points are offset by a few seconds
so the error bars stay legible; the counting comparison lives in the text.

Run:  python3 tools/plot_tstart.py [--t-starts 60,120,180]
"""
import argparse
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, crystal_label, scenario_out
from fit_activity_profile import BLUE, GRIDC, INK, MUTED, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
# The two scanners of the comparison: (label shown, scanner dir, crystal arm,
# wall mm, colour).
SCANNERS = [("BGO 195 K", "crysp_ring_1m_bgo_2x0", "bgo_195k", 22.36, RED),
            ("CsI", "crysp_ring_1m_csi_2x0", "csi", 37.20, BLUE)]
LEVEL = "all_ev"          # the working protocol


def series(scanner, crystal, wall_mm, t_starts):
    cfg = config_out(SCENARIO, TOPOLOGY, scanner,
                     crystal_label(crystal, wall_mm))
    with open(os.path.join(cfg, "ten_shards", "results.toml"), "rb") as f:
        base = tomllib.load(f)[LEVEL]["erfc"]
    t, mean, std = [0.0], [base["delta_R50_mean_mm"]], [base["delta_R50_std_mm"]]
    for ts in t_starts:
        with open(os.path.join(cfg, "ten_shards", f"tstart_{ts:g}.toml"),
                  "rb") as f:
            v = tomllib.load(f)[LEVEL]
        t.append(ts)
        mean.append(v["delta_R50_mean_mm"])
        std.append(v["delta_R50_std_mm"])
    return np.array(t), np.array(mean), np.array(std)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--t-starts", default="60,120,180")
    args = p.parse_args()
    t_starts = [float(v) for v in args.t_starts.split(",")]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 7.0), facecolor=SURFACE,
                                 sharex=True, height_ratios=[1.4, 1])
    for a in (a1, a2):
        style(a)
    n = 10
    for jitter, (name, scanner, crystal, wall, colour) in zip((-3.0, 3.0),
                                                              SCANNERS):
        t, mean, std = series(scanner, crystal, wall, t_starts)
        a1.errorbar(t + jitter, mean, yerr=std, fmt="o", ms=6, color=colour,
                    elinewidth=1.4, capsize=4, label=name)
        a1.plot(t + jitter, mean, color=colour, lw=0.9, alpha=0.5)
        a2.errorbar(t + jitter, std, yerr=std / np.sqrt(2 * (n - 1)),
                    fmt="o", ms=6, color=colour, elinewidth=1.4, capsize=4,
                    label=name)
    a2.set_ylim(0.0, None)
    a1.set_ylabel("$\\Delta R$ [mm]", color=INK, fontsize=13)
    a1.legend(frameon=False, fontsize=12, labelcolor=INK, loc="best")
    a2.set_ylabel("$\\sigma_R$ [mm]", color=INK, fontsize=13)
    a2.set_xlabel("acquisition start $t_{start}$ [s]", color=INK, fontsize=13)
    fig.tight_layout()

    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison",
                       "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "tstart_r50.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
