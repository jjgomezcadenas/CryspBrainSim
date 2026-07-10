#!/usr/bin/env python3
"""
plot_washout.py — the isotope-washout result across the two 1 m scanners: what
washout does to the range measurement, separating the part that calibrates away
(the edge shift) from the parts that do not (its parameter systematic, and the
precision cost).

Left: R50 per shard, nominal vs washed, with the mean ± σ_R band — shows both
the calibrated shift and that σ_R barely grows. Right: σ_R nominal vs washed
per arm against the naive counting expectation σ_R/√f, with the truth-level
shift and its parameter systematic annotated.

Reads each arm's washout/sigma_r_washout.toml (drivers/washout_sigma_r.jl) and
the scenario-level washout/washout.toml (tools/washout.py, the systematic band).
Writes the comparison figure into the topology comparison/ leaf.

Run:  python3 tools/plot_washout.py
"""
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out, crystal_label, scenario_out
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style

SCENARIO, TOPOLOGY = "uniform_headep_sobp_1e8", "closed"
ARMS = [("BGO 195 K", "crysp_ring_1m_bgo_2x0", "bgo_195k", 22.36, RED),
        ("CsI", "crysp_ring_1m_csi_2x0", "csi", 37.20, BLUE)]


def load(scanner, crystal, wall):
    cfg = config_out(SCENARIO, TOPOLOGY, scanner, crystal_label(crystal, wall))
    with open(os.path.join(cfg, "washout", "sigma_r_washout.toml"), "rb") as f:
        return tomllib.load(f)


def main():
    with open(os.path.join(scenario_out(SCENARIO), "washout", "washout.toml"),
              "rb") as f:
        truth = tomllib.load(f)
    dR_truth = truth["delta_R50_washout_mm"]
    dR_sys = truth["delta_R50_washout_sys_mm"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 5.0), facecolor=SURFACE)
    for a in (a1, a2):
        style(a)

    for j, (name, scanner, crystal, wall, colour) in enumerate(ARMS):
        d = load(scanner, crystal, wall)
        nom = d["nominal"]["R50_per_shard_mm"]
        wo = d["washed"]["R50_per_shard_mm"]
        x0, x1 = 2 * j, 2 * j + 1
        # left: per-shard R50, nominal vs washed, mean ± σ_R
        for x, vals, mk in ((x0, nom, "nominal"), (x1, wo, "washed")):
            a1.plot(np.full(len(vals), x), vals, "o", ms=5, mfc="none",
                    mec=INK, mew=1.0, alpha=0.6)
            a1.errorbar([x], [np.mean(vals)], yerr=[np.std(vals, ddof=1)],
                        fmt="_", ms=22, color=colour, elinewidth=1.8, capsize=5)
        # right: σ_R nominal, washed, and counting expectation
        sn, sw = d["nominal"]["sigma_R_mm"], d["washed"]["sigma_R_mm"]
        a2.plot([x0, x1], [sn, sw], "-o", color=colour, ms=7, lw=1.4, label=name)
        a2.plot([x1], [sn * d["counting_expectation"]], "x", ms=10,
                color=colour, mew=2.0,
                label="counting $\\sigma_R/\\sqrt{f}$" if j == 0 else None)

    a1.set_xticks([0, 1, 2, 3], ["BGO\nnom", "BGO\nwash", "CsI\nnom",
                                 "CsI\nwash"], color=INK, fontsize=11)
    a1.set_ylabel("$R_{50}$ [mm]", color=INK, fontsize=13)
    a1.set_title(f"Edge position: shift $+{dR_truth:.2f}$ mm calibrates away",
                 color=INK, fontsize=11, loc="left")

    a2.set_xticks([0.5, 2.5], ["BGO 195 K", "CsI"], color=INK, fontsize=11)
    a2.set_xlim(-0.5, 3.5)
    a2.set_ylabel("$\\sigma_R$ [mm]", color=INK, fontsize=13)
    a2.set_title("Precision: washout removes the noisiest counts",
                 color=INK, fontsize=11, loc="left")
    a2.legend(frameon=False, fontsize=11, labelcolor=INK, loc="best")
    a2.text(0.02, 0.02, f"parameter systematic on the shift: "
            f"$\\pm{dR_sys:.3f}$ mm  ($\\ll\\sigma_R$)", transform=a2.transAxes,
            color=MUTED, fontsize=10)

    fig.tight_layout()
    out = os.path.join(scenario_out(SCENARIO), TOPOLOGY, "comparison", "figures")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, "washout_sigma_r.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
