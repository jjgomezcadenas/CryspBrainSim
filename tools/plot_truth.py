#!/usr/bin/env python3
"""
plot_truth.py — the four starting-point figures from a scenario's truth/
bundle (validation ladder rung 1; the figure companion of
src/characterize.jl):

  activity.png       per-isotope + total activity(z)
  depth_dose.png     dose(z) with dose-R80 marked
  dose_activity.png  normalized dose and activity overlaid; R80, R50, offset
  sobp_plateau.png   dose zoomed to the SOBP plateau / target box

Run:  python3 tools/plot_truth.py --scenario-dir <.../PtCryspProds/<scenario>>
Writes the PNGs under out/characterize/figures/ (override with --out-dir).
"""
import argparse
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SCENARIO = os.path.join(os.path.dirname(REPO), "PtCryspProds",
                                "uniform_headep_sobp_1e8")

# dataviz reference palette, categorical slots 1–5 (CVD-validated order) for
# the five isotopes; ink/muted for total, dose and annotations.
SERIES = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7"]
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
C_DOSE, C_ACT = "#e34948", "#2a78d6"  # dose vs activity in the overlay
SURFACE = "#fcfcfb"


def read_table(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    cols = {}
    for j, name in enumerate(header):
        vals = [r[j] for r in data]
        try:
            cols[name.strip()] = np.array([float(v) for v in vals])
        except ValueError:
            cols[name.strip()] = np.array(vals)
    return cols


def distal_crossing(z, y, level, reference=None):
    """Last downward crossing of level*reference, linearly interpolated —
    the same reading as CryspBrainSim.distal_crossing."""
    thr = level * (reference if reference is not None else y.max())
    for i in range(len(y) - 2, -1, -1):
        if y[i] >= thr > y[i + 1]:
            return z[i] + (thr - y[i]) * (z[i + 1] - z[i]) / (y[i + 1] - y[i])
    return float("nan")


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def new_fig(title):
    fig, ax = plt.subplots(figsize=(8, 4.6), facecolor=SURFACE)
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    style(ax)
    return fig, ax


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name)
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario-dir", default=DEFAULT_SCENARIO)
    p.add_argument("--budget", default="fast")
    p.add_argument("--out-dir", default=os.path.join(REPO, "out", "characterize",
                                                     "figures"))
    args = p.parse_args()

    tdir = os.path.join(args.scenario_dir, "truth")
    scen = os.path.basename(os.path.abspath(args.scenario_dir))
    dose_t = read_table(os.path.join(tdir, "depth_dose.csv"))
    act_t = read_table(os.path.join(tdir, f"activity_profile_{args.budget}.csv"))
    run = read_table(os.path.join(tdir, "run_meta.csv"))
    z = dose_t["z_mm"]
    dose = dose_t["dose_core_Gy"]
    total = act_t["total"]
    isotopes = [k for k in act_t if k not in ("z_mm", "total")]

    R80 = distal_crossing(z, dose, 0.8)
    R50 = distal_crossing(act_t["z_mm"], total, 0.5)
    offset = R50 - R80
    print(f"dose-R80 = {R80:.3f} mm | activity-R50 = {R50:.3f} mm | "
          f"offset = {offset:.3f} mm")

    # 1. activity: per-isotope + total
    fig, ax = new_fig(f"{scen} — β⁺ activity(z), budget {args.budget}")
    for i, iso in enumerate(isotopes):
        ax.plot(act_t["z_mm"], act_t[iso], lw=1.6, color=SERIES[i % len(SERIES)],
                label=iso)
    ax.plot(act_t["z_mm"], total, lw=2.2, color=INK, label="total")
    ax.axvline(R50, color=MUTED, ls=":", lw=1.2)
    ax.text(R50, ax.get_ylim()[1] * 0.95, " activity-R50", color=MUTED, fontsize=8,
            ha="left", va="top")
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("expected decays / bin", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    save(fig, args.out_dir, "activity.png")

    # 2. depth dose with R80
    fig, ax = new_fig(f"{scen} — depth dose (core)")
    ax.plot(z, dose, lw=2, color=C_DOSE)
    ax.axvline(R80, color=MUTED, ls=":", lw=1.2)
    ax.text(R80, ax.get_ylim()[1] * 0.95, f" dose-R80 = {R80:.2f} mm", color=MUTED,
            fontsize=8, ha="left", va="top")
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("dose [Gy]", color=INK)
    save(fig, args.out_dir, "depth_dose.png")

    # 3. normalized overlay: the locked reference in one picture
    fig, ax = new_fig(f"{scen} — dose vs activity (normalized), "
                      f"offset = {offset:.2f} mm")
    ax.plot(z, dose / dose.max(), lw=2, color=C_DOSE, label="dose / max")
    ax.plot(act_t["z_mm"], total / total.max(), lw=2, color=C_ACT,
            label="activity / max")
    ax.axvline(R80, color=C_DOSE, ls=":", lw=1.2)
    ax.axvline(R50, color=C_ACT, ls=":", lw=1.2)
    ax.annotate(f"R50 {R50:.2f}", (R50, 0.5), textcoords="offset points",
                xytext=(-6, 4), ha="right", color=C_ACT, fontsize=8)
    ax.annotate(f"R80 {R80:.2f}", (R80, 0.8), textcoords="offset points",
                xytext=(6, 4), ha="left", color=C_DOSE, fontsize=8)
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("normalized", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    save(fig, args.out_dir, "dose_activity.png")

    # 4. SOBP plateau zoom over the target box (depths → z via the phantom's
    #    distal-edge convention already baked into the z-frame)
    fig, ax = new_fig(f"{scen} — SOBP plateau / target box")
    tgt = float(run["target_dose_Gy"][0])
    prox, dist = float(run["target_prox_depth_mm"][0]), float(run["target_dist_depth_mm"][0])
    # The target box in the z-frame ends at the nominal distal edge ≈ R80; its
    # length is (dist − prox).
    z_hi = R80
    z_lo = R80 - (dist - prox)
    ax.plot(z, dose, lw=2, color=C_DOSE)
    ax.axvspan(z_lo, z_hi, color=GRIDC, alpha=0.5, lw=0)
    ax.axhline(tgt, color=MUTED, ls="--", lw=1, label=f"target dose {tgt:.3g} Gy")
    ax.set_xlim(z_lo - 25, z_hi + 20)
    sel = (z > z_lo - 25) & (z < z_hi + 20)
    ax.set_ylim(0, 1.25 * dose[sel].max())
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("dose [Gy]", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    save(fig, args.out_dir, "sobp_plateau.png")


if __name__ == "__main__":
    main()
