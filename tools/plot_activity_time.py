#!/usr/bin/env python3
"""
plot_activity_time.py — the induced-activity time profile (note Fig. 8,
fig:decays): the activity of the produced beta+ emitters as a function of time
after the beam turns off, with the acquisition window marked. This is the
committed source for latex/figs/activity_time.png; the figure was previously an
orphan PNG with no producer.

Inputs, all committed / on disk (no re-simulation):
  config/washout_brain.toml                     half-lives + irradiation/acquisition timing
  PtCryspProds/<scenario>/truth/
      activity_profile_fast.csv                 per-isotope decay counts D_i in [t1, t2] at 1 Gy
      run_meta.csv                              phantom material + dose, for the title

Physics. Production runs at a constant rate over the irradiation [0, t_irr];
each species then decays with its physical constant lambda_i = ln2 / T_half_i.
The truth CSV gives D_i, the number of decays of species i inside the
acquisition window [t1, t2] (t measured from the start of irradiation). With
build-up folded in, the nuclei present at beam-off are

  N_i(t_irr) = D_i / ( exp(-lambda_i (t1 - t_irr)) - exp(-lambda_i (t2 - t_irr)) ),

and the activity at a time tau after beam-off is

  A_i(tau) = lambda_i N_i(t_irr) exp(-lambda_i tau)      [Bq -> MBq].

The acquisition window opens at tau = t1 - t_irr and closes at tau = t2 - t_irr
(i.e. tau measured from beam-off; the old orphan PNG mislabelled this by t_irr).

C-10 (T_half = 19.3 s) is spent within ~1 min of beam-off, entirely before the
window opens; it is excluded from the curves and the total for a clean plot of
the acquisition-relevant emitters (matching the original figure).

Run:  python3 tools/plot_activity_time.py [--scenario-dir <.../PtCryspProds/<scenario>>]
Writes out/<scenario>/truth/figures/activity_time.png (override with --out-dir).
"""
import argparse
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from crysp_paths import truth_out  # noqa: E402
from fit_activity_profile import BLUE, INK, MUTED, RED, SURFACE, style  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_SCENARIO = os.path.join(os.path.dirname(REPO), "PtCryspProds",
                                "uniform_headep_sobp_1e8")
LN2 = np.log(2.0)
C10_MIN_HALFLIFE_S = 30.0        # exclude emitters spent before the window opens
GRID_SHADE = "#eeeeea"           # the acquisition-window shading

# One colour per long-lived emitter, echoing the note palette (O15 blue,
# C11 amber, N13 green, O14 red); the total is ink.
COLOUR = {"O15": BLUE, "C11": "#eda100", "N13": "#1baf7a", "O14": RED}


def read_row(path):
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        val = f.readline().strip().split(",")
    return dict(zip(hdr, val))


def read_meta(truth_dir):
    """Phantom material from run_meta; the D_i normalisation dose from the
    activity-profile meta (the counts are 'expected at dose_Gy', 1 Gy)."""
    material = read_row(os.path.join(truth_dir, "run_meta.csv")).get(
        "phantom_material", "brain")
    dose = float(read_row(os.path.join(truth_dir, "activity_profile_fast_meta.csv"))
                 .get("dose_Gy", 1.0))
    return material, dose


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scenario-dir", default=DEFAULT_SCENARIO,
                   help="PtCryspProds/<scenario> (reads truth/)")
    p.add_argument("--out-dir", default=None,
                   help="default: out/<scenario>/truth/figures")
    p.add_argument("--tmax-min", type=float, default=45.0)
    args = p.parse_args()

    scenario = os.path.basename(os.path.normpath(args.scenario_dir))
    truth_dir = os.path.join(args.scenario_dir, "truth")
    out_dir = args.out_dir or os.path.join(truth_out(scenario), "figures")

    wp = tomllib.load(open(os.path.join(REPO, "config", "washout_brain.toml"), "rb"))
    isos = wp["physical"]["isotopes"]
    thalf = dict(zip(isos, wp["physical"]["T_half_s"]))
    t_irr, t1, t2 = (wp["timing"]["t_irr_s"], wp["timing"]["t1_s"], wp["timing"]["t2_s"])

    tab = np.genfromtxt(os.path.join(truth_dir, "activity_profile_fast.csv"),
                        delimiter=",", names=True)
    decays = {iso: float(tab[iso].sum()) for iso in isos}   # D_i in [t1, t2] at 1 Gy
    material, dose = read_meta(truth_dir)

    tau = np.linspace(0.0, args.tmax_min * 60.0, 1200)      # s after beam-off
    shown = [iso for iso in isos if thalf[iso] >= C10_MIN_HALFLIFE_S]
    activity = {}
    for iso in shown:
        lam = LN2 / thalf[iso]
        n_beam_off = decays[iso] / (np.exp(-lam * (t1 - t_irr))
                                    - np.exp(-lam * (t2 - t_irr)))
        activity[iso] = lam * n_beam_off * np.exp(-lam * tau) / 1.0e6   # MBq
    total = np.sum([activity[iso] for iso in shown], axis=0)

    tau_min = tau / 60.0
    win_lo, win_hi = (t1 - t_irr) / 60.0, (t2 - t_irr) / 60.0   # window in min after beam-off

    fig, ax = plt.subplots(figsize=(9.0, 5.6), facecolor=SURFACE)
    style(ax)
    ax.axvspan(win_lo, win_hi, color=GRID_SHADE, lw=0, label="acquisition")
    ax.axvline(win_lo, color=MUTED, lw=1.0, ls="--")
    for iso in shown:
        ax.plot(tau_min, activity[iso], color=COLOUR[iso], lw=1.8, label=iso)
    ax.plot(tau_min, total, color=INK, lw=2.2, label="total")

    ax.set_yscale("log")
    ax.set_xlim(0.0, args.tmax_min)
    ax.set_ylim(1e-3, 1.2 * float(total.max()))
    ax.set_xlabel("time after beam-off [min]", color=INK, fontsize=12)
    ax.set_ylabel("activity [MBq]", color=INK, fontsize=12)
    ax.set_title(f"Induced activity ({material}, {dose:g} Gy)",
                 color=INK, fontsize=12)
    # isotopes, then total, then the acquisition band last (as in the note figure)
    h, lab = ax.get_legend_handles_labels()
    order = [lab.index(k) for k in shown + ["total", "acquisition"] if k in lab]
    ax.legend([h[i] for i in order], [lab[i] for i in order],
              frameon=False, fontsize=10, labelcolor=INK, loc="upper right")

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "activity_time.png")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")
    for iso in shown:
        print(f"  {iso}: A(beam-off) = {activity[iso][0]:.4f} MBq  "
              f"(T1/2 {thalf[iso]:.1f} s, D {decays[iso]:.3g})")
    print(f"  total A(beam-off) = {total[0]:.4f} MBq; "
          f"window [{win_lo:.0f}, {win_hi:.0f}] min after beam-off")


if __name__ == "__main__":
    main()
