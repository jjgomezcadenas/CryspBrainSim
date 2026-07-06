#!/usr/bin/env python3
"""
plot_sigma_r.py — figures for the range-precision σ_R (validation ladder
rungs 6–7). σ_R is the spread of the fitted range endpoint R across many
acquisitions; these figures show that spread and how it behaves.

Reads what the σ_R drivers write into out/sigma_r/:
  from_shards.{npz,toml}   drivers/sigma_r_at_dose.jl --from-shards
  at_dose_1Gy.{npz,toml}   drivers/sigma_r_at_dose.jl --realizations N
  sweep.{npz,toml}         drivers/sigma_r_sweep_dose.jl

Modes:
  (default)   per-shard R50 in both conventions, with mean and ±σ_R band
              → figures/from_shards.png
  --at-dose   the thinned realizations at the nominal dose, with the gate
              comparing their σ_R to the shard reference → figures/at_dose.png
  --hist      the R50 distribution of the thinned realizations, with a
              Gaussian(mean, σ_R) overlay → figures/r50_hist.png
  --sweep     σ_R vs dose (the curve) → figures/sweep.png

Run:  python3 tools/plot_sigma_r.py [--at-dose | --hist | --sweep]
"""
import os
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "out", "sigma_r")
NOMINAL = "at_dose_1Gy"  # the nominal-dose thinned run

BLUE, AQUA = "#2a78d6", "#1baf7a"
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
SURFACE = "#fcfcfb"


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def load(name):
    d = np.load(os.path.join(OUT, name + ".npz"))
    with open(os.path.join(OUT, name + ".toml"), "rb") as f:
        res = tomllib.load(f)
    return d, res


def save(fig, name):
    figdir = os.path.join(OUT, "figures")
    os.makedirs(figdir, exist_ok=True)
    path = os.path.join(figdir, name)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


# --- default: σ_R from the ten independent shards --------------------------
def plot_from_shards():
    d, res = load("from_shards")

    def panel(ax, vals, errs, stats, color, title):
        mean, sig = stats["mean_mm"], stats["sigma_mm"]
        ax.axhspan(mean - sig, mean + sig, color=GRIDC, alpha=0.6, lw=0)
        ax.axhline(mean, color=MUTED, ls="--", lw=1)
        if errs is not None:
            ax.errorbar(d["shards"], vals, yerr=errs, fmt="o", ms=5, color=color,
                        capsize=2, lw=1.4)
        else:
            ax.plot(d["shards"], vals, "o", ms=5, color=color)
        ax.set_title(f"{title}: mean {mean:.3f} mm, σ_R {sig:.3f} mm",
                     color=INK, fontsize=10, loc="left")
        ax.set_xlabel("shard", color=INK)
        ax.set_xticks(d["shards"])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)
    for a in (a1, a2):
        style(a)
    panel(a1, d["r50_fit"], d["z0_err"], res["sigma_R_fit"], BLUE,
          "erfc fit (bars: each fit's own error)")
    a1.set_ylabel("R50 [mm]", color=INK)
    panel(a2, d["r50_cross"], None, res["sigma_R_crossing"], AQUA,
          "windowed crossing")
    fig.suptitle(f"σ_R from the {res['n_shards']} independent shards at 1 Gy "
                 f"(grey band = ±σ_R)", color=INK, fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "from_shards.png")


# --- --at-dose: thinned realizations + the gate ----------------------------
def plot_at_dose():
    d, res = load(NOMINAL)
    ref = None
    if os.path.exists(os.path.join(OUT, "from_shards.toml")):
        with open(os.path.join(OUT, "from_shards.toml"), "rb") as f:
            ref = tomllib.load(f)

    fits = d["r50_fit"]
    idx = np.arange(1, fits.size + 1)
    st = res["sigma_R_fit"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURFACE,
                                 gridspec_kw={"width_ratios": [2.2, 1]})
    for a in (a1, a2):
        style(a)

    a1.axhspan(st["mean_mm"] - st["sigma_mm"], st["mean_mm"] + st["sigma_mm"],
               color=GRIDC, alpha=0.6, lw=0)
    a1.axhline(st["mean_mm"], color=MUTED, ls="--", lw=1)
    a1.errorbar(idx, fits, yerr=d["z0_err"], fmt="o", ms=4, color=BLUE,
                capsize=2, lw=1.2)
    a1.set_title(f"{fits.size} thinned realizations at 1 Gy (erfc fit): "
                 f"mean {st['mean_mm']:.3f} mm, σ_R {st['sigma_mm']:.3f} mm",
                 color=INK, fontsize=10, loc="left")
    a1.set_xlabel("realization", color=INK)
    a1.set_ylabel("R50 [mm]", color=INK)

    a2.bar(["thinned\n(N=%d)" % fits.size], [st["sigma_mm"]], color=BLUE, width=0.5,
           label=f"thinned {st['sigma_mm']:.3f} mm")
    title = f"gate: σ_R = {st['sigma_mm']:.3f} mm"
    if ref is not None:
        s0 = ref["sigma_R_fit"]["sigma_mm"]
        band = s0 * np.sqrt(1 / (2 * (fits.size - 1)) + 1 / (2 * (ref["n_shards"] - 1)))
        a2.axhspan(s0 - 2 * band, s0 + 2 * band, color=AQUA, alpha=0.25, lw=0,
                   label="shard reference ±2σ")
        a2.axhline(s0, color=AQUA, lw=1.8, label=f"reference {s0:.3f} mm")
        g = res.get("gate", {})
        verdict = "AGREES" if g.get("pass") else "check"
        title += f"  (ratio {g.get('ratio', float('nan')):.2f} → {verdict})"
        a2.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
        a2.set_ylim(0, max(st["sigma_mm"], s0) * 1.6)
    a2.set_title(title, color=INK, fontsize=9.5, loc="left")
    a2.set_ylabel("σ_R [mm]", color=INK)

    fig.suptitle("Thinned σ_R at 1 Gy vs the shard reference (rung-6 gate)",
                 color=INK, fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "at_dose.png")


# --- --hist: the R50 distribution ------------------------------------------
def plot_hist():
    d, res = load(NOMINAL)

    def panel(ax, vals, stats, color, title):
        mean, sig = stats["mean_mm"], stats["sigma_mm"]
        ax.hist(vals, bins=12, color=color, alpha=0.75, edgecolor=SURFACE, lw=1.2)
        lo, hi = vals.min(), vals.max()
        pad = 0.25 * (hi - lo)
        xs = np.linspace(lo - pad, hi + pad, 300)
        binw = (hi - lo) / 12
        pdf = np.exp(-0.5 * ((xs - mean) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        ax.plot(xs, vals.size * binw * pdf, color=INK, lw=1.6,
                label=f"Gauss(μ={mean:.3f}, σ={sig:.3f})")
        ax.axvline(mean, color=MUTED, ls="--", lw=1)
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        ax.set_xlabel("R50 [mm]", color=INK)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK)

    n = d["r50_fit"].size
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)
    for a in (a1, a2):
        style(a)
    panel(a1, d["r50_fit"], res["sigma_R_fit"], BLUE, "erfc fit")
    a1.set_ylabel("realizations", color=INK)
    panel(a2, d["r50_cross"], res["sigma_R_crossing"], AQUA, "windowed crossing")
    fig.suptitle(f"R50 distribution — {n} thinned realizations at 1 Gy",
                 color=INK, fontsize=11, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, "r50_hist.png")


# --- --sweep: σ_R vs dose ---------------------------------------------------
def plot_sweep():
    d, res = load("sweep")
    dose = d["dose_Gy"]
    fig, ax = plt.subplots(figsize=(8, 4.8), facecolor=SURFACE)
    style(ax)
    ax.plot(dose, d["sigma_fit_mm"], "o-", ms=6, color=BLUE, lw=1.8,
            label="erfc fit")
    ax.plot(dose, d["sigma_crossing_mm"], "s--", ms=5, color=AQUA, lw=1.5,
            label="windowed crossing")
    # 1/sqrt(N) reference anchored on the top-dose fit point.
    top = np.argmax(dose)
    ax.plot(dose, d["sigma_fit_mm"][top] * np.sqrt(dose[top] / dose),
            color=MUTED, ls=":", lw=1.2, label="1/√N (anchored at top dose)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("dose [Gy]", color=INK)
    ax.set_ylabel("σ_R [mm]", color=INK)
    ax.set_title("σ_R vs dose — crysp_ring_1m, BGO (rung 7)",
                 color=INK, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()
    save(fig, "sweep.png")


def main():
    if "--at-dose" in sys.argv:
        plot_at_dose()
    elif "--hist" in sys.argv:
        plot_hist()
    elif "--sweep" in sys.argv:
        plot_sweep()
    else:
        plot_from_shards()


if __name__ == "__main__":
    main()
