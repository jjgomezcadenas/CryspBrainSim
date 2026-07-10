#!/usr/bin/env python3
"""
ten_shards.py — endpoint stability over the ten stored shards: fit the
distal edge of every shard at each rung of the degradation ladder and
tabulate the activity−dose endpoint differences with their shard spread.
With --dose-sweep, the dose axis instead: the same fit on every thinned
realization (drivers/ten_shards_dose.jl), grouped by dose.

The rungs, in order of added effects (all whole-plane, free baseline,
frozen window recipe — the protocol settled in the fit lab):
  origins   the detected trues' TRUE annihilation positions, histogrammed
            straight from the shard file — detector acceptance only, no
            reconstruction
  recon     the MLEM reconstruction of the trues (drivers/one_shard.jl)
  all_ev    the MLEM reconstruction of ALL events — scatters and randoms
            in, uncorrected (drivers/one_shard.jl --all-uncorr)

Each rung is fitted with both edge models (erfc and sigmoid, free
baseline), giving per shard the midpoint R50 and the tangent endpoint Rp.
The references — the truth activity profile and the depth-dose curve —
are shard-independent and fitted once. The headline per rung × model:
  Δ_R50 = R50 − R50(dose),  Δ_Rp = Rp − Rp(dose)
as the 10-shard mean ± sem, with the shard-to-shard std alongside. The
model whose Δs move least across the rungs is the more robust endpoint
convention.

The fit machinery is imported from fit_activity_profile.py — one place
defines the models, windows, and profiles. No per-shard figures are
written here (inspect any single case with fit_activity_profile.py);
the two summary figures show every shard as a point.

The Fano test (--fano) is the model-free check of the origins error
assignment: per z bin, the variance of the count across the ten shards
divided by its mean. Independent same-condition acquisitions give exactly
Fano = 1 for Poisson counts — no fit model enters anywhere. A window
average well above 1 would mean super-Poisson bin errors; ≈ 1 pins the
origins fit's large χ²/ndf on the erfc shape, not the errors.

The dose sweep (--dose-sweep) answers the test-dose question: how precisely
does a single low-dose acquisition locate the distal edge? Per dose
(1.0 Gy = the ten shards; 0.5/0.2/0.1 Gy = seeded thins of each shard) every
realization is fitted, the group mean of Δ_R50 checks that the calibration
shift stays put, and the group std IS σ_R(dose) — compared against the
1/√dose scaling anchored at 1 Gy. erfc in the figures, both models in the
TOML.

Writes  <config>/ten_shards/results.toml
        <config>/ten_shards/figures/delta_r50.png, delta_rp.png
  and with --dose-sweep
        <config>/ten_shards/dose_sweep.toml
        <config>/ten_shards/figures/delta_r50_vs_dose.png
  and with --fano
        <config>/ten_shards/fano.toml
        <config>/ten_shards/figures/fano_origins.png
"""
import argparse
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import REPO
from fit_activity_profile import (
    ACTIVE, BLUE, CFG, GRIDC, INK, MUTED, RED, SURFACE, TRUTH,
    analyze, edge_window, profile_from_image, profile_from_origins, style,
    toml_dump)

MODELS = ("erfc", "sigmoid")
SHARDS = range(10)
RUNGS = ("origins", "recon", "all_ev")
# Publication wording: the rungs are "levels of the measurement chain".
RUNG_TITLE = {"origins": "detected origins",
              "recon": "reconstruction\n(trues)",
              "all_ev": "reconstruction\n(all events)"}
OUT = os.path.join(CFG, "ten_shards")


def shard_h5(i):
    return os.path.join(ACTIVE.products_leaf, f"lors_shard{i:03d}.h5")


def rung_profile(rung, shard, params, centre):
    """The whole-plane depth profile of one shard at one rung."""
    if rung == "origins":
        return profile_from_origins(shard_h5(shard), params["grid"], centre,
                                    None)
    tag = f"shard{shard:03d}" + ("_all_uncorr" if rung == "all_ev" else "")
    d = np.load(os.path.join(CFG, "one_shard", f"recon_{tag}.npz"))
    return profile_from_image(d["image"], params["grid"], centre, None)


def summarize(vals):
    v = np.asarray(vals, float)
    return {"mean": float(v.mean()), "std": float(v.std(ddof=1)),
            "sem": float(v.std(ddof=1) / np.sqrt(v.size))}


def plot_deltas(results, dose, truth_act, obs, path):
    """Δ_obs per level of the measurement chain: every run as a point, the
    mean ± sem as a bar (erfc — the paper's model; the sigmoid cross-check
    stays in the TOML). The truth-activity level is the dashed reference
    line; the y-range follows the data."""
    fig, ax = plt.subplots(figsize=(8.5, 5.2), facecolor=SURFACE)
    style(ax)
    mk = "erfc"
    tline = truth_act[mk][f"{obs}_mm"] - dose[mk][f"{obs}_mm"]
    ax.axhline(tline, color=BLUE, ls="--", lw=1.1,
               label="truth activity")
    vals = [tline]
    for ix, rung in enumerate(RUNGS):
        r = results[rung][mk]
        d = np.asarray(r[f"{obs}_mm"]) - dose[mk][f"{obs}_mm"]
        ax.plot(np.full(d.size, ix), d, "o", ms=5, mfc="none",
                mec=INK, mew=1.1, alpha=0.7,
                label="single runs" if ix == 0 else None)
        s = summarize(d)
        ax.errorbar([ix], [s["mean"]], yerr=[s["std"]], fmt="_", ms=20,
                    color=RED, elinewidth=1.8, capsize=5,
                    label="mean ± σ" if ix == 0 else None)
        vals.extend(d)
    lo, hi = min(vals), max(vals)
    pad = 0.10 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlim(-0.5, len(RUNGS) - 0.5)
    ax.set_xticks(range(len(RUNGS)),
                  [RUNG_TITLE[r] for r in RUNGS], color=INK, fontsize=12)
    ax.set_ylabel(f"$\\Delta R$ [mm]", color=INK, fontsize=13)
    ax.legend(frameon=False, fontsize=12, labelcolor=INK, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def setup():
    """The run parameters and the two shard-independent references: the
    truth activity profile (frozen window) and the depth dose (the same
    margins around its own edge)."""
    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        params = tomllib.load(f)
    centre = params["roi"]["centre_mm"]
    window = (params["window"]["z_lo_mm"], params["window"]["z_hi_mm"])
    margins = (params["window"]["proximal_margin_mm"],
               params["window"]["distal_margin_mm"])
    os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
    ta = np.genfromtxt(os.path.join(TRUTH, "activity_profile_fast.csv"),
                       delimiter=",", names=True)
    dd = np.genfromtxt(os.path.join(TRUTH, "depth_dose.csv"),
                       delimiter=",", names=True)
    zd, yd = np.asarray(dd["z_mm"], float), np.asarray(dd["dose_core_Gy"], float)
    truth_act = analyze("truth_activity", np.asarray(ta["z_mm"], float),
                        np.asarray(ta["total"], float), window, 0.0, True,
                        models=MODELS)
    dose = analyze("truth_dose", zd, yd, edge_window(zd, yd, margins), 0.0,
                   False, models=MODELS)
    return params, centre, window, truth_act, dose


def ladder(params, centre, window, truth_act, dose):
    # The ladder: every shard at every rung, both models.
    results = {r: {mk: {"R50_mm": [], "R50_err_mm": [], "Rp_mm": [],
                        "Rp_err_mm": [], "chi2_ndf": []}
                   for mk in MODELS} for r in RUNGS}
    for rung in RUNGS:
        for shard in SHARDS:
            z, prof = rung_profile(rung, shard, params, centre)
            res = analyze(f"{rung}_shard{shard:03d}", z, prof, window, 0.0,
                          True, models=MODELS)
            for mk in MODELS:
                e = res[mk]
                r = results[rung][mk]
                r["R50_mm"].append(e["R50_mm"])
                r["R50_err_mm"].append(e["z0_err_mm"])
                r["Rp_mm"].append(e["Rp_mm"])
                r["Rp_err_mm"].append(e["Rp_err_mm"])
                r["chi2_ndf"].append(e["chi2_ndf"])

    # Summary: Δs against the dose endpoints, mean ± sem + shard std.
    tree = {"meta": {"shards": len(list(SHARDS)), "roi": "whole-plane",
                     "baseline": "free", "models": ",".join(MODELS),
                     "window_lo_mm": window[0], "window_hi_mm": window[1]},
            "dose": {mk: {k: dose[mk][k] for k in
                          ("R50_mm", "Rp_mm", "z0_err_mm", "Rp_err_mm")}
                     for mk in MODELS},
            "truth_activity": {
                mk: {"R50_mm": truth_act[mk]["R50_mm"],
                     "Rp_mm": truth_act[mk]["Rp_mm"],
                     "delta_R50_mm":
                         truth_act[mk]["R50_mm"] - dose[mk]["R50_mm"],
                     "delta_Rp_mm":
                         truth_act[mk]["Rp_mm"] - dose[mk]["Rp_mm"]}
                for mk in MODELS}}
    print(f"\n=== Δ = activity − dose, 10 shards (mean ± sem, [std]) ===")
    for rung in RUNGS:
        tree[rung] = {}
        for mk in MODELS:
            r = results[rung][mk]
            d50 = summarize(np.asarray(r["R50_mm"]) - dose[mk]["R50_mm"])
            dp = summarize(np.asarray(r["Rp_mm"]) - dose[mk]["Rp_mm"])
            tree[rung][mk] = dict(
                r,
                delta_R50_mean_mm=d50["mean"], delta_R50_sem_mm=d50["sem"],
                delta_R50_std_mm=d50["std"],
                delta_Rp_mean_mm=dp["mean"], delta_Rp_sem_mm=dp["sem"],
                delta_Rp_std_mm=dp["std"])
            print(f"  {rung:8s} {mk:8s} Δ_R50 {d50['mean']:+8.3f} ± "
                  f"{d50['sem']:.3f} [{d50['std']:.3f}]   Δ_Rp "
                  f"{dp['mean']:+8.3f} ± {dp['sem']:.3f} [{dp['std']:.3f}]")

    path = os.path.join(OUT, "results.toml")
    toml_dump(path, tree)
    print(f"wrote {path}")
    plot_deltas(results, dose, truth_act, "R50",
                os.path.join(OUT, "figures", "delta_r50.png"))
    plot_deltas(results, dose, truth_act, "Rp",
                os.path.join(OUT, "figures", "delta_rp.png"))


# --- the dose sweep -----------------------------------------------------------
# (dose_Gy, thin seeds); None = the ten full shards themselves.
DOSE_GROUPS = [(1.0, None), (0.5, 1), (0.2, 2), (0.1, 2)]


def dose_group_files(dose_gy, nseeds):
    if nseeds is None:
        return [os.path.join(CFG, "one_shard", f"recon_shard{i:03d}.npz")
                for i in SHARDS]
    return [os.path.join(OUT, "recons",
                         f"recon_shard{i:03d}_d{dose_gy:g}_s{s}.npz")
            for i in SHARDS for s in range(nseeds)]


def dose_sweep(params, centre, window, dose):
    groups = {}
    for dose_gy, nseeds in DOSE_GROUPS:
        r50 = {mk: [] for mk in MODELS}
        for path in dose_group_files(dose_gy, nseeds):
            d = np.load(path)
            z, prof = profile_from_image(d["image"], params["grid"], centre,
                                         None)
            res = analyze(os.path.basename(path)[:-4], z, prof, window, 0.0,
                          True, models=MODELS)
            for mk in MODELS:
                r50[mk].append(res[mk]["R50_mm"])
        groups[dose_gy] = r50

    sig1 = {mk: summarize(groups[1.0][mk])["std"] for mk in MODELS}
    tree = {"meta": {"roi": "whole-plane", "baseline": "free",
                     "models": ",".join(MODELS),
                     "window_lo_mm": window[0], "window_hi_mm": window[1]},
            "dose_fit": {mk: {"R50_mm": dose[mk]["R50_mm"]} for mk in MODELS}}
    print("\n=== Δ_R50 vs dose (mean ± sem, σ_R ± err, 1/√dose prediction) ===")
    for dose_gy, nseeds in DOSE_GROUPS:
        sec = {}
        for mk in MODELS:
            vals = np.asarray(groups[dose_gy][mk])
            s = summarize(vals - dose[mk]["R50_mm"])
            n = vals.size
            sig_err = s["std"] / np.sqrt(2 * (n - 1))
            pred = sig1[mk] * np.sqrt(1.0 / dose_gy)
            sec[mk] = {"n": n, "R50_mm": list(vals),
                       "delta_R50_mean_mm": s["mean"],
                       "delta_R50_sem_mm": s["sem"],
                       "sigma_R_mm": s["std"], "sigma_R_err_mm": sig_err,
                       "sigma_R_pred_mm": pred}
            if mk == "erfc":
                print(f"  {dose_gy:4.1f} Gy (n={n:2d})  Δ_R50 {s['mean']:+8.3f}"
                      f" ± {s['sem']:.3f}   σ_R {s['std']:.3f} ± {sig_err:.3f}"
                      f"   pred {pred:.3f}")
        tree[f"d{int(1000 * dose_gy)}mGy"] = sec

    path = os.path.join(OUT, "dose_sweep.toml")
    toml_dump(path, tree)
    print(f"wrote {path}")
    plot_dose_sweep(groups, dose, sig1,
                    os.path.join(OUT, "figures", "delta_r50_vs_dose.png"))


def plot_dose_sweep(groups, dose, sig1, path, mk="erfc"):
    """Top: Δ_R50 of every run vs dose (log x), mean ± sem per group, the
    1 Gy mean as the anchor line. Bottom: σ_R vs dose with the 1/√dose
    prediction anchored at 1 Gy as a thin red line."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.5, 6.8), facecolor=SURFACE,
                                 sharex=True, height_ratios=[2, 1])
    for a in (a1, a2):
        style(a)
    doses = [dg for dg, _ in DOSE_GROUPS]
    ref = dose[mk]["R50_mm"]
    anchor = None
    for dose_gy, _ in DOSE_GROUPS:
        d = np.asarray(groups[dose_gy][mk]) - ref
        s = summarize(d)
        if dose_gy == 1.0:
            anchor = s["mean"]
        a1.plot(np.full(d.size, dose_gy), d, "o", ms=5, mfc="none", mec=INK,
                mew=1.1, alpha=0.6,
                label="single runs" if dose_gy == doses[0] else None)
        a1.errorbar([dose_gy], [s["mean"]], yerr=[s["std"]], fmt="_", ms=20,
                    color=BLUE, elinewidth=1.8, capsize=5, zorder=5,
                    label="mean ± σ" if dose_gy == doses[0] else None)
        n = d.size
        a2.errorbar([dose_gy], [s["std"]],
                    yerr=[s["std"] / np.sqrt(2 * (n - 1))], fmt="o", ms=6,
                    color=INK, elinewidth=1.2, capsize=4)
    a1.axhline(anchor, color=BLUE, ls="--", lw=1.0, alpha=0.7,
               label="1 Gy mean")
    a1.set_xscale("log")
    a1.set_xticks(doses, [f"{d:g}" for d in doses])
    a1.set_ylabel("$\\Delta R$ [mm]", color=INK, fontsize=13)
    a1.legend(frameon=False, fontsize=12, labelcolor=INK, loc="best")
    dd = np.geomspace(min(doses), max(doses), 100)
    a2.plot(dd, sig1[mk] * np.sqrt(1.0 / dd), color=RED, lw=1.2,
            label="$\\sigma_R(1\\,\\mathrm{Gy})\\,\\sqrt{1\\,\\mathrm{Gy}/\\mathrm{dose}}$")
    a2.set_ylabel("$\\sigma_R$ [mm]", color=INK, fontsize=13)
    a2.set_xlabel("dose [Gy]", color=INK, fontsize=13)
    a2.legend(frameon=False, fontsize=12, labelcolor=INK, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


# --- the Fano test ------------------------------------------------------------
def fano(params, centre, window):
    """Variance/mean of the origins count per z bin across the ten shards.
    Poisson counts give Fano = 1 with per-bin sampling error √(2/(n−1));
    the window average tests the error assignment to a few percent."""
    counts = []
    for i in SHARDS:
        z, c = profile_from_origins(shard_h5(i), params["grid"], centre, None)
        counts.append(c)
    counts = np.stack(counts)
    n = counts.shape[0]
    mean = counts.mean(axis=0)
    var = counts.var(axis=0, ddof=1)
    ok = mean > 100          # bins with enough counts for a meaningful ratio
    f = np.where(ok, var / np.where(ok, mean, 1.0), np.nan)
    per_bin_se = np.sqrt(2.0 / (n - 1))

    def band(sel):
        v = f[sel & ok]
        return {"n_bins": int(v.size), "fano_mean": float(v.mean()),
                "fano_se": float(v.std(ddof=1) / np.sqrt(v.size)),
                "expected_se": float(per_bin_se / np.sqrt(v.size))}

    inw = (z >= window[0]) & (z <= window[1])
    res = {"meta": {"shards": n, "min_mean_counts": 100.0,
                    "per_bin_se": float(per_bin_se),
                    "window_lo_mm": window[0], "window_hi_mm": window[1]},
           "window": band(inw), "all_bins": band(np.ones_like(ok, bool))}
    w, a = res["window"], res["all_bins"]
    print(f"Fano (origins, {n} shards): window {w['fano_mean']:.3f} ± "
          f"{w['fano_se']:.3f} over {w['n_bins']} bins "
          f"(sampling se {w['expected_se']:.3f}); all bins "
          f"{a['fano_mean']:.3f} ± {a['fano_se']:.3f} ({a['n_bins']})")
    path = os.path.join(OUT, "fano.toml")
    toml_dump(path, res)
    print(f"wrote {path}")

    fig, ax = plt.subplots(figsize=(8.5, 4.6), facecolor=SURFACE)
    style(ax)
    ax.axvspan(window[0], window[1], color=GRIDC, alpha=0.5, lw=0)
    ax.errorbar(z[ok], f[ok], yerr=per_bin_se * f[ok], fmt="o", ms=4,
                color=INK, elinewidth=0.8, capsize=0, label="var/mean per bin")
    ax.axhline(1.0, color=RED, lw=1.0, label="Poisson (Fano = 1)")
    ax.set_xlabel("z [mm]", color=INK)
    ax.set_ylabel("Fano = var/mean", color=INK)
    ax.set_title("Origins counts across the ten shards: Fano factor per z bin",
                 color=INK, fontsize=11, loc="left")
    ax.legend(frameon=False, fontsize=9, labelcolor=INK, loc="best")
    fig.tight_layout()
    fpath = os.path.join(OUT, "figures", "fano_origins.png")
    fig.savefig(fpath, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {fpath}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dose-sweep", action="store_true",
                   help="fit the thinned dose realizations instead of the "
                        "rung ladder")
    p.add_argument("--fano", action="store_true",
                   help="Fano test of the origins bin errors instead of the "
                        "rung ladder")
    args = p.parse_args()
    params, centre, window, truth_act, dose = setup()
    if args.fano:
        fano(params, centre, window)
    elif args.dose_sweep:
        dose_sweep(params, centre, window, dose)
    else:
        ladder(params, centre, window, truth_act, dose)


if __name__ == "__main__":
    main()
