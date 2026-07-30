#!/usr/bin/env python3
"""Plot and export the reference BGO statistical-procedure result on the
data-driven source (nominal fitted cross sections). The figure and the paper's
washed macros are washed-thinning only; the nominal (no-washout) ensemble is
read solely to emit the StatNom* macros for the numerical finite-pool validation
(direct unwashed shards vs nominal thinned realisations)."""

import argparse
import math
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from crysp_paths import config_out
from fit_activity_profile import BLUE, INK, RED, SURFACE, style

SCENARIO = "uniform_headep_sobp_1e8_dd"
TOPOLOGY = "closed"
SCANNER = "crysp_ring_1m_bgo_2x0"
CRYSTAL = "bgo_195k_2X0"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def gaussian(x, mean, sigma):
    return np.exp(-0.5 * ((x - mean) / sigma) ** 2) / (
        sigma * math.sqrt(2 * math.pi)
    )


def plot_shard(path):
    with open(path, "rb") as stream:
        result = tomllib.load(stream)
    z = np.asarray(result["profile_z_mm"], dtype=float)
    profile = np.asarray(result["profile"], dtype=float)
    base, amplitude, r50, width = np.asarray(result["erfc_popt"], dtype=float)
    edge = base + 0.5 * amplitude * np.vectorize(math.erfc)(
        (z - r50) / (math.sqrt(2.0) * width)
    )
    fig, axis = plt.subplots(figsize=(6.4, 4.2), facecolor=SURFACE)
    style(axis)
    axis.plot(z, profile, "o", ms=3.5, color=INK, label="reconstructed profile")
    axis.plot(z, edge, color=BLUE, lw=2.0, label="bounded erfc fit")
    axis.axvspan(*result["fit_window_mm"], color=BLUE, alpha=0.08, label="fit window")
    axis.axvline(r50, color=RED, lw=1.5, ls="--", label=rf"$R_{{50}}={r50:.3f}$ mm")
    axis.set_xlabel("depth [mm]", color=INK)
    axis.set_ylabel("reconstructed activity [a.u.]", color=INK)
    axis.set_title(f"Shard {result['index']:03d}: $\\chi^2$/dof = {result['erfc_chi2_dof']:.2f}",
                   loc="left", color=INK)
    axis.legend(frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout()
    figure = os.path.splitext(path)[0] + ".png"
    fig.savefig(figure, dpi=220, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {figure}")


def main(scanner=SCANNER, crystal=CRYSTAL, out_name="statistical_procedure_bgo.png",
         write_macros=True, title_tag=""):
    root = os.path.join(
        config_out(SCENARIO, TOPOLOGY, scanner, crystal),
        "statistical_procedure", "del120s_ac300s_1Gy_D1p0Gy",
    )
    shard_files = sorted(
        os.path.join(root, "shards", name)
        for name in os.listdir(os.path.join(root, "shards"))
        if name.endswith(".toml")
    )
    shard_rows = [tomllib.load(open(path, "rb")) for path in shard_files]
    washed = tomllib.load(open(os.path.join(root, "combined", "washed_N100.toml"), "rb"))
    r_shards = np.asarray([row["R50_mm"] for row in shard_rows], dtype=float)
    r_washed = np.asarray(washed["R50_mm"], dtype=float)
    shard_mean, shard_sigma = r_shards.mean(), r_shards.std(ddof=1)
    shard_chi2 = np.asarray([row["erfc_chi2_dof"] for row in shard_rows], dtype=float)
    representative = shard_rows[np.argsort(shard_chi2)[len(shard_rows) // 2]]

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), facecolor=SURFACE)
    for axis in axes:
        style(axis)

    # (a) A representative individual edge fit.
    axis = axes[0]
    z = np.asarray(representative["profile_z_mm"], dtype=float)
    profile = np.asarray(representative["profile"], dtype=float)
    base, amplitude, r50, width = np.asarray(representative["erfc_popt"], dtype=float)
    edge = base + 0.5 * amplitude * np.vectorize(math.erfc)(
        (z - r50) / (math.sqrt(2.0) * width)
    )
    axis.plot(z, profile, "o", ms=3, color=INK, label="reconstructed profile")
    axis.plot(z, edge, color=BLUE, lw=2, label="bounded erfc fit")
    axis.axvspan(*representative["fit_window_mm"], color=BLUE, alpha=0.08)
    axis.axvline(r50, color=RED, lw=1.5, ls="--")
    axis.set_title(f"(a) One independent 1-Gy simulation{title_tag}", loc="left", color=INK)
    axis.set_xlabel("depth [mm]", color=INK)
    axis.set_ylabel("reconstructed activity [a.u.]", color=INK)
    axis.legend(frameon=False, fontsize=8, labelcolor=INK)

    # (b) The washed production ensemble.
    axis = axes[1]
    bins = np.linspace(r_washed.min(), r_washed.max(), 11)
    axis.hist(r_washed, bins=bins, density=True, histtype="stepfilled",
              color=RED, alpha=0.28, label="washed thinning")
    x = np.linspace(bins[0], bins[-1], 400)
    axis.plot(x, gaussian(x, washed["mean_R50_mm"], washed["raw_sigma_R_mm"]),
              color=RED, lw=1.8, label="Gaussian from raw spread")
    axis.set_title("(b) Washed fixed-pool thinning, $N=100$", loc="left", color=INK)
    axis.set_xlabel("$R_{50}$ [mm]", color=INK)
    axis.set_ylabel("density", color=INK)
    axis.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout()

    figure = os.path.join(REPO, "latex", "figs", out_name)
    fig.savefig(figure, dpi=220, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {figure}")
    if not write_macros:
        return

    # Nominal (no-washout) thinning ensemble: used ONLY for the numerical
    # finite-pool validation (direct unwashed shards vs nominal thinned). It does
    # not enter the washed figure or the washed macros. Same reference case and
    # [120,420] s window as the shards it is validated against.
    nominal = tomllib.load(open(os.path.join(root, "combined", "nominal_N100.toml"), "rb"))
    nom_chi = np.nanmedian([
        tomllib.load(open(os.path.join(root, "nominal", f"realization{i:04d}.toml"), "rb"))["erfc_chi2_dof"]
        for i in nominal["indices"]])

    macros = os.path.join(REPO, "latex", "statistical_procedure_results.tex")
    lines = [
        "% Generated by tools/plot_statistical_procedure.py; do not edit.",
        rf"\newcommand{{\StatShardMean}}{{{shard_mean:.3f}}}",
        rf"\newcommand{{\StatShardSigma}}{{{shard_sigma:.3f}}}",
        rf"\newcommand{{\StatShardFitErrMin}}{{{np.nanmin([r['R50_fit_error_mm'] for r in shard_rows]):.3f}}}",
        rf"\newcommand{{\StatShardFitErrMax}}{{{np.nanmax([r['R50_fit_error_mm'] for r in shard_rows]):.3f}}}",
        rf"\newcommand{{\StatShardChiMedian}}{{{np.nanmedian(shard_chi2):.2f}}}",
        rf"\newcommand{{\StatNomMean}}{{{nominal['mean_R50_mm']:.3f}}}",
        rf"\newcommand{{\StatNomRaw}}{{{nominal['raw_sigma_R_mm']:.3f}}}",
        rf"\newcommand{{\StatNomCorrection}}{{{nominal['finite_pool_correction']:.3f}}}",
        rf"\newcommand{{\StatNomCorrected}}{{{nominal['corrected_sigma_R_mm']:.3f}}}",
        rf"\newcommand{{\StatNomChiMedian}}{{{nom_chi:.2f}}}",
        rf"\newcommand{{\StatNomFailures}}{{{nominal['n_fail']}}}",
        rf"\newcommand{{\StatWashMean}}{{{washed['mean_R50_mm']:.3f}}}",
        rf"\newcommand{{\StatWashRaw}}{{{washed['raw_sigma_R_mm']:.3f}}}",
        rf"\newcommand{{\StatWashCorrection}}{{{washed['finite_pool_correction']:.3f}}}",
        rf"\newcommand{{\StatWashCorrected}}{{{washed['corrected_sigma_R_mm']:.3f}}}",
        rf"\newcommand{{\StatWashChiMedian}}{{{np.nanmedian([tomllib.load(open(os.path.join(root, 'washed', f'realization{i:04d}.toml'), 'rb'))['erfc_chi2_dof'] for i in washed['indices']]):.2f}}}",
        rf"\newcommand{{\StatWashFailures}}{{{washed['n_fail']}}}",
    ]
    with open(macros, "w", encoding="utf-8") as stream:
        stream.write("\n".join(lines) + "\n")
    print(f"wrote {macros}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", help="plot one completed shard TOML")
    parser.add_argument("--cafov", action="store_true",
                        help="compact-AFOV BGO (crysp_r40_35cm); figure only, no macros")
    args = parser.parse_args()
    if args.shard:
        plot_shard(args.shard)
    elif args.cafov:
        main(scanner="crysp_r40_35cm_bgo_2x0", crystal="bgo_195k_2X0",
             out_name="statistical_procedure_cafov.png", write_macros=False,
             title_tag=" (CAFOV)")
    else:
        main()
