#!/usr/bin/env python3
"""
fit_activity_profile.py — the endpoint-fit laboratory: rebuild the depth
profile P(z) from a reconstructed activity image, fit its distal falloff with
two models, and compare the predicted "activity dies here" point against the
proton dose edge, using the same constructions on both curves.

Inputs (all under the current configuration, via crysp_paths):
  one_shard/recon_shardNNN.npz     the reconstructed image λ(r) from
                                   drivers/one_shard.jl
  config/run_parameters.toml       grid, disc ROI, fixed fit window
  <crystal>/crystal.toml           the detector response (σ_xyz for the blur)
  PtCryspProds/<scenario>/truth/   activity_profile_fast.csv, depth_dose.csv

The profile is the transverse sum at each depth over the WHOLE plane (the
settled convention — a disc clips the depth-widening beam halo and shifts the
endpoint proximally); pass --roi R to clip to a disc instead. The constant
voxel-area factor is absorbed by the fit's scale parameters. Against the disc
ROI the recomputed profile is checked against the one the driver stored — a
wiring self-test.

Edge models fitted inside the window (choose with --model):
  erfc     P(z) = b + a/2 · erfc((z − z0)/(√2 w))
           z0 is R50 by construction; any level follows analytically. The
           "practical endpoint" extrapolates the tangent at z0 to the baseline:
           R_p = z0 + √(π/2)·w  (the electron-dosimetry construction; the erfc
           itself never reaches zero). Gaussian distal tail.
  sigmoid  P(z) = b + a / (1 + exp((z − z0)/r))
           the logistic of Zapien-Campos et al. (Med Phys 2025, their Eq. 3;
           they fit it without the baseline — reproduce that with
           --no-baseline). Same midpoint observable (z0 = R50 = their PAR);
           tangent endpoint R_p = z0 + 2r; crossing R_x = z0 + r·ln((1−x)/x).
           Exponential distal tail — heavier than the erfc's, so low-fraction
           crossings sit more distal while R50 and R_p barely move. At equal
           midpoint slope the shapes map as r = √(2π)/4 · w ≈ 0.627 w
           (stored as w_equiv_mm = r/0.627 for direct width comparison).
Plus, always, the TOML-only cross-check:
  ramp     a source with a true endpoint — activity constant, falling linearly
           to zero at z_end over length L — convolved with a Gaussian of FIXED
           width σ_blur (the known resolution; fixing it breaks the L–σ
           degeneracy). z_end is the fitted "activity = 0" crossing.

σ_blur defaults: the reconstruction uses the σ the simulation applied to the
hits (crystal.toml, 1.7 mm; the real-scanner 3.5 mm FWHM ≈ 1.49 mm applies to
future productions — override with --sigma-blur). The truth curves carry no
detector blur, σ_blur = 0 (the ramp is then piecewise linear). The positron
range (O15 FWHM 0.5 mm, cusp with heavy tails) is negligible in quadrature.

Options:
  --shard N          shard index (default 0)
  --all-uncorr       read the all-events-uncorrected reconstruction
  --source S         what carries the activity curve (default recon):
                       recon    the reconstructed image → disc-ROI profile
                       origins  the detected trues' TRUE annihilation
                                positions from the shard, histogrammed on the
                                same z grid — no reconstruction involved, so
                                recon-vs-origins comparisons isolate detector
                                + reconstruction effects from geometry ones
  --roi R            transverse ROI radius in mm to clip the profile to a
                     disc; the default is the WHOLE PLANE (the settled
                     convention). `none` is also accepted as an explicit
                     whole-plane request
  --window LO,HI     fit window in mm (default: the frozen run parameter;
                     the distal tail past the frozen edge constrains z_end,
                     so widening distally is a useful experiment)
  --model M          edge model(s): erfc | sigmoid | both (default erfc)
  --no-baseline      fix b = 0 in the edge models (the paper's 3-parameter
                     protocol; the ramp cross-check keeps its own baseline)
  --sigma-blur MM    Gaussian σ for the recon ramp model (default crystal.toml)
  --float-sigma      let σ float in the recon ramp fit (shows the L–σ
                     degeneracy; 5 parameters)
  --no-pulls         single-panel figures (data, fit, markers) without the
                     pull panel — the publication variant
  --show             display each figure and wait for return (PNGs are
                     written in every mode)

Writes one_shard/fits/fit_<tag>.toml + fits/figures/*.png, where the tag
carries the shard and any non-default source/ROI choice.
"""
import argparse
import os
import sys
import tomllib

import matplotlib
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc, erfcinv, expit

from crysp_paths import REPO, config_out

SCENARIO, TOPOLOGY, RING, CRYSTAL = (
    "uniform_headep_sobp_1e8", "closed", "crysp_ring_1m", "bgo_3X0")
CFG = config_out(SCENARIO, TOPOLOGY, RING, CRYSTAL)
TRUTH = os.path.join(os.path.dirname(REPO), "PtCryspProds", SCENARIO, "truth")

BLUE, AQUA, RED = "#2a78d6", "#1baf7a", "#e34948"
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
SURFACE = "#fcfcfb"

SQ2 = np.sqrt(2.0)
RP_C = np.sqrt(np.pi / 2.0)          # tangent-extrapolation coefficient
SLOPE_MATCH = np.sqrt(2.0 * np.pi) / 4.0   # r = SLOPE_MATCH·w at equal midpoint slope


# --- models ------------------------------------------------------------------
def norm_cdf(x):
    return 0.5 * erfc(-x / SQ2)


def norm_pdf(x):
    return np.exp(-0.5 * x * x) / np.sqrt(2.0 * np.pi)


def erfc_model(z, b, a, z0, w):
    return b + a * 0.5 * erfc((z - z0) / (SQ2 * w))


def sigmoid_model(z, b, a, z0, r):
    """The logistic edge (Zapien-Campos Eq. 3 plus a baseline): plateau `a`
    falling through half-height at `z0` with fall-off rate `r`."""
    return b + a * expit(-(z - z0) / r)


EDGE_FN = {"erfc": erfc_model, "sigmoid": sigmoid_model}


def edge_coeffs(model_key, zero_fraction):
    """Distances of R20, the tangent endpoint Rp, and the crossing Rx from
    the midpoint z0, in units of the shape parameter (w or r)."""
    if model_key == "erfc":
        return (SQ2 * erfcinv(0.4), RP_C, SQ2 * erfcinv(2 * zero_fraction))
    return (np.log(4.0), 2.0, np.log((1 - zero_fraction) / zero_fraction))


def ramp_blur(z, b, a, z_end, L, sigma):
    """Plateau `a` falling linearly to zero at `z_end` over `L`, convolved
    with a Gaussian of width `sigma`, plus baseline `b`. At sigma = 0 this is
    the bare piecewise-linear source."""
    L = max(L, 1e-6)
    if sigma <= 0:
        return b + a * np.clip((z_end - z) / L, 0.0, 1.0)
    t1 = z_end - L
    al = (t1 - z) / sigma
    be = (z_end - z) / sigma
    ramp = ((z_end - z) * (norm_cdf(be) - norm_cdf(al))
            + sigma * (norm_pdf(be) - norm_pdf(al))) / L
    return b + a * (norm_cdf(al) + ramp)


def _selftest_ramp_blur():
    """The analytic convolution against a brute-force numeric one."""
    z = np.linspace(-40, 40, 4001)
    src = ramp_blur(z, 0.0, 1.0, 10.0, 12.0, 0.0)
    g = norm_pdf((z - z.mean()) / 2.5) / 2.5
    num = np.convolve(src, g / g.sum(), mode="same")
    ana = ramp_blur(z, 0.0, 1.0, 10.0, 12.0, 2.5)
    core = (z > -25) & (z < 30)
    assert np.allclose(num[core], ana[core], atol=2e-3), "ramp_blur analytic form"


# --- fitting -----------------------------------------------------------------
def fit_curve(z, y, window, model, p0, bounds, weighted):
    """Fit `model` inside `window`; return (popt, pcov, chi2, ndf)."""
    sel = (z >= window[0]) & (z <= window[1])
    zf, yf = z[sel], y[sel]
    sig = np.sqrt(np.clip(yf, 1.0, None)) if weighted else None
    popt, pcov = curve_fit(model, zf, yf, p0=p0, sigma=sig,
                           absolute_sigma=weighted, bounds=bounds, maxfev=40000)
    resid = yf - model(zf, *popt)
    chi2 = float(np.sum((resid / (sig if weighted else 1.0)) ** 2))
    return popt, pcov, chi2, zf.size - len(popt)


def seeds(z, y, window):
    """Common initial values: baseline from the distal fifth of the window,
    plateau from the proximal fifth, position from the half-height crossing."""
    sel = (z >= window[0]) & (z <= window[1])
    zf, yf = z[sel], y[sel]
    k = max(2, yf.size // 5)
    b0, p0 = float(np.median(yf[-k:])), float(np.median(yf[:k]))
    a0 = max(p0 - b0, 1e-9)
    above = np.where(yf >= b0 + 0.5 * a0)[0]
    zc = float(zf[above[-1]]) if above.size else float(zf[zf.size // 2])
    return b0, a0, zc, float(0.15 * (zf[-1] - zf[0]))


def fit_edge(model_key, z, y, window, weighted, zero_fraction, baseline, seed):
    """Fit one edge model (with or without the baseline term) and return its
    parameters plus the derived endpoints with propagated errors."""
    b0, a0, zc, w0 = seed
    s0 = w0 if model_key == "erfc" else SLOPE_MATCH * w0
    fn = EDGE_FN[model_key]
    if baseline:
        model, p0 = fn, [b0, a0, zc, s0]
        lo = [-np.inf, 0.0, window[0] - 30, 0.05]
        hi = [np.inf, np.inf, window[1] + 30, 60.0]
        iz, iw = 2, 3
    else:
        model = lambda zz, a_, z0_, s_: fn(zz, 0.0, a_, z0_, s_)
        p0 = [max(a0 + b0, 1e-9), zc, s0]
        lo = [0.0, window[0] - 30, 0.05]
        hi = [np.inf, window[1] + 30, 60.0]
        iz, iw = 1, 2
    popt, pcov, chi2, ndf = fit_curve(z, y, window, model, p0, (lo, hi), weighted)
    b = float(popt[0]) if baseline else 0.0
    a, z0, s = float(popt[iz - 1]), float(popt[iz]), float(popt[iw])
    c20, cp, cx = edge_coeffs(model_key, zero_fraction)

    def derived(c):
        err = float(np.sqrt(pcov[iz, iz] + c**2 * pcov[iw, iw]
                            + 2 * c * pcov[iz, iw]))
        return z0 + c * s, err

    out = {"b": b, "a": a, "z0_mm": z0,
           ("w_mm" if model_key == "erfc" else "r_mm"): s,
           "z0_err_mm": float(np.sqrt(pcov[iz, iz])),
           "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf,
           "baseline": baseline, "zero_fraction": zero_fraction,
           "R50_mm": z0}
    out["R20_mm"], out["R20_err_mm"] = derived(c20)
    out["Rp_mm"], out["Rp_err_mm"] = derived(cp)
    # The crossing: the fitted signal falls to `zero_fraction` of the plateau.
    # Any fraction below 10.5% (the level the tangent construction touches)
    # sits distal of Rp.
    out["Rx_mm"], out["Rx_err_mm"] = derived(cx)
    if model_key == "sigmoid":
        out["w_equiv_mm"] = s / SLOPE_MATCH   # the equal-slope dictionary
    # The strict P = 0 crossing exists only for a negative fitted baseline.
    if -a < b < 0:
        out["z_zero_mm"] = (z0 + SQ2 * s * erfcinv(-2 * b / a)
                            if model_key == "erfc"
                            else z0 + s * np.log(-a / b - 1.0))
    else:
        out["z_zero_mm"] = float("nan")
    return out


def analyze(name, z, y, window, sigma_blur, weighted, models=("erfc",),
            baseline=True, zero_fraction=0.01, float_sigma=False):
    """Fit the requested edge models plus the ramp cross-check to one curve;
    return a flat result dict."""
    seed = seeds(z, y, window)
    b0, a0, zc, w0 = seed
    res = {"window_lo_mm": window[0], "window_hi_mm": window[1],
           "sigma_blur_mm": sigma_blur}
    for mk in models:
        res[mk] = fit_edge(mk, z, y, window, weighted, zero_fraction,
                           baseline, seed)

    # ramp ⊗ Gaussian: z_end is the fitted zero-crossing.
    if float_sigma:
        model = lambda zz, b_, a_, ze, L, s: ramp_blur(zz, b_, a_, ze, L, s)
        p0 = [b0, a0, zc + w0, 2.5 * w0, max(sigma_blur, 1.0)]
        lo = [-np.inf, 0, window[0] - 30, 0.5, 0.05]
        hi = [np.inf, np.inf, window[1] + 30, 80, 15]
    else:
        model = lambda zz, b_, a_, ze, L: ramp_blur(zz, b_, a_, ze, L, sigma_blur)
        p0 = [b0, a0, zc + w0, 2.5 * w0]
        lo = [-np.inf, 0, window[0] - 30, 0.5]
        hi = [np.inf, np.inf, window[1] + 30, 80]
    popt, pcov, chi2, ndf = fit_curve(z, y, window, model, p0, (lo, hi), weighted)
    ramp = {"b": popt[0], "a": popt[1], "z_end_mm": popt[2], "L_mm": popt[3],
            "z_end_err_mm": float(np.sqrt(pcov[2, 2])),
            "chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf}
    if float_sigma:
        ramp["sigma_mm"] = popt[4]
        ramp["sigma_err_mm"] = float(np.sqrt(pcov[4, 4]))
    res["ramp"] = ramp

    r = res["ramp"]
    print(f"\n{name} (window {window[0]:.2f}..{window[1]:.2f}, "
          f"σ_blur {sigma_blur:g} mm, baseline {'free' if baseline else 'b=0'}):")
    for mk in models:
        e = res[mk]
        sname = "w" if mk == "erfc" else "r"
        print(f"  {mk:7s} χ²/ndf {e['chi2_ndf']:10.3g} | R50 {e['R50_mm']:8.3f}"
              f" ± {e['z0_err_mm']:.3f}"
              f" | Rp {e['Rp_mm']:8.3f} ± {e['Rp_err_mm']:.3f}"
              f" | R{100 * e['zero_fraction']:g}% {e['Rx_mm']:8.3f} ± "
              f"{e['Rx_err_mm']:.3f} | {sname} {e[sname + '_mm']:.2f}")
    extra = (f" | σ {ramp['sigma_mm']:.2f} ± {ramp['sigma_err_mm']:.2f}"
             if float_sigma else "")
    print(f"  ramp    χ²/ndf {r['chi2_ndf']:10.3g} | z_end {r['z_end_mm']:8.3f} ± "
          f"{r['z_end_err_mm']:.3f} | L {r['L_mm']:.2f}{extra}")
    return res


def edge_window(z, y, margins):
    """The paper's window construction applied to THIS curve: the frozen
    proximal/distal margins around the curve's own half-height edge (last
    downward crossing of half-max). The dose edge sits ~13 mm distal of the
    activity edge, so each curve gets a window bracketing its own falloff."""
    thr = 0.5 * y.max()
    idx = np.where((y[:-1] >= thr) & (y[1:] < thr))[0]
    i = idx[-1]
    zc = z[i] + (thr - y[i]) * (z[i + 1] - z[i]) / (y[i + 1] - y[i])
    return (zc - margins[0], zc + margins[1])


# --- profile from the image --------------------------------------------------
def profile_from_image(image, grid, centre, roi_radius):
    """P(z): transverse sum at each depth slice — over the disc of
    `roi_radius` (mm) around `centre`, or the whole plane when None."""
    n = grid["n"]
    org, vs = grid["img_origin_mm"], grid["voxsize_mm"]
    assert tuple(image.shape) == tuple(n), f"image shape {image.shape} ≠ grid {n}"
    z = org[2] + vs[2] * np.arange(n[2])
    if roi_radius is None:
        return z, image.sum(axis=(0, 1))
    x = org[0] + vs[0] * np.arange(n[0])
    y = org[1] + vs[1] * np.arange(n[1])
    X, Y = np.meshgrid(x, y, indexing="ij")
    mask = (X - centre[0]) ** 2 + (Y - centre[1]) ** 2 <= roi_radius ** 2
    return z, image[mask, :].sum(axis=0)


def profile_from_origins(shard_file, grid, centre, roi_radius):
    """P(z) from the detected trues' TRUE annihilation positions: histogram
    z0 on the grid's z bins, within the transverse disc (or the whole plane).
    Genuinely Poisson counts, and no reconstruction in the chain."""
    import h5py
    with h5py.File(shard_file, "r") as f:
        s = float(f.attrs["xyz_scale_mm"])
        keep = f["truth"][:] == 0
        x0 = f["x0_mm"][:].astype(np.float64)[keep] * s
        y0 = f["y0_mm"][:].astype(np.float64)[keep] * s
        z0 = f["z0_mm"][:].astype(np.float64)[keep] * s
    if roi_radius is not None:
        sel = (x0 - centre[0]) ** 2 + (y0 - centre[1]) ** 2 <= roi_radius ** 2
        z0 = z0[sel]
    org, vs, n = grid["img_origin_mm"], grid["voxsize_mm"], grid["n"]
    z = org[2] + vs[2] * np.arange(n[2])
    edges = np.concatenate([z - vs[2] / 2, [z[-1] + vs[2] / 2]])
    counts, _ = np.histogram(z0, bins=edges)
    return z, counts.astype(float)


# --- plotting ----------------------------------------------------------------
def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def plot_curve(plt, name, z, y, res, model_key, weighted, figdir, show,
               pulls=True):
    """Data with error bars, one edge-model fit as a thin red line, the fixed
    window shaded, R_p marked; lower panel shows the fit pulls (weighted
    curves) or the raw residuals (dose). `pulls=False` drops the lower panel
    (the publication variant)."""
    win = (res["window_lo_mm"], res["window_hi_mm"])
    e = res[model_key]
    fn = EDGE_FN[model_key]
    s = e["w_mm"] if model_key == "erfc" else e["r_mm"]
    zf = np.linspace(win[0] - 12, win[1] + 10, 800)
    yerr = np.sqrt(np.clip(y, 1.0, None)) if weighted else None

    if pulls:
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 6.4),
                                     facecolor=SURFACE, sharex=True,
                                     height_ratios=[3, 1])
        axes = (a1, a2)
    else:
        fig, a1 = plt.subplots(figsize=(9.5, 5.0), facecolor=SURFACE)
        a2, axes = None, (a1,)
    for a in axes:
        style(a)
    a1.axvspan(*win, color=GRIDC, alpha=0.5, lw=0)
    a1.errorbar(z, y, yerr=yerr, fmt="o", ms=3, color=INK, mec="none",
                elinewidth=0.8, capsize=0, label="data")
    a1.plot(zf, fn(zf, e["b"], e["a"], e["z0_mm"], s),
            color=RED, lw=1.0,
            label=f"{model_key} fit  χ²/ndf {e['chi2_ndf']:.3g}")
    a1.axvline(e["Rp_mm"], color=BLUE, ls=":", lw=1.2,
               label=f"Rp {e['Rp_mm']:.2f} ± {e['Rp_err_mm']:.2f} mm")
    a1.axvline(e["Rx_mm"], color=AQUA, ls=":", lw=1.2,
               label=f"R{100 * e['zero_fraction']:g}% {e['Rx_mm']:.2f} ± "
                     f"{e['Rx_err_mm']:.2f} mm")
    a1.set_xlim(win[0] - 12, max(win[1] + 10, e["Rx_mm"] + 4))
    sel = (z > win[0] - 12) & (z < win[1] + 10)
    a1.set_ylim(min(0, 1.05 * y[sel].min()), 1.15 * y[sel].max())
    a1.set_ylabel("P(z)", color=INK)
    a1.set_title(f"{name}: distal-edge fit and endpoint ({model_key})",
                 color=INK, fontsize=11, loc="left")
    a1.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc="upper right")

    if a2 is not None:
        inw = (z >= win[0]) & (z <= win[1])
        fit = fn(z[inw], e["b"], e["a"], e["z0_mm"], s)
        a2.axhline(0, color=MUTED, lw=1)
        if weighted:
            pull = (y[inw] - fit) / np.sqrt(np.clip(y[inw], 1.0, None))
            a2.errorbar(z[inw], pull, yerr=1.0, fmt="o", ms=2.5, color=INK,
                        elinewidth=0.8, capsize=0)
            a2.set_ylabel("pull", color=INK)
        else:
            a2.plot(z[inw], y[inw] - fit, "o", ms=2.5, color=INK)
            a2.set_ylabel("residual", color=INK)
        a2.set_xlabel("z [mm]", color=INK)
    else:
        a1.set_xlabel("z [mm]", color=INK)

    fig.tight_layout()
    path = os.path.join(figdir, f"{name}.png")
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")
    if show:
        plt.show(block=False)
        plt.pause(0.15)
        input("[return to continue] ")
    plt.close(fig)


# --- TOML output -------------------------------------------------------------
def toml_dump(path, tree):
    """Minimal TOML writer for nested dicts of scalars and flat arrays."""
    def fmt(v):
        if isinstance(v, str):
            return f'"{v}"'
        if isinstance(v, (bool, np.bool_)):
            return str(bool(v)).lower()
        if isinstance(v, (int, np.integer)):
            return str(int(v))
        if isinstance(v, (list, tuple, np.ndarray)):
            return "[" + ", ".join(fmt(x) for x in v) + "]"
        return repr(float(v))

    with open(path, "w") as io:
        def emit(d, prefix=""):
            subs = []
            for k, v in d.items():
                if isinstance(v, dict):
                    subs.append((k, v))
                else:
                    io.write(f"{k} = {fmt(v)}\n")
            for k, v in subs:
                name = f"{prefix}{k}"
                io.write(f"\n[{name}]\n")
                emit(v, prefix=f"{name}.")
        emit(tree)


# --- main --------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--all-uncorr", action="store_true")
    p.add_argument("--source", choices=["recon", "origins"], default="recon")
    p.add_argument("--roi", default=None,
                   help="disc ROI radius in mm; default is the whole plane")
    p.add_argument("--window", default=None, help="LO,HI in mm")
    p.add_argument("--model", choices=["erfc", "sigmoid", "both"],
                   default="erfc", help="edge model(s) to fit (default erfc)")
    p.add_argument("--no-baseline", action="store_true",
                   help="fix b = 0 in the edge models (the paper's protocol)")
    p.add_argument("--sigma-blur", type=float, default=None,
                   help="recon blur σ in mm (default: crystal.toml)")
    p.add_argument("--zero-fraction", type=float, default=0.01,
                   help="plateau fraction defining the crossing Rx (default 1%%)")
    p.add_argument("--float-sigma", action="store_true")
    p.add_argument("--no-pulls", action="store_true",
                   help="single-panel figures without the pull panel")
    p.add_argument("--show", action="store_true")
    args = p.parse_args()

    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    _selftest_ramp_blur()

    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        params = tomllib.load(f)
    with open(os.path.join(CFG, "crystal.toml"), "rb") as f:
        crystal = tomllib.load(f)
    sigma_recon = (args.sigma_blur if args.sigma_blur is not None
                   else crystal["detector"]["sigma_xyz_mm"])
    window = (tuple(float(v) for v in args.window.split(",")) if args.window
              else (params["window"]["z_lo_mm"], params["window"]["z_hi_mm"]))

    frozen_roi = params["roi"]["radius_mm"]
    centre = params["roi"]["centre_mm"]
    # Whole plane by default (the settled convention); a numeric --roi clips
    # to a disc, `none` is an explicit whole-plane request.
    if args.roi is None or args.roi.lower() == "none":
        roi_radius = None
    else:
        roi_radius = float(args.roi)

    shard_tag = f"shard{args.shard:03d}" + ("_all_uncorr" if args.all_uncorr else "")
    if args.source == "recon":
        d = np.load(os.path.join(CFG, "one_shard", f"recon_{shard_tag}.npz"))
        z, prof = profile_from_image(d["image"], params["grid"], centre, roi_radius)
        if roi_radius == frozen_roi and not args.all_uncorr:
            assert np.allclose(prof, d["profile"], rtol=1e-6), \
                "recomputed profile ≠ driver profile (ROI/grid wiring)"
        act_name = ("recon_all_events_activity" if args.all_uncorr
                    else "recon_activity")
    else:
        shard_file = os.path.join(
            os.path.dirname(REPO), "PtCryspProds", SCENARIO, RING, "bgo",
            "fast_1Gy", f"lors_shard{args.shard:03d}.h5")
        z, prof = profile_from_origins(shard_file, params["grid"], centre,
                                       roi_radius)
        act_name = "origins_activity"

    models = (("erfc", "sigmoid") if args.model == "both" else (args.model,))
    baseline = not args.no_baseline
    roi_tag = "" if roi_radius is None else f"_roi{roi_radius:g}"
    nob_tag = "_nob" if args.no_baseline else ""
    tag = (shard_tag + ("" if args.source == "recon" else "_origins")
           + roi_tag + nob_tag)

    ta = np.genfromtxt(os.path.join(TRUTH, "activity_profile_fast.csv"),
                       delimiter=",", names=True)
    dd = np.genfromtxt(os.path.join(TRUTH, "depth_dose.csv"),
                       delimiter=",", names=True)

    out = os.path.join(CFG, "one_shard", "fits")
    figdir = os.path.join(out, "figures")
    os.makedirs(figdir, exist_ok=True)

    # Activity curves use the frozen window; the dose falloff sits ~13 mm
    # distal of the activity edge, so it gets the same margins around its own
    # edge (else the window truncates it and z_end is unconstrained). An
    # explicit --window overrides everywhere.
    margins = (params["window"]["proximal_margin_mm"],
               params["window"]["distal_margin_mm"])
    zd = np.asarray(dd["z_mm"], float)
    yd = np.asarray(dd["dose_core_Gy"], float)
    dose_window = window if args.window else edge_window(zd, yd, margins)

    results = {"meta": {"tag": tag, "source": args.source,
                        "roi_mm": (-1.0 if roi_radius is None else roi_radius),
                        "sigma_blur_recon_mm": sigma_recon,
                        "window_lo_mm": window[0], "window_hi_mm": window[1],
                        "models": ",".join(models), "baseline": baseline,
                        "float_sigma": args.float_sigma}}
    curves = [
        (act_name, z, prof, window, sigma_recon if args.source == "recon" else 0.0,
         True),
        ("truth_activity", ta["z_mm"], ta["total"], window, 0.0, True),
        ("truth_dose", zd, yd, dose_window, 0.0, False),
    ]
    for name, zc, yc, win, sb, weighted in curves:
        res = analyze(name, np.asarray(zc, float), np.asarray(yc, float),
                      win, sb, weighted, models=models, baseline=baseline,
                      zero_fraction=args.zero_fraction,
                      float_sigma=(args.float_sigma and name == act_name))
        results[name] = res
        base_name = (name + roi_tag + nob_tag) if name == act_name \
            else (name + nob_tag)
        for mk in models:
            fig_name = base_name + ("" if mk == "erfc" else "_sigmoid")
            plot_curve(plt, fig_name, np.asarray(zc, float),
                       np.asarray(yc, float), res, mk, weighted, figdir,
                       args.show, pulls=not args.no_pulls)

    # The headline: activity endpoints against the dose endpoints, same
    # construction (same model, same window recipe) on both curves.
    ra, rd = results[act_name], results["truth_dose"]
    tact = results["truth_activity"]
    print(f"\nendpoint comparison (activity vs dose, same construction; "
          f"activity = {act_name}{roi_tag}{nob_tag}):")
    comparison = {}
    for mk in models:
        pre = "" if mk == "erfc" else f"{mk}_"
        for obs in ("R50", "Rp"):
            act, ta_, ds = (ra[mk][f"{obs}_mm"], tact[mk][f"{obs}_mm"],
                            rd[mk][f"{obs}_mm"])
            print(f"  {mk:7s} {obs:3s} : act {act:8.3f} | truth act {ta_:8.3f}"
                  f" | dose {ds:8.3f}  → act−dose = {act - ds:+.3f} mm")
            comparison[f"{pre}{obs}_act_minus_dose_mm"] = act - ds
            comparison[f"{pre}{obs}_truthact_minus_dose_mm"] = ta_ - ds
    print(f"  ramp z_end  : act {ra['ramp']['z_end_mm']:8.3f} | truth act "
          f"{tact['ramp']['z_end_mm']:8.3f} | dose {rd['ramp']['z_end_mm']:8.3f}"
          f"  → act−dose = {ra['ramp']['z_end_mm'] - rd['ramp']['z_end_mm']:+.3f} mm")
    comparison["z_end_act_minus_dose_mm"] = (
        ra["ramp"]["z_end_mm"] - rd["ramp"]["z_end_mm"])
    comparison["z_end_truthact_minus_dose_mm"] = (
        tact["ramp"]["z_end_mm"] - rd["ramp"]["z_end_mm"])
    results["comparison"] = comparison

    path = os.path.join(out, f"fit_{tag}.toml")
    toml_dump(path, results)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
