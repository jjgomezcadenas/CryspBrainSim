#!/usr/bin/env python3
"""
washout_sensitivity.py — per-parameter Mizuno sensitivity of the activity edge.

The P1 deliverable of latex/docs/systematics.tex: how much does R50 move when the
Mizuno washout parameters move within their uncertainties? The central washout
shift itself calibrates away (it is a fixed offset per setup), so the systematic
is the *spread about the central washed value*, not the shift from physical:

    dR50(p) = R50(washed at p) - R50(washed at central p).

Two views are produced for one acquisition window:

  (a) one-at-a-time, each Mizuno parameter moved by +/-1 sigma. Sum M_k = 1 is a
      constraint of the model, so a shifted weight is renormalised; the reported
      derivative is therefore "M_k up, the other two down pro rata".
  (b) the joint band, N_MC draws over all six parameters (the number quoted as
      the washout parameter systematic).

The mechanism is reported alongside: washout moves the edge only by reweighting
the isotope mixture, so g(O15)/g(C11) and the O15 window fraction are tabulated
with each dR50. A common scale on every g_i cancels exactly in R50.

Truth level (origins): the per-isotope truth depth columns are reweighted and
refitted, no reconstruction. Bias transfers scanner-independently (md/results.md:
truth +0.218 mm vs detected +0.20..0.25 mm), so this is the right first pass.

Frame note: truth/activity_profile_fast.csv is in the native frame, while the v2
run-parameter files carry the tumour-centred window (legacy + source_z_offset_mm
= 25.584, a shard attribute this tool does not read). The window is therefore
rebuilt here around the truth profile's own fitted edge using the config margins
-- [edge - proximal_margin, edge + distal_margin], which reproduces the config
window geometry exactly. Every reported quantity is a difference of R50 values,
so the choice of frame cancels.

Reads  config/washout_brain.toml (model, errors, physical half-lives, timing)
       config/run_parameters.toml (the fit margins)
       truth/activity_profile_fast.csv (per-isotope depth columns)
Writes out/<scenario>/washout/washout_sensitivity_t<t1>_<t2>.toml

Run:  python3 tools/washout_sensitivity.py                 # window from config
      python3 tools/washout_sensitivity.py --t1 120 --t2 420   # cbs.tex protocol
"""
import argparse
import os
import tomllib

import numpy as np

from crysp_paths import REPO, scenario_out
from fit_activity_profile import ACTIVE, TRUTH, fit_edge, seeds, toml_dump
from washout import g_factor

LN2 = np.log(2.0)
SCENARIO = ACTIVE.scenario
OUT = os.path.join(scenario_out(SCENARIO), "washout")
MC_SEED = 20260723
trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# The six Mizuno parameters, in report order: (label, index, kind).
PARAMS = [("M_fast", 0, "M"), ("M_medium", 1, "M"), ("M_slow", 2, "M"),
          ("T_fast", 0, "T"), ("T_medium", 1, "T"), ("T_slow", 2, "T")]


def load(config_name):
    """Model + errors + physical half-lives + timing, and the fit margins."""
    with open(os.path.join(REPO, "config", "washout_brain.toml"), "rb") as f:
        w = tomllib.load(f)
    with open(os.path.join(REPO, "config", config_name), "rb") as f:
        rp = tomllib.load(f)
    return w, (rp["window"]["proximal_margin_mm"],
               rp["window"]["distal_margin_mm"])


def r50(z, y, window):
    """R50 of one activity curve, quietly (whole-plane erfc, free baseline)."""
    return float(fit_edge("erfc", z, y, window, True, 0.01, True,
                          seeds(z, y, window))["R50_mm"])


def truth_window(z, P0, margins):
    """The config fit window rebuilt in the truth profile's own frame:
    [edge - proximal, edge + distal]. Seeded from the half-maximum crossing,
    then re-centred once on the fitted physical edge."""
    prox, dist = margins
    above = np.where(P0 >= 0.5 * P0.max())[0]
    edge = float(z[above[-1]])
    for _ in range(2):
        window = (edge - prox, edge + dist)
        edge = r50(z, P0, window)
    return (edge - prox, edge + dist)


def perturb(M0, T0, idx, kind, delta):
    """One parameter moved by `delta`. Weights are renormalised to keep the
    model constraint sum M_k = 1; half-lives are kept positive."""
    M, T = list(M0), list(T0)
    if kind == "M":
        M[idx] = max(M[idx] + delta, 0.0)
        s = sum(M)
        M = [x / s for x in M]
    else:
        T[idx] = max(T[idx] + delta, 0.1)
    return M, T


def g_all(M, T, lam, isotopes, t_irr, t1, t2):
    """Per-isotope survival factors g_i for one parameter point."""
    mu = [LN2 / x for x in T]
    return {i: float(g_factor(lam[i], M, mu, t_irr, t1, t2)) for i in isotopes}


def washed_r50(g, z, P, isotopes, window):
    """R50 of the mixture reweighted by g_i."""
    return r50(z, sum(g[i] * P[i] for i in isotopes), window)


def o15_fraction(g, N, isotopes):
    """O15 share of the washed window signal — the mechanism behind dR50."""
    tot = sum(g[i] * N[i] for i in isotopes)
    return g["O15"] * N["O15"] / tot


def sensitivity(mod, z, P, N, P0, window, t1, t2, n_mc):
    """The full per-window result: central g_i, the one-at-a-time table, and
    the joint band. `mod` bundles the model arrays read from config."""
    M0, T0, Merr, Terr, lam, isotopes, t_irr = mod
    g0 = g_all(M0, T0, lam, isotopes, t_irr, t1, t2)
    r50_phys = r50(z, P0, window)
    r50_0 = washed_r50(g0, z, P, isotopes, window)

    oat = {}
    for label, idx, kind in PARAMS:
        err = (Merr if kind == "M" else Terr)[idx]
        shifts = []
        for sign in (+1.0, -1.0):
            M, T = perturb(M0, T0, idx, kind, sign * err)
            g = g_all(M, T, lam, isotopes, t_irr, t1, t2)
            shifts.append((washed_r50(g, z, P, isotopes, window) - r50_0,
                           g["O15"] / g["C11"], o15_fraction(g, N, isotopes)))
        up, dn = shifts
        oat[label] = {"sigma": err, "dR50_up_mm": up[0], "dR50_dn_mm": dn[0],
                      "dR50_sym_mm": 0.5 * abs(up[0] - dn[0]),
                      "g_ratio_up": up[1], "f_O15_up": up[2]}

    rng = np.random.default_rng(MC_SEED)
    M0a, Merra = np.array(M0), np.array(Merr)
    T0a, Terra = np.array(T0), np.array(Terr)
    dr = []
    for _ in range(n_mc):
        M = np.clip(rng.normal(M0a, Merra), 0.0, None)
        M = M / M.sum()
        T = np.clip(rng.normal(T0a, Terra), 0.1, None)
        g = g_all(list(M), list(T), lam, isotopes, t_irr, t1, t2)
        dr.append(washed_r50(g, z, P, isotopes, window) - r50_0)
    dr = np.array(dr)

    return {"window_t1_s": t1, "window_t2_s": t2,
            "g_factor": g0, "g_ratio_O15_C11": g0["O15"] / g0["C11"],
            "f_O15": o15_fraction(g0, N, isotopes),
            "R50_physical_mm": r50_phys, "R50_washed_mm": r50_0,
            "delta_R50_washout_mm": r50_0 - r50_phys,
            "one_at_a_time": oat,
            "quadrature_mm": float(np.sqrt(sum(v["dR50_sym_mm"] ** 2
                                               for v in oat.values()))),
            "joint_band_mm": float(np.std(dr)),
            "joint_mean_mm": float(np.mean(dr))}


def nsigma_bound(mod, z, P, N, P0, window, t1, t2, oat, n):
    """Worst-case coherent n-sigma excursion: every Mizuno parameter displaced
    by n sigma in the direction that moved R50 the same way at 1 sigma, all at
    once. This is the blind zero-case -- the calibration stays at the mean while
    the pseudo-data is generated at the displaced point -- so the returned dR50
    is a bias bound, not a band.

    Weights move additively (renormalised, and they stay inside [0,1] out to
    5 sigma). Half-lives move MULTIPLICATIVELY, T -> T*exp(±n*sigma/T): they are
    positive quantities with large relative errors, and an additive Gaussian
    excursion drives T_slow (10191 +/- 2200) negative beyond ~4.6 sigma, which
    clips to a degenerate instantaneous component and fakes a large shift."""
    M0, T0, Merr, Terr, lam, isotopes, t_irr = mod
    M, T = list(M0), list(T0)
    for label, idx, kind in PARAMS:
        sign = 1.0 if oat[label]["dR50_up_mm"] >= 0 else -1.0
        err = (Merr if kind == "M" else Terr)[idx]
        if kind == "M":
            M[idx] = min(max(M[idx] + sign * n * err, 0.0), 1.0)
        else:
            T[idx] = T[idx] * float(np.exp(sign * n * err / T[idx]))
    M = [x / sum(M) for x in M]
    g = g_all(M, T, lam, isotopes, t_irr, t1, t2)
    r50_0 = washed_r50(g_all(M0, T0, lam, isotopes, t_irr, t1, t2),
                       z, P, isotopes, window)
    return {"n": n, "M": M, "T_s": T,
            "g_factor": g, "g_ratio_O15_C11": g["O15"] / g["C11"],
            "f_O15": o15_fraction(g, N, isotopes),
            "dR50_mm": washed_r50(g, z, P, isotopes, window) - r50_0}


def report(res, isotopes):
    """Print one window's result."""
    print("\ncentral g_i: "
          + "  ".join(f"{i} {res['g_factor'][i]:.4f}" for i in isotopes))
    print(f"g(O15)/g(C11) = {res['g_ratio_O15_C11']:.4f} | "
          f"f(O15) = {res['f_O15']:.4f}")
    print(f"R50 physical {res['R50_physical_mm']:+.3f} mm -> washed "
          f"{res['R50_washed_mm']:+.3f} mm "
          f"(shift {res['delta_R50_washout_mm']:+.3f} mm, calibrates away)")
    print(f"\n{'parameter':>10} {'1sigma':>9} {'g(O15)/g(C11)':>14} "
          f"{'f(O15)':>8} {'dR50 [mm]':>10}")
    for label, _, _ in PARAMS:
        v = res["one_at_a_time"][label]
        print(f"{label:>10} {v['sigma']:>9g} {v['g_ratio_up']:>14.4f} "
              f"{v['f_O15_up']:>8.4f} {v['dR50_up_mm']:>+10.4f}")
    print(f"{'quadrature':>10} {'':>9} {'':>14} {'':>8} "
          f"{res['quadrature_mm']:>10.4f}")
    print(f"\njoint band: dR50 = {res['joint_mean_mm']:+.4f} "
          f"+/- {res['joint_band_mm']:.4f} mm   [does NOT calibrate away]")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--t1", type=float, default=None,
                    help="acquisition opens [s] (default: config timing)")
    ap.add_argument("--t2", type=float, default=None,
                    help="acquisition closes [s] (default: config timing)")
    ap.add_argument("--scan", default=None,
                    help="comma-separated t1 values [s] scanned at fixed "
                         "--duration, instead of a single window")
    ap.add_argument("--duration", type=float, default=300.0,
                    help="acquisition length [s] used by --scan (default 300)")
    ap.add_argument("--n-mc", type=int, default=1000,
                    help="draws for the joint band (default 1000)")
    ap.add_argument("--nsigma", default="1,3,5",
                    help="coherent worst-case excursions to evaluate on the "
                         "single-window path (default 1,3,5)")
    ap.add_argument("--config", default="run_parameters.toml",
                    help="run-parameter file supplying the fit margins")
    args = ap.parse_args()

    w, margins = load(args.config)
    t_irr = w["timing"]["t_irr_s"]
    M0, T0 = list(w["model"]["M"]), list(w["model"]["T_s"])
    Merr, Terr = w["model"]["M_err"], w["model"]["T_s_err"]
    isotopes = w["physical"]["isotopes"]
    lam = {i: LN2 / T for i, T in zip(isotopes, w["physical"]["T_half_s"])}
    mod = (M0, T0, Merr, Terr, lam, isotopes, t_irr)

    d = np.genfromtxt(os.path.join(TRUTH, "activity_profile_fast.csv"),
                      delimiter=",", names=True)
    z = np.asarray(d["z_mm"], float)
    P = {i: np.asarray(d[i], float) for i in isotopes}
    N = {i: float(trapz(P[i], z)) for i in isotopes}
    P0 = sum(P.values())
    window = truth_window(z, P0, margins)

    print(f"t_irr {t_irr:g} s | fit window [{window[0]:.3f}, {window[1]:.3f}] mm "
          f"(truth frame, margins {margins[0]:g}/{margins[1]:g}) | "
          f"scenario {SCENARIO}")

    meta = {"t_irr_s": t_irr, "fit_window_lo_mm": window[0],
            "fit_window_hi_mm": window[1],
            "proximal_margin_mm": margins[0], "distal_margin_mm": margins[1],
            "n_mc": args.n_mc, "mc_seed": MC_SEED, "config": args.config,
            "frame": "native truth frame; all outputs are R50 differences",
            "level": "truth (origins); no reconstruction"}
    os.makedirs(OUT, exist_ok=True)

    if args.scan:
        t1s = [float(x) for x in args.scan.split(",")]
        print(f"\nscan at fixed {args.duration:g} s duration "
              f"(t1 measured on the config clock, t_irr = {t_irr:g} s)\n")
        print(f"{'t1':>6} {'t2':>6} {'g(O15)/g(C11)':>14} {'f(O15)':>8} "
              f"{'shift':>9} {'M_med':>9} {'band [mm]':>10}")
        rows = {}
        for t1 in t1s:
            t2 = t1 + args.duration
            res = sensitivity(mod, z, P, N, P0, window, t1, t2, args.n_mc)
            rows[f"t{int(round(t1))}"] = res
            print(f"{t1:>6g} {t2:>6g} {res['g_ratio_O15_C11']:>14.4f} "
                  f"{res['f_O15']:>8.4f} "
                  f"{res['delta_R50_washout_mm']:>+9.4f} "
                  f"{res['one_at_a_time']['M_medium']['dR50_up_mm']:>+9.4f} "
                  f"{res['joint_band_mm']:>10.4f}")
        meta["duration_s"] = args.duration
        path = os.path.join(OUT,
                            f"washout_sensitivity_scan_d{int(args.duration)}.toml")
        toml_dump(path, {"meta": meta, "scan": rows})
    else:
        t1 = args.t1 if args.t1 is not None else w["timing"]["t1_s"]
        t2 = args.t2 if args.t2 is not None else w["timing"]["t2_s"]
        res = sensitivity(mod, z, P, N, P0, window, t1, t2, args.n_mc)
        report(res, isotopes)
        bounds = {}
        print(f"\n{'excursion':>10} {'g(O15)/g(C11)':>14} {'f(O15)':>8} "
              f"{'dR50 [mm]':>10}")
        for n in [float(x) for x in args.nsigma.split(",")]:
            b = nsigma_bound(mod, z, P, N, P0, window, t1, t2,
                             res["one_at_a_time"], n)
            bounds[f"n{int(n)}"] = b
            print(f"{n:>9g}s {b['g_ratio_O15_C11']:>14.4f} "
                  f"{b['f_O15']:>8.4f} {b['dR50_mm']:>+10.4f}")
        res["nsigma_bounds"] = bounds
        tag = f"t{int(round(t1))}_{int(round(t2))}"
        path = os.path.join(OUT, f"washout_sensitivity_{tag}.toml")
        toml_dump(path, {"meta": meta, "result": res})
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
