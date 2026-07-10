#!/usr/bin/env python3
"""
washout.py — the truth-level isotope-washout (IW) study: reweight the per-isotope
truth activity columns by each isotope's window-integrated survival factor g_i
(latex/washout_brain.tex Eq. 7), refit the distal edge, and report the shift.

Washout-as-loss on a spatially-uniform brain (config/washout_brain.toml) reduces
to a per-isotope scalar

    g_i = Σ_k M_k Φ(λ_i^phys + λ_k^bio) / Φ(λ_i^phys),
    Φ(a) = (exp(a·t_irr) − 1)(exp(−a·t1) − exp(−a·t2)) / a² ,

computed here from central Mizuno values (no error propagation this pass). The
washed profile is P^wo(z) = Σ_i g_i P_i(z); washout distorts no single isotope's
shape, so the edge moves only through the mix. The gating number is
ΔR50^wo = R50(P^wo) − R50(P^0): if it is small the detected study is moot.

g_i is cross-checked against a direct 2-D numerical integration of S_i/S_i^0.

Reads  config/washout_brain.toml, config/run_parameters.toml (the fit window),
       truth/activity_profile_fast.csv (per-isotope depth columns).
Writes out/<scenario>/washout/washout.toml (g_i, kept fraction, ΔR50^wo)
       out/<scenario>/washout/figures/washout_profile.png, g_factors.png

Run:  python3 tools/washout.py
"""
import os
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import dblquad

from crysp_paths import REPO, scenario_out
from fit_activity_profile import (
    ACTIVE, BLUE, INK, MUTED, RED, SURFACE, TRUTH, analyze, fit_edge, seeds,
    style, toml_dump)

LN2 = np.log(2.0)
SCENARIO = ACTIVE.scenario
OUT = os.path.join(scenario_out(SCENARIO), "washout")
N_MC = 1000                # Monte-Carlo draws for the washout-parameter band
MC_SEED = 20260710
trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def load_params():
    with open(os.path.join(REPO, "config", "washout_brain.toml"), "rb") as f:
        w = tomllib.load(f)
    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        rp = tomllib.load(f)
    window = (rp["window"]["z_lo_mm"], rp["window"]["z_hi_mm"])
    return w, window


def phi(a, t_irr, t1, t2):
    """The elementary window/build-up integral of note Eq. 6."""
    a = np.asarray(a, float)
    return (np.expm1(a * t_irr) * (np.exp(-a * t1) - np.exp(-a * t2))) / a**2


def g_factor(lam_phys, M, mu_bio, t_irr, t1, t2):
    """Per-isotope survival factor g_i (note Eq. 7), central values."""
    num = sum(Mk * phi(lam_phys + muk, t_irr, t1, t2)
              for Mk, muk in zip(M, mu_bio))
    return num / phi(lam_phys, t_irr, t1, t2)


def g_factor_numeric(lam_phys, M, mu_bio, t_irr, t1, t2):
    """g_i by direct 2-D integration of S_i/S_i^0 — the closed-form check.
    Uniform production R_i=1 on [0, t_irr]; inner t' integral runs to
    min(t, t_irr) = t_irr since t ≥ t1 > t_irr."""
    def S(kernel):
        val, _ = dblquad(lambda tp, t: kernel(t - tp),
                         t1, t2, lambda t: 0.0, lambda t: t_irr)
        return val
    s0 = S(lambda dt: lam_phys * np.exp(-lam_phys * dt))
    sw = S(lambda dt: lam_phys * sum(Mk * np.exp(-(lam_phys + muk) * dt)
                                     for Mk, muk in zip(M, mu_bio)))
    return sw / s0


def refit_r50(z, y, window):
    """R50 of an activity curve, quietly (no per-fit printout for the MC loop)."""
    return fit_edge("erfc", z, y, window, True, 0.01, True,
                    seeds(z, y, window))["R50_mm"]


def mc_band(z, P, isotopes, lam, window, w, R50_phys):
    """The washout-parameter systematic: draw the Mizuno (M_k, T_k) from their
    1σ uncertainties, recompute g_i and refit R50 each time, and return the
    spread of g_i, the kept fraction f, and the calibration shift ΔR50^wo. The
    shift band is the systematic that does NOT calibrate away."""
    t_irr, t1, t2 = w["timing"]["t_irr_s"], w["timing"]["t1_s"], w["timing"]["t2_s"]
    M0, Merr = np.array(w["model"]["M"]), np.array(w["model"]["M_err"])
    T0, Terr = np.array(w["model"]["T_s"]), np.array(w["model"]["T_s_err"])
    N = {iso: float(trapz(P[iso], z)) for iso in isotopes}
    Ntot = sum(N.values())
    rng = np.random.default_rng(MC_SEED)
    g_s = {iso: [] for iso in isotopes}
    dR50_s, f_s = [], []
    for _ in range(N_MC):
        M = np.clip(rng.normal(M0, Merr), 0.0, None)
        M = M / M.sum()                                  # renormalise Σ M = 1
        T = np.clip(rng.normal(T0, Terr), 0.1, None)     # positive half-lives
        mu = LN2 / T
        g = {iso: float(g_factor(lam[iso], M, mu, t_irr, t1, t2))
             for iso in isotopes}
        Pwo = sum(g[iso] * P[iso] for iso in isotopes)
        dR50_s.append(refit_r50(z, Pwo, window) - R50_phys)
        f_s.append(sum(g[iso] * N[iso] for iso in isotopes) / Ntot)
        for iso in isotopes:
            g_s[iso].append(g[iso])
    return ({iso: (float(np.mean(g_s[iso])), float(np.std(g_s[iso])))
             for iso in isotopes},
            np.array(dR50_s), np.array(f_s))


def main():
    w, window = load_params()
    t_irr = w["timing"]["t_irr_s"]
    t1, t2 = w["timing"]["t1_s"], w["timing"]["t2_s"]
    M = w["model"]["M"]
    mu_bio = [LN2 / T for T in w["model"]["T_s"]]
    isotopes = w["physical"]["isotopes"]
    lam = {iso: LN2 / T for iso, T in
           zip(isotopes, w["physical"]["T_half_s"])}

    # --- g_i, closed form + numeric cross-check
    g = {}
    print("isotope   T½[s]     g_i (closed)   g_i (numeric)   Δ")
    for iso in isotopes:
        gc = float(g_factor(lam[iso], M, mu_bio, t_irr, t1, t2))
        gn = float(g_factor_numeric(lam[iso], M, mu_bio, t_irr, t1, t2))
        g[iso] = gc
        assert abs(gc - gn) < 1e-6, f"{iso}: closed {gc} vs numeric {gn}"
        print(f"  {iso:4s}  {LN2/lam[iso]:8.1f}   {gc:9.4f}      "
              f"{gn:9.4f}   {gc-gn:+.2e}")

    # --- reweight the truth columns
    d = np.genfromtxt(os.path.join(TRUTH, "activity_profile_fast.csv"),
                      delimiter=",", names=True)
    z = np.asarray(d["z_mm"], float)
    P = {iso: np.asarray(d[iso], float) for iso in isotopes}
    P0 = sum(P.values())                       # physical-decay-only total
    Pwo = sum(g[iso] * P[iso] for iso in isotopes)   # washed total

    N = {iso: float(trapz(P[iso], z)) for iso in isotopes}
    Ntot = sum(N.values())
    f_keep = sum(g[iso] * N[iso] for iso in isotopes) / Ntot

    # --- refit the distal edge, physical vs washed (whole-plane erfc, free b)
    r0 = analyze("physical", z, P0, window, 0.0, True)["erfc"]
    rw = analyze("washed", z, Pwo, window, 0.0, True)["erfc"]
    dR50 = rw["R50_mm"] - r0["R50_mm"]

    # --- washout-parameter systematic (the band that does NOT calibrate away)
    g_band, dR50_mc, f_mc = mc_band(z, P, isotopes, lam, window, w, r0["R50_mm"])
    dR50_sys = float(np.std(dR50_mc))
    f_sys = float(np.std(f_mc))

    print(f"\nkept fraction f = {f_keep:.4f} ± {f_sys:.4f} (param systematic)")
    print(f"R50 physical = {r0['R50_mm']:.3f} mm | washed = {rw['R50_mm']:.3f} mm")
    print(f"ΔR50^wo (central) = {dR50:+.3f} mm  [calibrates away]")
    print(f"ΔR50^wo systematic band = ± {dR50_sys:.3f} mm "
          f"(from the washout parameters — this does NOT calibrate away; "
          f"compare to σ_R ≈ 0.11 mm)")

    os.makedirs(os.path.join(OUT, "figures"), exist_ok=True)
    toml_dump(os.path.join(OUT, "washout.toml"), {
        "meta": {"model": "Mizuno brain 3-exp, uniform; central + MC param band",
                 "n_mc": N_MC, "t_irr_s": t_irr, "t1_s": t1, "t2_s": t2,
                 "window_lo_mm": window[0], "window_hi_mm": window[1]},
        "g_factor": {iso: g[iso] for iso in isotopes},
        "g_factor_sys": {iso: g_band[iso][1] for iso in isotopes},
        "abundance_fraction": {iso: N[iso] / Ntot for iso in isotopes},
        "kept_fraction": f_keep, "kept_fraction_sys": f_sys,
        "R50_physical_mm": r0["R50_mm"], "R50_washed_mm": rw["R50_mm"],
        "delta_R50_washout_mm": dR50,
        "delta_R50_washout_sys_mm": dR50_sys,
        "w_physical_mm": r0["w_mm"], "w_washed_mm": rw["w_mm"]})
    print(f"wrote {os.path.join(OUT, 'washout.toml')}")

    plot_profiles(z, P0, Pwo, window, r0, rw,
                  os.path.join(OUT, "figures", "washout_profile.png"))
    plot_gfactors(isotopes, g, g_band, N, Ntot,
                  os.path.join(OUT, "figures", "g_factors.png"))


def plot_profiles(z, P0, Pwo, window, r0, rw, path):
    fig, ax = plt.subplots(figsize=(8.5, 5.0), facecolor=SURFACE)
    style(ax)
    ax.axvspan(*window, color=MUTED, alpha=0.12, lw=0)
    ax.plot(z, P0, "o", ms=3, mfc="none", mec=INK, mew=0.9,
            label="physical decay only")
    ax.plot(z, Pwo, "o", ms=3, mfc=RED, mec=RED, mew=0.9,
            label="with washout")
    for r, c in ((r0, INK), (rw, RED)):
        ax.axvline(r["R50_mm"], color=c, lw=1.0, ls="--", alpha=0.7)
    ax.set_xlabel("depth z [mm]", color=INK, fontsize=13)
    ax.set_ylabel("activity  A(z)  [decays]", color=INK, fontsize=13)
    ax.set_title(f"Truth activity edge: $\\Delta R_{{50}}^{{wo}}$ = "
                 f"{rw['R50_mm']-r0['R50_mm']:+.3f} mm", color=INK,
                 fontsize=12, loc="left")
    ax.legend(frameon=False, fontsize=12, labelcolor=INK, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


def plot_gfactors(isotopes, g, g_band, N, Ntot, path):
    fig, ax = plt.subplots(figsize=(7.0, 4.6), facecolor=SURFACE)
    style(ax)
    x = np.arange(len(isotopes))
    frac = [N[iso] / Ntot for iso in isotopes]
    err = [g_band[iso][1] for iso in isotopes]
    ax.bar(x, [g[iso] for iso in isotopes], width=0.6, color=BLUE, alpha=0.85,
           yerr=err, capsize=5, ecolor=INK)
    for i, iso in enumerate(isotopes):
        ax.text(i, g[iso] + err[i] + 0.02, f"{g[iso]:.2f}", ha="center",
                color=INK, fontsize=10)
        ax.text(i, 0.02, f"{100*frac[i]:.0f}%", ha="center", color="white",
                fontsize=9)
    ax.set_xticks(x, isotopes, color=INK, fontsize=12)
    ax.set_ylabel("survival factor $g_i$", color=INK, fontsize=13)
    ax.set_ylim(0, 1)
    ax.set_title("Per-isotope washout survival (bar) and abundance (label)",
                 color=INK, fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
