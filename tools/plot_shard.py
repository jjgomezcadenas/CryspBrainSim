#!/usr/bin/env python3
"""
plot_shard.py — the 3×3 detector QA panel for one stored LOR shard
(validation ladder rung 2; the figure companion of src/qa.jl's ShardQA).
Adapted from PTCryspMC py/plot_prod.py for the PtCryspProds tree: the input
is a lors_shardNNN.h5 config-leaf shard, the scanner geometry is found by
walking up to the scanner's scanner_geometry.json, and the PNG lands under
this repo's out/ tree.

Panels: energy spectra e1/e2 · e1 vs e2 · energy by truth · ring hit map ·
LOR axial midpoint · source fill (transverse, axial) · coincidence Δt ·
composition + acceptance.

Run:  python3 tools/plot_shard.py --shard <.../fast_1Gy/lors_shard000.h5>
Writes out/shard_qa/figures/<scenario>_<crystal>_<leaf>_shardNNN.png
(override with --out).
"""
import argparse
import json
import os

import h5py
import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Truth-class colors, consistent across every panel (colorblind-safe Tableau,
# as in the upstream plotter).
C_TRUE, C_SCAT, C_RAND = "#2ca02c", "#ff7f0e", "#d62728"
C_SINGLE, C_MULTI = "#fdae6b", "#e6550d"  # single vs multiple scatter
TRUTH_TRUE, TRUTH_SCATTER, TRUTH_RANDOM = 0, 1, 2


def find_geometry(shard_path):
    """Walk up from the shard to the first scanner_geometry.json (two levels
    for a homogeneous scanner, one for a heterogeneous one)."""
    d = os.path.dirname(os.path.abspath(shard_path))
    for _ in range(4):
        p = os.path.join(d, "scanner_geometry.json")
        if os.path.isfile(p):
            return p
        d = os.path.dirname(d)
    return None


def scanner_dims(geom_file):
    """(r_inner_mm, wall_mm, half_length_mm) from the geometry JSON."""
    with open(geom_file) as f:
        s = json.load(f)["scanner"]
    return 10.0 * s["r_inner_cm"], 10.0 * s["wall_thickness_cm"], 10.0 * s["half_length_cm"]


def load_shard(path):
    """Shard columns as float arrays (de-quantized via the stored scales)
    plus the provenance attrs."""
    with h5py.File(path, "r") as f:
        a = dict(f.attrs)
        xs = float(a["xyz_scale_mm"])
        es = float(a["e_scale_keV"])
        col = lambda k, s=1.0: f[k][:].astype(np.float64) * s
        d = {}
        for k in ("x1_mm", "y1_mm", "z1_mm", "x2_mm", "y2_mm", "z2_mm",
                  "x0_mm", "y0_mm", "z0_mm"):
            d[k] = col(k, xs)
        for k in ("e1_keV", "e2_keV"):
            d[k] = col(k, es)
        for k in ("t1_ns", "t2_ns", "dt_ns"):
            d[k] = col(k)
        for k in ("truth", "nscat1", "nscat2"):
            d[k] = f[k][:].astype(np.int64)
    return d, a


def dec(v):
    return v.decode() if isinstance(v, bytes) else v


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--shard", required=True, help="path to lors_shardNNN.h5")
    p.add_argument("--geometry", default=None,
                   help="scanner_geometry.json (default: found above the shard)")
    p.add_argument("--out", default=None, help="output PNG path")
    args = p.parse_args()

    geom = args.geometry or find_geometry(args.shard)
    Ri, wall, H = scanner_dims(geom) if geom else (387.0, 37.0, 512.0)

    d, attrs = load_shard(args.shard)
    scenario = dec(attrs.get("scenario", "?"))
    crystal = dec(attrs.get("crystal", "?"))
    shard_ix = int(attrs.get("realization", -1))
    leaf = os.path.basename(os.path.dirname(os.path.abspath(args.shard)))
    if not args.out:
        name = f"{scenario}_{crystal.lower()}_{leaf}_shard{shard_ix:03d}.png"
        args.out = os.path.join(REPO, "out", "shard_qa", "figures", name)

    truth = d["truth"]
    n = truth.size
    is_t = truth == TRUTH_TRUE
    is_s = truth == TRUTH_SCATTER
    is_r = truth == TRUTH_RANDOM
    n_t, n_s, n_r = int(is_t.sum()), int(is_s.sum()), int(is_r.sum())
    nev = int(attrs.get("nevents", 0))
    eres = float(attrs.get("eres", 0.0))
    emin = float(attrs.get("emin_keV", 0.0))
    tau = float(attrs.get("tau_ns", 0.0))
    sig = float(attrs.get("sigma_xyz_mm", 0.0))

    # both-hit arrays (per-hit panels), tagged with each pair's truth
    e_all = np.concatenate([d["e1_keV"], d["e2_keV"]])
    x_all = np.concatenate([d["x1_mm"], d["x2_mm"]])
    y_all = np.concatenate([d["y1_mm"], d["y2_mm"]])
    z_all = np.concatenate([d["z1_mm"], d["z2_mm"]])
    truth_all = np.concatenate([truth, truth])
    phi_all = np.degrees(np.arctan2(y_all, x_all) % (2 * np.pi))
    zmid = 0.5 * (d["z1_mm"] + d["z2_mm"])  # LOR axial midpoint
    nsc = d["nscat1"] + d["nscat2"]         # total scatter multiplicity

    fig = plt.figure(figsize=(15, 13))
    fig.suptitle(f"{scenario} / {crystal} / {leaf} — shard {shard_ix}: {n:,} LORs "
                 f"(eres {eres:.0%}, Emin {emin:.0f} keV, τ {tau:g} ns)", fontsize=13)
    ax = lambda i: fig.add_subplot(3, 3, i)

    # 1. energy spectra e1, e2 + cut + 511
    a = ax(1)
    a.hist(d["e1_keV"], bins=80, range=(0, 560), histtype="step", color="C0", label="e1")
    a.hist(d["e2_keV"], bins=80, range=(0, 560), histtype="step", color="C1", label="e2")
    a.axvline(511, color="gray", ls="--", lw=1)
    if emin:
        a.axvline(emin, color=C_RAND, ls=":", lw=1.2, label=f"Emin {emin:.0f}")
    a.set_yscale("log"); a.set_xlabel("energy [keV]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Energy spectra (e1, e2)")

    # 2. e1 vs e2
    a = ax(2)
    h = a.hist2d(d["e1_keV"], d["e2_keV"], bins=70, range=[[300, 560], [300, 560]],
                 cmap="viridis")
    fig.colorbar(h[3], ax=a, label="LORs")
    a.set_xlabel("e1 [keV]"); a.set_ylabel("e2 [keV]"); a.set_aspect("equal")
    a.set_title("e1 vs e2 (photopeak 511,511)")

    # 3. energy by truth class
    a = ax(3)
    for mask, c, lab in ((truth_all == TRUTH_TRUE, C_TRUE, "true"),
                         (truth_all == TRUTH_SCATTER, C_SCAT, "scatter"),
                         (truth_all == TRUTH_RANDOM, C_RAND, "random")):
        if mask.any():
            a.hist(e_all[mask], bins=80, range=(0, 560), histtype="step", color=c, label=lab)
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_yscale("log"); a.set_xlabel("energy [keV]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Energy by truth")

    # 4. ring hit map (unrolled φ–z)
    a = ax(4)
    h = a.hist2d(phi_all, z_all, bins=[72, 40], range=[[0, 360], [-H, H]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="hits")
    a.set_xlabel("φ [deg]"); a.set_ylabel("z [mm]")
    a.set_title("Ring hit map (unrolled)")

    # 5. LOR axial midpoint — for a compact central source this traces the
    #    SOURCE axial extent, not the detector sensitivity.
    a = ax(5)
    a.hist(zmid, bins=60, range=(-H, H), histtype="step", color="C0", label="all")
    if n_t:
        a.hist(zmid[is_t], bins=60, range=(-H, H), histtype="step", color=C_TRUE, label="true")
    a.set_xlabel("LOR axial midpoint z [mm]"); a.set_ylabel("LORs"); a.legend()
    a.set_title("LOR axial midpoint (≈ source axial extent)")

    # 6. source fill, transverse
    a = ax(6)
    lim = 1.05 * max(np.abs(d["x0_mm"]).max(), np.abs(d["y0_mm"]).max(), 1.0)
    h = a.hist2d(d["x0_mm"], d["y0_mm"], bins=60, range=[[-lim, lim], [-lim, lim]],
                 cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("y0 [mm]"); a.set_aspect("equal")
    a.set_title("Source fill (transverse)")

    # 7. source fill, axial
    a = ax(7)
    limz = 1.05 * max(np.abs(d["z0_mm"]).max(), 1.0)
    h = a.hist2d(d["x0_mm"], d["z0_mm"], bins=60, range=[[-lim, lim], [-limz, limz]],
                 cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("z0 [mm]")
    a.set_title("Source fill (axial)")

    # 8. coincidence Δt = t1 − t2 (raw, the τ-window quantity); the
    #    TOF-corrected residual annotated as one number.
    a = ax(8)
    draw = d["t1_ns"] - d["t2_ns"]
    xr = tau if tau > 0 else float(np.nanpercentile(np.abs(draw), 99.5) or 1.0)
    for mask, c, lab in ((is_t, C_TRUE, "true"), (is_s, C_SCAT, "scatter"),
                         (is_r, C_RAND, "random")):
        if mask.any():
            a.hist(draw[mask], bins=70, range=(-xr, xr), histtype="step", color=c, label=lab)
    a.axvline(0, color="gray", ls="--", lw=1)
    for s in (-tau, tau):
        if tau > 0:
            a.axvline(s, color="gray", ls=":", lw=1)
    a.set_yscale("log"); a.set_xlabel("Δt = t1 − t2 [ns]"); a.set_ylabel("LORs"); a.legend()
    a.set_title("Coincidence Δt (raw, τ window)")
    dt = d["dt_ns"]
    fin = np.isfinite(dt)
    med = np.median(np.abs(dt[is_t & fin])) if (is_t & fin).any() else float("nan")
    a.text(0.03, 0.97, f"TOF-corr.\nmedian|dt|\n{med:.3f} ns", transform=a.transAxes,
           ha="left", va="top", fontsize=8, family="monospace",
           bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    # 9. composition + acceptance
    a = ax(9)
    n_single = int(((nsc == 1) & is_s).sum())
    n_multi = int(((nsc >= 2) & is_s).sum())
    a.bar("true", n_t, color=C_TRUE)
    a.bar("scatter", n_single, color=C_SINGLE, label="single")
    a.bar("scatter", n_multi, bottom=n_single, color=C_MULTI, label="multiple")
    a.bar("random", n_r, color=C_RAND)
    a.set_ylabel("LORs"); a.legend(loc="upper right", fontsize=8)
    a.set_title("Composition")
    acc = n / nev if nev else float("nan")
    txt = (f"N annihilations = {nev:,}\n"
           f"n LORs = {n:,}\n"
           f"acceptance = {acc:.2%}\n"
           f"true    {n_t/max(n,1):.1%}\n"
           f"scatter {n_s/max(n,1):.1%}  (S {n_single/max(n,1):.1%} / M {n_multi/max(n,1):.1%})\n"
           f"random  {n_r/max(n,1):.1%}\n"
           f"R/(T+S) = {n_r/max(n_t+n_s,1):.2%}\n"
           f"σ_xyz {sig:g} mm")
    a.text(0.97, 0.97, txt, transform=a.transAxes, ha="right", va="top", fontsize=8,
           family="monospace", bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")
    print(f"  {n:,} LORs;  acceptance {acc:.2%};  "
          f"true {n_t/max(n,1):.1%} / scatter {n_s/max(n,1):.1%} / random {n_r/max(n,1):.1%}")


if __name__ == "__main__":
    main()
