#!/usr/bin/env python3
"""
plot_endpoint_port.py

Render the endpoint-port cross-check figure (validation ladder rung 3): the
frozen Python erfc fit vs the CryspBrainSim Julia port, overlaid on the two
shared self-test datasets from test/data/endpoint_reference.npz.

The script is self-contained: it invokes Julia on this repo's project to fit
the shared arrays with the ported estimator, then draws both fits.

Run from anywhere:  python3 tools/plot_endpoint_port.py
Writes:             out/endpoint_port/julia_fits.npz
                    out/endpoint_port/figures/fit_crosscheck.png
"""
import pathlib
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.special import erfc  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
REF_PATH = REPO / "test" / "data" / "endpoint_reference.npz"
OUT = REPO / "out" / "endpoint_port"
JL_PATH = OUT / "julia_fits.npz"

# Categorical slots 1-2 of the validated dataviz reference palette, plus
# neutral inks; data points stay grey so the two fits carry the identity.
BLUE, AQUA = "#2a78d6", "#1baf7a"  # Python fit, Julia fit
INK, MUTED, GRIDC = "#1a1a19", "#8a897f", "#e8e7e2"
SURFACE = "#fcfcfb"

JULIA_SNIPPET = f"""
using CryspBrainSim, NPZ
ref = npzread("{REF_PATH}")
fA = fit_endpoint(ref["zA"], ref["profA"];
                  window=(ref["winA"][1], ref["winA"][2]), weighted=true)
img = Float64.(ref["image"])
zB, profB = depth_profile(img; voxel_size_mm=Tuple(ref["image_vox"]),
                          beam_axis=3, roi_radius_mm=ref["image_roi_radius"][1])
fB = fit_endpoint(zB, profB; window=(ref["winB"][1], ref["winB"][2]),
                  weighted=true)
npzwrite("{JL_PATH}",
    Dict("A_popt" => fA.popt, "A_z0_err" => [fA.z0_err],
         "B_popt" => fB.popt, "B_z0_err" => [fB.z0_err]))
"""


def julia_fits() -> np.lib.npyio.NpzFile:
    OUT.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["julia", f"--project={REPO}", "-e", JULIA_SNIPPET], check=True
    )
    return np.load(JL_PATH)


def model(z, p):
    return p[0] + p[1] * 0.5 * erfc((z - p[2]) / (np.sqrt(2.0) * p[3]))


def panel(ax, z, prof, win, py_popt, jl_popt, py_err, jl_err, title):
    zf = np.linspace(z[0], z[-1], 1200)
    ax.axvspan(*win, color=GRIDC, alpha=0.5, lw=0, zorder=0)
    ax.plot(z, prof, "o", ms=2.5, color=MUTED, mec="none", zorder=2,
            label="profile counts")
    ax.plot(zf, model(zf, py_popt), "-", lw=2, color=BLUE, zorder=3,
            label=f"Python fit  z0 = {py_popt[2]:.2f} ± {py_err:.2f} mm")
    ax.plot(zf, model(zf, jl_popt), "--", lw=2, color=AQUA, zorder=4,
            label=f"Julia fit   z0 = {jl_popt[2]:.2f} ± {jl_err:.2f} mm")
    ax.set_title(title, color=INK, fontsize=11, loc="left")
    ax.set_xlabel("depth z (mm)", color=INK)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK)
    ax.grid(axis="y", color=GRIDC, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8.5)
    ax.set_xlim(win[0] - 45, win[1] + 15)


def main():
    ref = np.load(REF_PATH)
    jl = julia_fits()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)
    for ax in (axA, axB):
        ax.set_facecolor(SURFACE)

    panel(axA, ref["zA"], ref["profA"], ref["winA"],
          ref["A_weighted_popt"], jl["A_popt"],
          ref["A_weighted_z0_err"][0], jl["A_z0_err"][0],
          "Test A — 1-D erfc edge + Poisson noise (true R50 = 150 mm)")
    axA.set_ylabel("counts per bin", color=INK)

    panel(axB, ref["zB"], ref["profB"], ref["winB"],
          ref["B_weighted_popt"], jl["B_popt"],
          ref["B_weighted_z0_err"][0], jl["B_z0_err"][0],
          "Test B — image → disc-ROI profile (true R50 = 150 mm)")
    axB.set_ylabel("ROI-summed counts per slice", color=INK)

    fig.suptitle("Endpoint estimator port: scipy reference vs Julia (shared "
                 "arrays; grey band = fixed distal fit window)", color=INK,
                 fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    figdir = OUT / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    path = figdir / "fit_crosscheck.png"
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
