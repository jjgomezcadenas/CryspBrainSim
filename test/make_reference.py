#!/usr/bin/env python3
"""
make_reference.py

Generate the shared cross-validation data for the endpoint port (dev/PLAN.md
validation ladder rung 3): the two py self-test datasets plus the frozen
Python fit outputs, all in one NPZ that `test/runtests.jl` reads via NPZ.jl.

Run from anywhere:  python3 test/make_reference.py
Writes:             test/data/endpoint_reference.npz  (committed)

Re-run only after a deliberate change to the frozen estimator in py/ — the
committed NPZ is the reference the Julia tests compare against, so the Julia
test suite runs without a Python environment.
"""
from __future__ import annotations
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "py"))

from depth_profile import depth_profile, distal_window, fit_endpoint  # noqa: E402
from range_endpoint import _erfc_model, sigma_R  # noqa: E402

out: dict[str, np.ndarray] = {}
LEVELS = (0.5, 0.8, 0.2)


def record_fit(prefix: str, fit: dict) -> None:
    out[f"{prefix}_popt"] = np.asarray(fit["popt"], float)
    out[f"{prefix}_pcov_diag"] = np.diag(fit["pcov"]).astype(float)
    out[f"{prefix}_z0_err"] = np.array([fit["z0_err"]])
    out[f"{prefix}_R"] = np.array([fit["R"][lv] for lv in LEVELS])
    out[f"{prefix}_n_points"] = np.array([fit["n_points"]])


# ---------------------------------------------------------------------------
# Test A — 1-D erfc edge with Poisson noise (the range_endpoint.py self-test
# data), fitted with the windowed estimator from depth_profile.py in both
# weighted (absolute-sigma) and unweighted (MSE-scaled) covariance modes.
# ---------------------------------------------------------------------------
rng = np.random.default_rng(0)
zA = np.linspace(0.0, 200.0, 400)
truthA = _erfc_model(zA, base=2.0, amp=100.0, z0=150.0, w=3.0)
profA = rng.poisson(truthA).astype(float)
winA = distal_window(150.0, proximal_margin_mm=25.0, distal_margin_mm=20.0)

out["zA"] = zA
out["profA"] = profA
out["winA"] = np.asarray(winA)
record_fit("A_weighted", fit_endpoint(zA, profA, window=winA, weighted=True))
record_fit("A_unweighted", fit_endpoint(zA, profA, window=winA, weighted=False))

# ---------------------------------------------------------------------------
# Test B — synthetic beam-aligned 3-D image (the depth_profile.py self-test):
# Gaussian pencil beam × (gentle rise + erfc falloff at R50 = 150 mm), Poisson
# noise, profiled over the fixed disc ROI and fitted in the fixed window.
# ---------------------------------------------------------------------------
rng = np.random.default_rng(1)
nx, ny, nz = 41, 41, 160
vox = (1.5, 1.5, 1.5)
x = (np.arange(nx) - (nx - 1) / 2.0) * vox[0]
y = (np.arange(ny) - (ny - 1) / 2.0) * vox[1]
z = np.arange(nz) * vox[2]
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

from scipy.special import erfc  # noqa: E402

sigma_perp = 2.0
transverse = np.exp(-(X**2 + Y**2) / (2 * sigma_perp**2))
z0_true = 150.0
depth = 0.5 * erfc((Z - z0_true) / (np.sqrt(2.0) * 3.0))
depth *= 1.0 + 0.15 * np.tanh((Z - 20.0) / 15.0)
mean_img = 40.0 * transverse * depth
image = rng.poisson(np.clip(mean_img, 0, None)).astype(np.int16)

roi_radius = 3 * sigma_perp + 2.0
zB, profB = depth_profile(
    image.astype(float), voxel_size_mm=vox, beam_axis=2, roi_radius_mm=roi_radius
)
winB = distal_window(z0_true, proximal_margin_mm=25.0, distal_margin_mm=20.0)

out["image"] = image
out["image_vox"] = np.asarray(vox)
out["image_roi_radius"] = np.array([roi_radius])
out["zB"] = zB
out["profB"] = profB
out["winB"] = np.asarray(winB)
record_fit("B_weighted", fit_endpoint(zB, profB, window=winB, weighted=True))

# ---------------------------------------------------------------------------
# sigma_R — deterministic endpoint set with two failed fits (NaN) and a
# dose-R80 offset, collapsed by the py aggregator.
# ---------------------------------------------------------------------------
endpoints = np.array([150.1, 149.7, np.nan, 150.4, 149.9, np.nan, 150.2, 149.8])
res = sigma_R(endpoints, dose_bragg_peak=152.0)
out["sig_endpoints"] = endpoints
out["sig_dose_bragg_peak"] = np.array([152.0])
out["sig_ref"] = np.array(
    [res["n_ok"], res["n_fail"], res["mean"], res["sigma"], res["sem"], res["offset"]]
)

data_dir = HERE / "data"
data_dir.mkdir(exist_ok=True)
path = data_dir / "endpoint_reference.npz"
np.savez_compressed(path, **out)

fa, fb = out["A_weighted_popt"], out["B_weighted_popt"]
print(f"wrote {path}")
print(
    f"A weighted : z0 = {fa[2]:.4f} +/- {out['A_weighted_z0_err'][0]:.4f} mm, "
    f"w = {fa[3]:.4f} mm, n = {out['A_weighted_n_points'][0]}"
)
print(
    f"B weighted : z0 = {fb[2]:.4f} +/- {out['B_weighted_z0_err'][0]:.4f} mm, "
    f"w = {fb[3]:.4f} mm, n = {out['B_weighted_n_points'][0]}"
)
print(
    f"sigma_R    : mean = {res['mean']:.4f}, sigma = {res['sigma']:.4f}, "
    f"sem = {res['sem']:.4f}, offset = {res['offset']:.4f}, "
    f"n_ok/n_fail = {res['n_ok']}/{res['n_fail']}"
)
