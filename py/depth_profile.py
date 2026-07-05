"""
depth_profile.py

Steps 4-5 of the CRYSP range-verification pipeline (see range_endpoint.py for
steps 2, 6). This module supersedes the fit_endpoint in range_endpoint.py:
it adds a FIXED distal analysis window and Poisson (sqrt-counts) weighting.

    reconstructed image  --depth_profile-->  (z, P(z))       [step 4]
    (z, P(z))            --fit_endpoint-->   R50 (+R80,R20)   [step 5]

Design invariants (identical across the three geometries):
  * beam-aligned grid: the beam/depth axis is a grid axis, so profiling is a
    plain transverse sum with no oblique resampling;
  * the transverse ROI is a fixed disc centred on the BEAM AXIS, never on the
    per-realisation activity centroid (a data-driven centre would inject a noisy
    shift into the endpoint and inflate sigma_R);
  * integrate (sum), do not average or normalise --- this preserves the count
    statistics the endpoint fit and sigma_R depend on;
  * the fit window is fixed from the NOMINAL range, not found per realisation.

Honest caveat on weighting: after MLEM the voxels in a slice are correlated, so
P(z) is only approximately Poisson. The sqrt-counts weighting is therefore a
practical heuristic that improves the fit in the sparse tail; the authoritative
uncertainty on the endpoint is the ensemble spread sigma_R (step 6), NOT the
per-realisation fit error returned here.

British spelling in comments for consistency with the wider project.
"""

from __future__ import annotations
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc, erfcinv


# ----------------------------------------------------------------------
# Step 4: 1-D depth-activity profile from a beam-aligned image
# ----------------------------------------------------------------------
def depth_profile(image, voxel_size_mm, beam_axis=2,
                  roi_radius_mm=None, roi_centre_mm=(0.0, 0.0),
                  z_origin_mm=0.0):
    """
    Sum a beam-aligned reconstructed image over a fixed transverse disc ROI at
    each depth slice.

    Parameters
    ----------
    image : ndarray (3-D)
        Reconstructed activity, already resampled onto a beam-aligned grid so
        that `beam_axis` is the depth (z) direction.
    voxel_size_mm : sequence of 3 floats
        Voxel edge lengths (mm) along the three image axes.
    beam_axis : int
        Which image axis is the beam/depth axis (0, 1 or 2).
    roi_radius_mm : float or None
        Radius of the transverse disc ROI (mm). If None, the whole transverse
        plane is summed (NOT recommended --- dilutes the distal contrast with
        noise; a warning-worthy default kept only for quick looks).
    roi_centre_mm : (float, float)
        Transverse centre of the ROI relative to the grid centre, in the order
        of the two non-beam axes. Keep at (0, 0) to sit on the beam axis.
    z_origin_mm : float
        Depth coordinate assigned to the first slice (voxel-centre convention).

    Returns
    -------
    z : 1-D ndarray
        Depth coordinate of each slice (mm).
    prof : 1-D ndarray
        Integrated (summed) activity counts within the ROI at each depth.

    Notes
    -----
    Same `roi_radius_mm`, `roi_centre_mm` and grid MUST be used for every
    geometry. In particular the ROI must be generous enough not to clip the
    limited-angle elongation of the open dual-head arm differently from the
    closed arm --- the aim is to MEASURE that elongation's effect on the
    endpoint, not to hide it behind an ROI wall. Rule of thumb:
        roi_radius_mm  >~  n * sqrt(sigma_spot^2 + sigma_PSF^2) + R_positron
    with n ~ 3.
    """
    image = np.asarray(image, float)
    if image.ndim != 3:
        raise ValueError("depth_profile expects a 3-D image")
    if len(voxel_size_mm) != 3:
        raise ValueError("voxel_size_mm must have 3 entries")

    # Bring the beam axis to the front: img has shape (nz, nu, nv).
    img = np.moveaxis(image, beam_axis, 0)
    nz, nu, nv = img.shape

    other_axes = [a for a in (0, 1, 2) if a != beam_axis]
    dz = float(voxel_size_mm[beam_axis])
    du = float(voxel_size_mm[other_axes[0]])
    dv = float(voxel_size_mm[other_axes[1]])

    # Transverse voxel-centre coordinates relative to the grid centre.
    u = (np.arange(nu) - (nu - 1) / 2.0) * du
    v = (np.arange(nv) - (nv - 1) / 2.0) * dv
    U, V = np.meshgrid(u, v, indexing="ij")

    cu, cv = roi_centre_mm
    if roi_radius_mm is None:
        mask = np.ones((nu, nv), dtype=bool)
    else:
        mask = ((U - cu) ** 2 + (V - cv) ** 2) <= float(roi_radius_mm) ** 2

    # Integrate (sum) over the ROI at each depth --- preserves Poisson statistics.
    prof = (img * mask[None, :, :]).sum(axis=(1, 2))
    z = z_origin_mm + np.arange(nz) * dz
    return z, prof


# ----------------------------------------------------------------------
# Helper: fixed distal analysis window from the nominal range
# ----------------------------------------------------------------------
def distal_window(z_edge_nominal_mm, proximal_margin_mm=20.0,
                  distal_margin_mm=15.0):
    """
    Build the FIXED fit window bracketing the distal falloff. Defined from the
    nominal (expected) activity edge, identical across arms and realisations ---
    never from a per-realisation peak search, which would correlate the window
    with the noise being measured.

    Returns (z_lo, z_hi): start on the plateau proximal to the edge, end in the
    tail distal to it.
    """
    return (z_edge_nominal_mm - proximal_margin_mm,
            z_edge_nominal_mm + distal_margin_mm)


# ----------------------------------------------------------------------
# Step 5: endpoint by windowed, Poisson-weighted erfc fit of the falloff
# ----------------------------------------------------------------------
def _erfc_model(z, base, amp, z0, w):
    """Monotonic distal falloff; 0.5*erfc(0)=0.5 so z0 is the R50 point."""
    return base + amp * 0.5 * erfc((z - z0) / (np.sqrt(2.0) * w))


def fit_endpoint(z, profile, window, levels=(0.5, 0.8, 0.2),
                 weighted=True, p0=None):
    """
    Fit the distal falloff inside a FIXED window and read the R-levels from the
    fit. Extends the range_endpoint.py version with the window and weighting.

    Parameters
    ----------
    z, profile : 1-D arrays
        Output of depth_profile (proximal->distal in increasing z).
    window : (z_lo, z_hi)
        Fixed distal analysis window (e.g. from distal_window()). The proximal
        rise and plateau outside the edge region are excluded so the single
        erfc edge is an appropriate model.
    levels : tuple of float
        Falloff fractions of the fitted plateau to report (0.5 primary).
    weighted : bool
        If True, weight each point by its Poisson sigma = sqrt(counts) (floored
        at 1) with absolute_sigma=True, so the returned z0_err reflects counting
        statistics. See module docstring for the correlation caveat.
    p0 : tuple or None
        Initial (base, amp, z0, w); auto-estimated if None.

    Returns
    -------
    dict: R (level->z_mm), z0, w, z0_err, n_points, popt, pcov.
    Endpoints are np.nan on fit failure or too few points (caller should count
    and exclude; a high failure rate at a dose level is itself a result).
    """
    z = np.asarray(z, float)
    profile = np.asarray(profile, float)

    zlo, zhi = window
    sel = (z >= zlo) & (z <= zhi)
    zf, pf = z[sel], profile[sel]
    nan = dict(R={lv: np.nan for lv in levels}, z0=np.nan, w=np.nan,
               z0_err=np.nan, n_points=int(sel.sum()), popt=None, pcov=None)
    if zf.size < 4:
        return nan

    if p0 is None:
        base0 = np.median(pf[-max(2, pf.size // 5):])      # distal tail
        plateau0 = np.median(pf[:max(2, pf.size // 5)])    # proximal plateau
        amp0 = max(plateau0 - base0, 1e-9)
        half = base0 + 0.5 * amp0
        above = np.where(pf >= half)[0]
        z0_0 = zf[above[-1]] if above.size else zf[zf.size // 2]
        w0 = 0.1 * (zf[-1] - zf[0]) + 1e-6
        p0 = (base0, amp0, z0_0, w0)

    sigma = np.sqrt(np.clip(pf, 1.0, None)) if weighted else None
    try:
        popt, pcov = curve_fit(_erfc_model, zf, pf, p0=p0, sigma=sigma,
                               absolute_sigma=bool(weighted), maxfev=20000)
        base, amp, z0, w = popt
        perr = np.sqrt(np.clip(np.diag(pcov), 0, np.inf))
        R = {lv: z0 + np.sqrt(2.0) * w * erfcinv(2.0 * lv) for lv in levels}
        return dict(R=R, z0=float(z0), w=abs(float(w)), z0_err=float(perr[2]),
                    n_points=int(zf.size), popt=popt, pcov=pcov)
    except Exception:
        return nan


# ----------------------------------------------------------------------
# Self-test: synthetic beam-aligned image -> profile -> endpoint
# ----------------------------------------------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(1)

    # Grid: 1.5 mm isotropic voxels; z = depth axis (axis 2).
    nx, ny, nz = 41, 41, 160
    vox = (1.5, 1.5, 1.5)
    x = (np.arange(nx) - (nx - 1) / 2.0) * vox[0]
    y = (np.arange(ny) - (ny - 1) / 2.0) * vox[1]
    z = np.arange(nz) * vox[2]
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Transverse Gaussian beam (sigma ~2 mm) x depth activity: plateau then
    # erfc falloff with true R50 = 150 mm.
    sigma_perp = 2.0
    transverse = np.exp(-(X**2 + Y**2) / (2 * sigma_perp**2))
    z0_true = 150.0
    depth = 0.5 * erfc((Z - z0_true) / (np.sqrt(2.0) * 3.0))       # falloff
    depth *= (1.0 + 0.15 * np.tanh((Z - 20.0) / 15.0))            # gentle rise
    mean_img = 40.0 * transverse * depth
    noisy_img = rng.poisson(np.clip(mean_img, 0, None)).astype(float)

    zc, prof = depth_profile(noisy_img, voxel_size_mm=vox, beam_axis=2,
                             roi_radius_mm=3 * sigma_perp + 2.0)
    win = distal_window(z0_true, proximal_margin_mm=25, distal_margin_mm=20)
    out = fit_endpoint(zc, prof, window=win, weighted=True)

    print(f"true R50 = {z0_true:.1f} mm | fitted z0 = {out['z0']:.2f} "
          f"+/- {out['z0_err']:.2f} mm  (window pts: {out['n_points']})")
    print(f"R80 = {out['R'][0.8]:.2f} mm | R20 = {out['R'][0.2]:.2f} mm | "
          f"edge width w = {out['w']:.2f} mm")
