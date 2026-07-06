# Shared setup and the reconstruction chain for the σ_R drivers
# (sigma_r_at_dose.jl and sigma_r_sweep_dose.jl). Both include this file.
#
# Vocabulary used across both drivers and their comments:
#
#   coincidence (LOR)  one recorded event — a pair of 511 keV annihilation
#                      photons caught in time coincidence, which defines a line
#                      of response (LOR) through the head.
#   shard              one complete simulated acquisition: every coincidence the
#                      scanner records for the proton field at the nominal (top)
#                      dose. Ten shards exist, each an independent repeat of the
#                      same source with a different random seed.
#   pooled master      the ten shards concatenated — the largest count sample
#                      available, and the source each realization is drawn from.
#   dose               the delivered radiation dose, in Gy. The β⁺ activity, and
#                      therefore the coincidence count, is proportional to it.
#   realization        one emulated acquisition, built by randomly keeping a
#                      fraction of the pooled coincidences. A valid stand-in for
#                      a fresh acquisition that needs no new simulation; the kept
#                      fraction sets the emulated dose.
#   R (R50)            the range endpoint — the depth of the distal activity
#                      edge, read from the reconstructed depth profile.
#   σ_R                the standard deviation of R across realizations: the
#                      statistical precision of the range measurement at one dose.
#
# The reconstruction chain is frozen in config/run_parameters.toml and identical
# for every run: keep the true coincidences → an attenuation factor per LOR →
# MLEM at the frozen iteration count → sum the image over the fixed disc ROI at
# each depth → fit the distal edge with an erfc. R is read two ways — the fit's
# half-height (the primary, noise-averaged value) and the profile's half-height
# crossing inside the window (a direct read, which is noisier).

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf
using Statistics: mean
using TOML

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const SCEN = joinpath(ROOT, "uniform_headep_sobp_1e8")
const OUT  = joinpath(dirname(@__DIR__), "out", "sigma_r")

const PARAMS = load_run_parameters()
const N   = PARAMS.grid.n
const VS  = PARAMS.grid.voxsize
const ORG = PARAMS.grid.img_origin
const SENS_CACHE = joinpath(dirname(@__DIR__), "out", "sensitivity",
    "crysp_ring_1m_grid64x64x96_1.5mm_orgm47.25_m47.25_m119.25_n$(PARAMS.n_sens)")

# A short filename tag for a dose, e.g. 1.0 → "1Gy", 0.5 → "0p5Gy".
dose_tag(d) = replace((@sprintf "%.3gGy" d), "." => "p")

"""
    setup() -> NamedTuple

Load the inputs shared by every run: the truth reference (dose-R80, the true
activity edge, the fixed fit window), the sensitivity image (checked to match
the frozen grid), the head's attenuation parameters, the compute device (Metal
GPU when available), and the ten shard files.
"""
function setup()
    ref = characterize(SCEN)
    base, meta = load_sensitivity(SENS_CACHE)
    g = meta["grid"]
    (g["n"] == collect(N) && Float32.(g["img_origin"]) == collect(ORG) &&
     Float32.(g["voxsize"]) == collect(VS)) ||
        error("sensitivity cache grid ≠ frozen run-parameter grid — rebuild the cache")
    ph = phantom_attenuation(SCEN)
    dev = Metal.functional() ? MtlArray : identity
    leaf = leaf_dir(ROOT; scenario="uniform_headep_sobp_1e8",
                    scanner="crysp_ring_1m", crystal="bgo", leaf="fast_1Gy")
    return (ref=ref, base=base, meta=meta, ph=ph, dev=dev,
            files=shard_files(leaf))
end

"Per-LOR attenuation factors exp(-μ · chord) through the uniform head ellipsoid."
attenuation(xs, xe, ph) =
    attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes, centre=ph.centre,
                          mu_mm_inv=ph.mu_mm_inv)

"""
    recon_endpoint(xs, xe, mult, ref, base, dev) -> NamedTuple

Run the frozen reconstruction chain on one list of coincidences (`xs`/`xe`
endpoints, `mult` the per-LOR attenuation) and read the range endpoint. Scales
the sensitivity to this list's count, runs MLEM, sums the image over the disc
ROI into a depth profile, and fits the distal edge. Returns the count, R in
both conventions (`r50_fit`, `r50_cross`), the fitted edge width, and the fit's
own error estimate.
"""
function recon_endpoint(xs, xe, mult, ref, base, dev)
    nev = size(xs, 2)
    sens = scaled_sensitivity(base, nev, PARAMS.n_sens)
    model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                 img_origin=ORG, voxsize=VS, mult=dev(mult))
    x = mlem(model, dev(Float32.(sens .> 0)); niter=PARAMS.niter)
    z, prof = depth_profile(Array(x); voxel_size_mm=VS, beam_axis=3,
                            roi_radius_mm=PARAMS.roi.radius_mm, z_origin_mm=ORG[3])
    fit = fit_endpoint(z, prof; window=ref.window, weighted=true)
    return (nev=nev, r50_fit=fit.z0, z0_err=fit.z0_err, w=fit.w,
            r50_cross=windowed_crossing(z, prof, ref.window))
end

"Print the two-convention σ_R summary for a set of realizations."
function report_sigma(label, sf, sc)
    @printf("\nσ_R %s (n = %d, the spread itself is known to ±%.0f%%):\n",
            label, sf.n_ok, 100 / sqrt(2 * (sf.n_ok - 1)))
    @printf("  fit read:      mean R50 %8.3f mm | σ_R %.3f mm | offset to dose-R80 %.3f mm\n",
            sf.mean, sf.sigma, sf.offset)
    @printf("  crossing read: mean R50 %8.3f mm | σ_R %.3f mm\n", sc.mean, sc.sigma)
end
