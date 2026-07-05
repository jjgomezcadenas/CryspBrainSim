# endpoint.jl — steps 5–6 of the range-verification pipeline: fit the distal
# falloff of a depth profile inside the fixed window (windowed,
# Poisson-weighted erfc fit) and collapse per-realisation endpoints into
# (mean offset, σ_R). Julia port of py/depth_profile.py (fit_endpoint) and
# py/range_endpoint.py (sigma_R); the frozen Python modules are the
# cross-validation reference in test/runtests.jl.
#
# Covariance convention (matches scipy `absolute_sigma=True`): each point is
# weighted by its Poisson sigma √max(P, 1), passed to LsqFit as
# `PrecisionWeights(1 ./ max.(P, 1))` — the known-variance weight type, for
# which `vcov` returns inv(JᵀWJ) directly with no MSE rescaling, so `z0_err`
# carries the counting-statistics scale. The unweighted fit gets the
# MSE-scaled covariance, matching scipy's default (`absolute_sigma=False`).
# The scipy-vs-LsqFit agreement on both conventions is pinned by
# test/runtests.jl against the frozen py reference.
#
# Caveat on weighting: after MLEM the voxels in a slice are correlated, so
# P(z) is only approximately Poisson. The √counts weighting is a practical
# heuristic that improves the fit in the sparse tail; the authoritative
# uncertainty on the endpoint is the ensemble spread σ_R, NOT the
# per-realisation `z0_err` returned here.

"""
    erfc_edge(z, p)

Monotonic distal falloff: plateau `p[2]` above base `p[1]`, dropping through
its half-height at `p[3]` with width `p[4]`. `0.5*erfc(0) = 0.5`, so `p[3]`
is the R50 point of the fitted edge by construction.
"""
erfc_edge(z, p) = @. p[1] + p[2] * 0.5 * erfc((z - p[3]) / (sqrt(2.0) * p[4]))

"""
    fit_endpoint(z, profile; window, levels=(0.5, 0.8, 0.2),
                 weighted=true, p0=nothing) -> NamedTuple

Fit the distal falloff inside a FIXED window and read the R-levels from the
fit — the fit noise-averages the edge and yields a per-realisation
uncertainty, the key to σ_R in the sparse regime.

# Arguments
- `z`, `profile`: output of [`depth_profile`](@ref), proximal → distal in
  increasing `z`; the distal edge is the last falling edge.
- `window`: `(z_lo, z_hi)` fixed distal analysis window (from
  [`distal_window`](@ref)). It excludes the proximal rise and plateau so the
  single erfc edge is an appropriate model.
- `levels`: falloff fractions of the fitted plateau to report (0.5 primary).
- `weighted`: weight each point by its Poisson sigma `√max(counts, 1)` with
  the absolute-sigma convention, so `z0_err` reflects counting statistics.
- `p0`: initial `(base, amp, z0, w)`; auto-estimated when `nothing`.

# Returns
NamedTuple `(R, z0, w, z0_err, n_points, popt, pcov)` where `R` maps each
level to its depth (mm). Endpoints are `NaN` on fit failure or fewer than 4
window points — count and exclude them downstream; a high failure rate at a
dose level is itself a result.
"""
function fit_endpoint(z::AbstractVector{<:Real}, profile::AbstractVector{<:Real};
                      window, levels=(0.5, 0.8, 0.2),
                      weighted::Bool=true, p0=nothing)
    zlo, zhi = window
    sel = (z .>= zlo) .& (z .<= zhi)
    zf = Float64.(z[sel])
    pf = Float64.(profile[sel])
    n = length(zf)
    nanres = (R=Dict(Float64(lv) => NaN for lv in levels),
              z0=NaN, w=NaN, z0_err=NaN, n_points=n, popt=nothing, pcov=nothing)
    n < 4 && return nanres

    if p0 === nothing
        k = max(2, n ÷ 5)
        base0 = median(pf[end-k+1:end])       # distal tail
        plateau0 = median(pf[1:k])            # proximal plateau
        amp0 = max(plateau0 - base0, 1e-9)
        half = base0 + 0.5 * amp0
        above = findall(>=(half), pf)
        z0_0 = isempty(above) ? zf[n ÷ 2 + 1] : zf[last(above)]
        w0 = 0.1 * (zf[end] - zf[1]) + 1e-6
        p0 = (base0, amp0, z0_0, w0)
    end
    p0v = collect(Float64, p0)

    try
        fit = weighted ?
            curve_fit(erfc_edge, zf, pf, PrecisionWeights(1 ./ max.(pf, 1.0)),
                      p0v; maxIter=20_000) :
            curve_fit(erfc_edge, zf, pf, p0v; maxIter=20_000)
        fit.converged || return nanres
        base, amp, z0, w = coef(fit)
        pcov = vcov(fit)
        perr = sqrt.(max.(diag(pcov), 0.0))
        R = Dict(Float64(lv) => z0 + sqrt(2.0) * w * erfcinv(2.0 * lv)
                 for lv in levels)
        return (R=R, z0=z0, w=abs(w), z0_err=perr[3], n_points=n,
                popt=coef(fit), pcov=pcov)
    catch
        return nanres
    end
end

"""
    sigma_R(endpoints; dose_bragg_peak=nothing) -> NamedTuple

Collapse the per-realisation endpoints at ONE dose level into the two numbers
that matter: the mean (systematic activity-edge position) and the spread
(statistical precision σ_R, the geometry discriminator).

`endpoints` holds the R50 (or chosen level) from each realisation; `NaN`s
(failed fits) are dropped and counted. When `dose_bragg_peak` (dose-R80 from
the truth bundle, same axis and units) is given, the mean endpoint is also
returned as an offset from it — the activity-edge-to-range proxy distance,
largely common-mode across geometries.

# Returns
NamedTuple `(n_ok, n_fail, mean, sigma, sem, offset)`: mean endpoint, its
sample std (`== σ_R`, ddof = 1), the standard error on the mean, and
`offset = mean - dose_bragg_peak` (`nothing` when no peak is supplied).
Fewer than 2 finite endpoints yield `NaN` statistics.
"""
function sigma_R(endpoints; dose_bragg_peak=nothing)
    e = collect(Float64, endpoints)
    ok = filter(isfinite, e)
    n_ok, n_fail = length(ok), length(e) - length(ok)
    n_ok < 2 && return (n_ok=n_ok, n_fail=n_fail, mean=NaN, sigma=NaN,
                        sem=NaN, offset=nothing)
    m, s = mean(ok), std(ok)
    return (n_ok=n_ok, n_fail=n_fail, mean=m, sigma=s, sem=s / sqrt(n_ok),
            offset=dose_bragg_peak === nothing ? nothing : m - dose_bragg_peak)
end
