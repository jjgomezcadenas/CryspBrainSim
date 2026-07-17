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
    gaussian_smooth(profile, dz_mm, fwhm_mm) -> Vector{Float64}

1-D Gaussian smoothing of a uniformly-sampled profile (`dz_mm` spacing) with a
kernel of `fwhm_mm` full width at half maximum, edge-normalised (each output
divided by the in-range kernel weight, so the plateau and tail are not pulled
down at the ends). `fwhm_mm <= 0` returns the profile unchanged. This is the
Grogg et al. PET-resolution smoothing (7 mm FWHM) applied before the linear
distal fit.
"""
function gaussian_smooth(profile::AbstractVector{<:Real}, dz_mm::Real, fwhm_mm::Real)
    p = Float64.(profile)
    fwhm_mm <= 0 && return p
    σ = fwhm_mm / (2.0 * sqrt(2.0 * log(2.0)))
    r = max(1, ceil(Int, 3σ / dz_mm))
    k = [exp(-0.5 * (j * dz_mm / σ)^2) for j in -r:r]
    n = length(p)
    out = similar(p)
    for i in 1:n
        acc = 0.0; wsum = 0.0
        for (m, j) in enumerate(-r:r)
            ii = i + j
            1 <= ii <= n || continue
            acc += k[m] * p[ii]; wsum += k[m]
        end
        out[i] = acc / wsum
    end
    return out
end

"""
    fit_endpoint_grogg(z, profile; window, weighted=true,
                       extent_mm=25.0, min_points=8, smooth_fwhm_mm=0.0) -> NamedTuple

Distal endpoint by the Grogg et al. estimator (IEEE TNS 60 (2013) 3290): the
x-intercept of a linear fit to the distal falloff. The comparison estimator to
the erfc R50 of [`fit_endpoint`](@ref), run on the same profile and window.

The fit region starts at the last distal maximum inside the window and extends
`extent_mm` beyond it; the fit range keeps that start and ends where the
(weighted) residual sum of squares per degree of freedom is smallest, over
candidates of at least `min_points` points — the paper's range selection with
the RSS normalised per dof, so ranges of different lengths compare. Two
adaptations to a reconstructed whole-plane profile, both stated here because
the paper leaves them open: the last distal maximum is the deepest local
maximum at ≥ 50% of the window maximum (a plain "last local max" can be a
noise bump in the empty tail), and `weighted=true` applies the same Poisson
precision weights `1/max(P,1)` as the erfc fit, so the two estimators differ
only in model, not in weighting convention (`weighted=false` is the paper's
literal unweighted least squares, with the MSE-scaled covariance).

`smooth_fwhm_mm > 0` Gaussian-smooths the profile first
([`gaussian_smooth`](@ref)) — the paper's PET-resolution smoothing (7 mm),
which stabilises the distal-maximum start selection at the cost of resolution.

# Returns
NamedTuple `(x0, x0_err, slope, n_points, z_first, z_last, rss_dof)`: the
x-intercept `-a/b` (mm), its propagated error, the fitted slope, the chosen
range (point count and z limits), and the selection score. `x0` is NaN when
the window holds fewer than `min_points` points past the start, the profile
has no positive maximum, or the best-range slope is not negative.
"""
function fit_endpoint_grogg(z::AbstractVector{<:Real}, profile::AbstractVector{<:Real};
                            window, weighted::Bool=true,
                            extent_mm::Real=25.0, min_points::Int=8,
                            smooth_fwhm_mm::Real=0.0)
    zlo, zhi = window
    prof = smooth_fwhm_mm > 0 ?
        gaussian_smooth(profile, length(z) > 1 ? Float64(z[2] - z[1]) : 1.0,
                        smooth_fwhm_mm) : profile
    sel = (z .>= zlo) .& (z .<= zhi)
    zf = Float64.(z[sel])
    pf = Float64.(prof[sel])
    n = length(zf)
    nanres = (x0=NaN, x0_err=NaN, slope=NaN, n_points=0,
              z_first=NaN, z_last=NaN, rss_dof=NaN)
    pmax = n == 0 ? 0.0 : maximum(pf)
    (n < min_points || pmax <= 0.0) && return nanres

    # Fit-region start: the last local maximum at ≥ 50% of the window maximum.
    islocmax(i) = (i == 1 || pf[i] >= pf[i-1]) && (i == n || pf[i] >= pf[i+1])
    starts = [i for i in 1:n if pf[i] >= 0.5 * pmax && islocmax(i)]
    isempty(starts) && return nanres
    i0 = last(starts)
    jmax = searchsortedlast(zf, zf[i0] + extent_mm)
    jmax - i0 + 1 < min_points && return nanres

    # Fixed start, variable end: weighted linear regression on each candidate,
    # keep the smallest RSS/dof.
    best = nothing
    for j in (i0 + min_points - 1):jmax
        zs = view(zf, i0:j); ps = view(pf, i0:j)
        w = weighted ? 1.0 ./ max.(ps, 1.0) : ones(length(ps))
        Sw = sum(w); Sx = sum(w .* zs); Sy = sum(w .* ps)
        Sxx = sum(w .* zs .^ 2); Sxy = sum(w .* zs .* ps)
        Δ = Sw * Sxx - Sx^2
        Δ <= 0 && continue
        b = (Sw * Sxy - Sx * Sy) / Δ
        a = (Sxx * Sy - Sx * Sxy) / Δ
        dof = length(ps) - 2
        rss_dof = sum(w .* (ps .- a .- b .* zs) .^ 2) / dof
        if best === nothing || rss_dof < best.rss_dof
            # Covariance: precision-weighted (absolute sigma) as-is; unweighted
            # MSE-scaled — the same two conventions as fit_endpoint.
            s2 = weighted ? 1.0 : rss_dof
            best = (a=a, b=b, rss_dof=rss_dof, j=j,
                    var_a=s2 * Sxx / Δ, var_b=s2 * Sw / Δ, cov_ab=-s2 * Sx / Δ)
        end
    end
    (best === nothing || best.b >= 0.0) && return nanres

    a, b = best.a, best.b
    x0 = -a / b
    var_x0 = best.var_a / b^2 + a^2 * best.var_b / b^4 - 2a * best.cov_ab / b^3
    return (x0=x0, x0_err=sqrt(max(var_x0, 0.0)), slope=b,
            n_points=best.j - i0 + 1, z_first=zf[i0], z_last=zf[best.j],
            rss_dof=best.rss_dof)
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
