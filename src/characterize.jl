# characterize.jl — validation ladder rung 1: read the truth/ bundle and lock
# the reference before any reconstruction. The reference is three numbers on
# the shared z-frame: **dose-R80** (distal 80%-of-max of the dose),
# **true activity-R50** (half-height of the total-activity distal edge), and
# their **offset** — the ground-truth activity-to-range distance every
# reconstructed edge is scored against. The bundle contract is
# dev/upstream_response_truth_bundle.md; the four truth figures come from
# tools/plot_truth.py.

"""
    distal_crossing(z, y; level, reference=maximum(y)) -> Float64

Depth of the LAST downward crossing of `level·reference`, by linear
interpolation between the bracketing samples — the assumption-free reading
of a distal edge (dose-R80 at `level = 0.8`, activity-R50 at `level = 0.5`).
`z` runs proximal → distal in increasing order. Returns `NaN` when the curve
never crosses the level.
"""
function distal_crossing(z::AbstractVector{<:Real}, y::AbstractVector{<:Real};
                         level::Real, reference::Real=maximum(y))
    length(z) == length(y) || throw(ArgumentError("z and y lengths differ"))
    thr = level * reference
    for i in length(y)-1:-1:1
        if y[i] >= thr > y[i+1]
            return z[i] + (thr - y[i]) * (z[i+1] - z[i]) / (y[i+1] - y[i])
        end
    end
    return NaN
end

"""
    windowed_crossing(z, profile, window) -> Float64

Half-height crossing read INSIDE the fixed window, against the window's own
proximal-plateau and distal-tail medians (first/last fifth of the window
points) — the crossing-convention reading of a reconstructed or detected
profile, comparable to the truth activity-R50 crossing. Reading locally
matters: a global-max threshold misreads profiles tilted by attenuation or
sloped by reconstruction. Returns `NaN` on fewer than 4 window points.
"""
function windowed_crossing(z::AbstractVector{<:Real},
                           profile::AbstractVector{<:Real}, window)
    lo, hi = window
    sel = (z .>= lo) .& (z .<= hi)
    zf, pf = Float64.(z[sel]), Float64.(profile[sel])
    length(zf) < 4 && return NaN
    k = max(2, length(pf) ÷ 5)
    plateau, base = median(pf[1:k]), median(pf[end-k+1:end])
    return distal_crossing(zf, pf; level=1.0,
                           reference=base + 0.5 * (plateau - base))
end

"""
    read_depth_dose(truth_dir) -> NamedTuple

`truth/depth_dose.csv` as `(z_mm, dose_core_Gy, edep_total_MeV)` — the core
dose column carries the SOBP the range is defined on.
"""
function read_depth_dose(truth_dir::AbstractString)
    t = read_csv_table(joinpath(truth_dir, "depth_dose.csv"))
    return (z_mm=t["z_mm"], dose_core_Gy=t["dose_core_Gy"],
            edep_total_MeV=t["edep_total_MeV"])
end

"""
    read_activity_profile(truth_dir; budget="fast") -> NamedTuple

`truth/activity_profile_<budget>.csv` as
`(z_mm, isotopes::Vector{String}, per_isotope::Matrix, total)` — expected
decay counts at the leaf's dose, on the exact z-frame of `depth_dose.csv`,
escaped positrons excluded (the same source the LOR shards materialize
from). `per_isotope` columns follow `isotopes` order.
"""
function read_activity_profile(truth_dir::AbstractString; budget::AbstractString="fast")
    path = joinpath(truth_dir, "activity_profile_$(budget).csv")
    t = read_csv_table(path)
    # Isotope columns in file order (isotope-id order per the bundle contract).
    header = String.(strip.(split(strip(readline(path)), ',')))
    isotopes = filter(k -> k ∉ ("z_mm", "total"), header)
    per_isotope = reduce(hcat, (t[k] for k in isotopes))
    return (z_mm=t["z_mm"], isotopes=isotopes, per_isotope=per_isotope,
            total=t["total"])
end

"""
    TruthReference

The locked rung-1 reference (all depths mm, shared z-frame):

- `dose_R80` — distal 80%-of-max crossing of `dose_core_Gy`.
- `activity_R50` — distal half-height crossing of the total activity
  (interpolated, assumption-free: the locked value).
- `activity_R50_fit`, `activity_w_fit` — the same edge through
  [`fit_endpoint`](@ref) in `window` (the estimator's reading; a
  cross-check against the crossing, not the reference).
- `offset = activity_R50 - dose_R80` — the activity-edge-to-range distance;
  largely common-mode across geometries.
- `window` — the fixed distal fit window implied by the nominal edge.
- `scenario`, `budget`.
"""
struct TruthReference
    scenario::String
    budget::String
    dose_R80::Float64
    activity_R50::Float64
    activity_R50_fit::Float64
    activity_w_fit::Float64
    offset::Float64
    window::Tuple{Float64,Float64}
end

"""
    characterize(scenario_dir; budget="fast", proximal_margin_mm=20.0,
                 distal_margin_mm=15.0) -> TruthReference

Lock the starting-point reference from the `truth/` bundle: dose-R80,
true activity-R50 (interpolated crossing, with the windowed erfc fit as a
cross-check), their offset, and the fixed distal window built from the
interpolated activity edge. Runs before any reconstruction.
"""
function characterize(scenario_dir::AbstractString; budget::AbstractString="fast",
                      proximal_margin_mm::Real=20.0, distal_margin_mm::Real=15.0)
    tdir = truth_dir(scenario_dir)
    dose = read_depth_dose(tdir)
    act = read_activity_profile(tdir; budget=budget)
    act.z_mm == dose.z_mm ||
        error("characterize: activity and dose z-frames differ — broken bundle")

    dose_R80 = distal_crossing(dose.z_mm, dose.dose_core_Gy; level=0.8)
    activity_R50 = distal_crossing(act.z_mm, act.total; level=0.5)
    (isnan(dose_R80) || isnan(activity_R50)) &&
        error("characterize: no distal crossing found (dose_R80=$dose_R80, " *
              "activity_R50=$activity_R50)")

    window = distal_window(activity_R50; proximal_margin_mm=proximal_margin_mm,
                           distal_margin_mm=distal_margin_mm)
    fit = fit_endpoint(act.z_mm, act.total; window=window, weighted=true)

    return TruthReference(basename(abspath(scenario_dir)), String(budget),
                          dose_R80, activity_R50, fit.z0, fit.w,
                          activity_R50 - dose_R80, window)
end

"""
    write_reference(ref, path)

Write a `TruthReference` as TOML — the regenerable rung-1 artifact (e.g.
`out/characterize/truth_reference.toml`) later runs and the `config/` run parameters
read from.
"""
function write_reference(ref::TruthReference, path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, Dict(
            "scenario" => ref.scenario, "budget" => ref.budget,
            "dose_R80_mm" => ref.dose_R80,
            "activity_R50_mm" => ref.activity_R50,
            "activity_R50_fit_mm" => ref.activity_R50_fit,
            "activity_w_fit_mm" => ref.activity_w_fit,
            "offset_mm" => ref.offset,
            "window_mm" => collect(ref.window)))
    end
    return path
end

function Base.show(io::IO, ::MIME"text/plain", r::TruthReference)
    fmt(x) = string(round(x; digits=3))
    println(io, "TruthReference  $(r.scenario) / $(r.budget)")
    println(io, "  dose-R80        $(fmt(r.dose_R80)) mm")
    println(io, "  activity-R50    $(fmt(r.activity_R50)) mm (crossing; " *
                "erfc fit $(fmt(r.activity_R50_fit)) mm, w $(fmt(r.activity_w_fit)) mm)")
    println(io, "  offset          $(fmt(r.offset)) mm (activity-R50 − dose-R80)")
    print(io,   "  distal window   ($(fmt(r.window[1])), $(fmt(r.window[2]))) mm")
end
