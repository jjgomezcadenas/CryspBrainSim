# profile.jl — steps 4–5a of the range-verification pipeline: reduce a
# beam-aligned reconstructed image to a 1-D depth-activity profile and build
# the fixed distal analysis window. Julia port of py/depth_profile.py
# (depth_profile, distal_window); the frozen Python module is the
# cross-validation reference in test/runtests.jl.
#
# Design invariants, identical across the scanner arms:
#   * beam-aligned grid: the beam/depth axis is a grid axis, so profiling is a
#     plain transverse sum with no oblique resampling;
#   * the transverse ROI is a fixed disc centred on the BEAM AXIS — a
#     data-driven centre would inject a noisy shift into the endpoint and
#     inflate σ_R;
#   * integrate (sum) over the ROI, preserving the count statistics the
#     endpoint fit and σ_R ride on;
#   * the fit window is fixed from the NOMINAL range, common-mode across arms
#     and realisations.

"""
    depth_profile(image; voxel_size_mm, beam_axis=3,
                  roi_radius_mm=nothing, roi_centre_mm=(0.0, 0.0),
                  z_origin_mm=0.0) -> (z, prof)

Sum a beam-aligned reconstructed image over a fixed transverse disc ROI at
each depth slice.

# Arguments
- `image`: 3-D reconstructed activity on a beam-aligned grid, so that
  `beam_axis` is the depth direction.
- `voxel_size_mm`: voxel edge lengths (mm) along the three image axes.
- `beam_axis`: which image axis (1, 2 or 3) is the beam/depth axis.
- `roi_radius_mm`: radius of the transverse disc ROI (mm). `nothing` sums the
  whole transverse plane — keep that for quick looks only; it dilutes the
  distal contrast with noise. Size it generously,
  `roi_radius_mm ≳ n·√(σ_spot² + σ_PSF²) + R_β⁺` with n ≈ 3, so the
  limited-angle elongation of an open arm is measured rather than clipped.
- `roi_centre_mm`: transverse ROI centre relative to the grid centre, in the
  order of the two non-beam axes. Keep at (0, 0) to sit on the beam axis.
- `z_origin_mm`: depth coordinate of the first slice (voxel-centre convention).

# Returns
- `z::Vector{Float64}`: depth coordinate of each slice (mm).
- `prof::Vector{Float64}`: activity summed within the ROI at each depth.

Use the same `roi_radius_mm`, `roi_centre_mm` and grid for every scanner
geometry — the ROI is one of the frozen common-mode knobs (`config/`).
"""
function depth_profile(image::AbstractArray{<:Real,3};
                       voxel_size_mm,
                       beam_axis::Integer=3,
                       roi_radius_mm::Union{Nothing,Real}=nothing,
                       roi_centre_mm=(0.0, 0.0),
                       z_origin_mm::Real=0.0)
    length(voxel_size_mm) == 3 ||
        throw(ArgumentError("voxel_size_mm must have 3 entries"))
    beam_axis in (1, 2, 3) ||
        throw(ArgumentError("beam_axis must be 1, 2 or 3"))

    # Bring the beam axis to the front: img has shape (nz, nu, nv).
    other = filter(!=(beam_axis), (1, 2, 3))
    img = PermutedDimsArray(image, (beam_axis, other...))
    nz, nu, nv = size(img)

    dz = Float64(voxel_size_mm[beam_axis])
    du = Float64(voxel_size_mm[other[1]])
    dv = Float64(voxel_size_mm[other[2]])

    # Transverse voxel-centre coordinates relative to the grid centre.
    u = ((0:nu-1) .- (nu - 1) / 2) .* du
    v = ((0:nv-1) .- (nv - 1) / 2) .* dv
    cu, cv = Float64.(roi_centre_mm)
    mask = roi_radius_mm === nothing ? trues(nu, nv) :
        [(ui - cu)^2 + (vi - cv)^2 <= Float64(roi_radius_mm)^2
         for ui in u, vi in v]

    prof = Vector{Float64}(undef, nz)
    for k in 1:nz
        acc = 0.0
        @inbounds for j in 1:nv, i in 1:nu
            mask[i, j] && (acc += Float64(img[k, i, j]))
        end
        prof[k] = acc
    end
    z = z_origin_mm .+ (0:nz-1) .* dz
    return collect(Float64, z), prof
end

"""
    distal_window(z_edge_nominal_mm; proximal_margin_mm=20.0,
                  distal_margin_mm=15.0) -> (z_lo, z_hi)

Build the FIXED fit window bracketing the distal falloff, defined from the
nominal (expected) activity edge and identical across arms and realisations —
a per-realisation peak search would correlate the window with the noise being
measured. The window starts on the plateau proximal to the edge and ends in
the tail distal to it.
"""
distal_window(z_edge_nominal_mm::Real;
              proximal_margin_mm::Real=20.0, distal_margin_mm::Real=15.0) =
    (z_edge_nominal_mm - proximal_margin_mm,
     z_edge_nominal_mm + distal_margin_mm)
