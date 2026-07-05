# mumap.jl — attenuation from the scenario phantom files (dev/PLAN.md
# write-item 4). Two routes to the per-LOR survival factor a = exp(−μ·chord):
#
#   * analytic — `ellipsoid_chord` with the scenario's semi-axes, centre and
#     μ(511). Exact for the uniform head; the headline route.
#   * voxel μ-map — `exp.(-joseph3d_fwd(...))` over a voxelized μ. Exact for
#     any voxelized μ; the path the multi-region head phantoms take later.
#
# In PET a coincidence's attenuation depends only on the TOTAL chord its LOR
# cuts through the medium, never on where along the line the decay sat — so
# both routes weight LORs, not emission points. Everything is mm and Float32
# (the RecoCrysp conventions).

"""
    phantom_attenuation(scenario_dir; region=1) -> NamedTuple

The uniform-body attenuation parameters read from the scenario files:
`(semi_axes, centre, mu_mm_inv, material, solid)` — `phantom_regions.csv`
supplies the ellipsoid, `phantom_material_*_meta.csv` supplies μ at 511 keV.
Errors on a non-ellipsoid region: route those through the voxel μ-map.
"""
function phantom_attenuation(scenario_dir::AbstractString; region::Integer=1)
    reg = phantom_region(scenario_dir; region=region)
    reg.solid == "ellipsoid" ||
        error("phantom_attenuation: region solid is $(reg.solid), not ellipsoid — " *
              "build a voxel μ-map instead")
    mu = material_mu(scenario_dir, reg.material)
    return (semi_axes=Float32.(reg.semi_axes), centre=Float32.(reg.centre),
            mu_mm_inv=Float32(mu.mu_mm_inv), material=reg.material,
            solid=reg.solid)
end

"""
    attenuation_ellipsoid(xs, xe; semi_axes, centre, mu_mm_inv) -> Vector{Float32}

Per-LOR survival factors `exp(−μ·chord)` for the uniform ellipsoid — the
analytic route, exact for the uniform head. `xs`/`xe` are `(3, N)` endpoint
matrices (mm, CPU); the loop threads across LORs.
"""
function attenuation_ellipsoid(xs::AbstractMatrix, xe::AbstractMatrix;
                               semi_axes, centre, mu_mm_inv::Real)
    n = size(xs, 2)
    size(xe, 2) == n || throw(ArgumentError("xs and xe column counts differ"))
    a = Vector{Float32}(undef, n)
    μ = Float32(mu_mm_inv)
    Threads.@threads for i in 1:n
        chord = ellipsoid_chord((xs[1, i], xs[2, i], xs[3, i]),
                                (xe[1, i], xe[2, i], xe[3, i]);
                                axes=semi_axes, center=centre)
        a[i] = exp(-μ * chord)
    end
    return a
end

"""
    centered_grid(n, voxsize) -> NTuple{3,Float32}

The `img_origin` (world coordinate of the FIRST voxel CENTRE) that centres an
`n = (n1, n2, n3)` grid of `voxsize` voxels on the world origin — the
RecoCrysp image convention.
"""
centered_grid(n::NTuple{3,<:Integer}, voxsize) =
    ntuple(i -> Float32(-(n[i] - 1) / 2 * voxsize[i]), 3)

"""
    build_mumap(; n, img_origin, voxsize, semi_axes, centre, mu_mm_inv)
        -> Array{Float32,3}

Voxelized μ-map of the uniform ellipsoid: voxels whose CENTRE lies inside get
`mu_mm_inv`, the rest 0. Size the grid to cover the whole body (the head runs
to y = centre_y − b ≈ −117 mm) — the μ-map grid is independent of the activity
grid and each is sized to its job. The voxel route for this uniform body is
the cross-check of the analytic one; it becomes the primary route with the
multi-region phantoms.
"""
function build_mumap(; n::NTuple{3,<:Integer}, img_origin, voxsize,
                     semi_axes, centre, mu_mm_inv::Real)
    mumap = zeros(Float32, n)
    a2, b2, c2 = Float32.(semi_axes) .^ 2
    cx, cy, cz = Float32.(centre)
    μ = Float32(mu_mm_inv)
    for k in 1:n[3], j in 1:n[2], i in 1:n[1]
        x = img_origin[1] + (i - 1) * voxsize[1] - cx
        y = img_origin[2] + (j - 1) * voxsize[2] - cy
        z = img_origin[3] + (k - 1) * voxsize[3] - cz
        (x^2 / a2 + y^2 / b2 + z^2 / c2 <= 1.0f0) && (mumap[i, j, k] = μ)
    end
    return mumap
end

"""
    attenuation_mumap(xs, xe, mumap, img_origin, voxsize) -> Vector{Float32}

Per-LOR survival factors `exp(−∫μ dl)` through a voxelized μ-map via the
Joseph forward projector — exact for any voxelized μ; the multi-region route.
Arrays follow the projector's backend (CPU or GPU alike).
"""
attenuation_mumap(xs, xe, mumap, img_origin, voxsize) =
    exp.(.-joseph3d_fwd(xs, xe, mumap, img_origin, voxsize))
