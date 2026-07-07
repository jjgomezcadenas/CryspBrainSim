# output.jl — the output directory layout, mirroring the products-tree axes so
# results decouple by what they depend on (dev/PLAN.md → "Output layout"):
#
#   out/<scenario>/
#     truth/   mumap/                                      [scenario tier]
#     <topology>/<ring>/sensitivity/                        [scenario × ring tier]
#     <topology>/<ring>/<crystal>/                          [config tier]
#       shard_stats/  one_shard/  sigma_r/  origin_profile/
#     out/validation/                                       package cross-checks
#
# The depth of a result equals the scope it depends on: truth and the μ-map are
# scenario-general; the sensitivity image depends on the ring geometry and the
# object attenuation but not the crystal, so one base serves every crystal on a
# ring; everything else is per-configuration.

# Radiation length X0 (mm) per crystal material — the material-independent unit
# for detector thickness. Expressing a wall in X0 lets BGO and CsI be compared
# by stopping power rather than by millimetres.
const CRYSTAL_X0_MM = Dict("BGO" => 11.18, "CSI" => 18.6)

"""
    crystal_label(material, wall_mm) -> String

The crystal directory name, folding material and thickness together:
`crystal_label("BGO", 37.0) == "bgo_3X0"` (37 mm ÷ 11.18 mm/X0 ≈ 3X0).
"""
function crystal_label(material::AbstractString, wall_mm::Real)
    key = uppercase(material)
    haskey(CRYSTAL_X0_MM, key) ||
        error("crystal_label: no radiation length on file for material $material")
    return "$(lowercase(material))_$(round(Int, wall_mm / CRYSTAL_X0_MM[key]))X0"
end

"The out/ root (the repo's output directory)."
out_root() = joinpath(pkgdir(CryspBrainSim), "out")

"`out/<scenario>/` — the root of one scenario's results."
scenario_out(scenario; root=out_root()) = joinpath(root, scenario)

"`out/<scenario>/truth/` — the truth reference and its figures (scenario tier)."
truth_out(scenario; root=out_root()) = joinpath(scenario_out(scenario; root), "truth")

"`out/<scenario>/mumap/` — the voxel μ-map (scenario tier)."
mumap_out(scenario; root=out_root()) = joinpath(scenario_out(scenario; root), "mumap")

"`out/<scenario>/<topology>/<ring>/` — one ring geometry's results."
ring_out(scenario, topology, ring; root=out_root()) =
    joinpath(scenario_out(scenario; root), topology, ring)

"`out/<scenario>/<topology>/<ring>/sensitivity/` — the sensitivity base cache."
sensitivity_out(scenario, topology, ring; root=out_root()) =
    joinpath(ring_out(scenario, topology, ring; root), "sensitivity")

"`out/<scenario>/<topology>/<ring>/<crystal>/` — one PET configuration's results."
config_out(scenario, topology, ring, crystal; root=out_root()) =
    joinpath(scenario_out(scenario; root), topology, ring, crystal)

"`out/validation/` — package cross-checks that belong to no scenario."
validation_out(; root=out_root()) = joinpath(root, "validation")

# ---------------------------------------------------------------------------
# Self-describing scanner descriptors, so a reader of out/ can answer "what
# scanner produced this?" without reaching back to the products tree.
# ---------------------------------------------------------------------------

"""
    scanner_spec(geo) -> Dict

The full numeric ring description for a [`scanner_geometry`](@ref) result:
the raw geometry plus the derived quantities — outer radius, total length,
crystal count (`n_phi × n_z`), and the crystal block size (transverse from the
azimuthal pitch at the inner face, axial from the length per wheel, radial =
wall depth). Values in mm.
"""
function scanner_spec(geo)
    r_outer = geo.r_inner_mm + geo.wall_mm
    length_mm = 2 * geo.half_length_mm
    return Dict(
        "name" => geo.name, "shape" => geo.shape,
        "r_inner_mm" => geo.r_inner_mm, "wall_mm" => geo.wall_mm,
        "r_outer_mm" => r_outer,
        "half_length_mm" => geo.half_length_mm, "length_mm" => length_mm,
        "n_phi" => geo.n_phi, "n_z" => geo.n_z,
        "n_crystals" => geo.n_phi * geo.n_z,
        "crystal_transverse_mm" => 2π * geo.r_inner_mm / geo.n_phi,
        "crystal_axial_mm" => length_mm / geo.n_z,
        "crystal_radial_mm" => geo.wall_mm)
end

"""
    write_ring_geometry(geo; scenario, topology, ring, root=out_root()) -> path

Write `geometry.toml` (the full [`scanner_spec`](@ref)) at the ring tier
`out/<scenario>/<topology>/<ring>/`. Shared across the crystals on the ring.
"""
function write_ring_geometry(geo; scenario, topology, ring, root=out_root())
    dir = ring_out(scenario, topology, ring; root)
    mkpath(dir)
    path = joinpath(dir, "geometry.toml")
    open(io -> TOML.print(io, scanner_spec(geo)), path, "w")
    return path
end

"""
    write_crystal_spec(; scenario, topology, ring, crystal, material, wall_mm,
                       detector=nothing, root=out_root()) -> path

Write `crystal.toml` at the crystal tier: material, wall depth, radiation
length, and the exact thickness in X0 (the directory name carries only the
rounded value). When `detector` is given (a NamedTuple with fields
`energy_resolution_fwhm`, `sigma_xyz_mm`, `emin_keV`, `tau_ns` — from the shard
attributes), a `[detector]` section records the response. `sigma_xyz_mm` is a
single isotropic value: a monolithic crystal (as in CRYSP) resolves the three
coordinates equally.
"""
function write_crystal_spec(; scenario, topology, ring, crystal, material,
                            wall_mm, detector=nothing, root=out_root())
    key = uppercase(material)
    haskey(CRYSTAL_X0_MM, key) ||
        error("write_crystal_spec: no radiation length on file for $material")
    x0 = CRYSTAL_X0_MM[key]
    dir = config_out(scenario, topology, ring, crystal; root)
    mkpath(dir)
    d = Dict{String,Any}("material" => uppercase(material), "label" => crystal,
                         "wall_mm" => Float64(wall_mm), "X0_mm" => x0,
                         "thickness_X0" => wall_mm / x0)
    detector === nothing || (d["detector"] = Dict(
        "energy_resolution_fwhm" => Float64(detector.energy_resolution_fwhm),
        "sigma_xyz_mm" => Float64(detector.sigma_xyz_mm),
        "emin_keV" => Float64(detector.emin_keV),
        "tau_ns" => Float64(detector.tau_ns)))
    path = joinpath(dir, "crystal.toml")
    open(io -> TOML.print(io, d), path, "w")
    return path
end
