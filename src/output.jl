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

"`out/<scenario>/<topology>/<ring>/sensitivity/` — the sensitivity base cache."
sensitivity_out(scenario, topology, ring; root=out_root()) =
    joinpath(scenario_out(scenario; root), topology, ring, "sensitivity")

"`out/<scenario>/<topology>/<ring>/<crystal>/` — one PET configuration's results."
config_out(scenario, topology, ring, crystal; root=out_root()) =
    joinpath(scenario_out(scenario; root), topology, ring, crystal)

"`out/validation/` — package cross-checks that belong to no scenario."
validation_out(; root=out_root()) = joinpath(root, "validation")
