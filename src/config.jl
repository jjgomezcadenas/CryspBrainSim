# config.jl — the frozen common-mode knobs (config/knobs.toml). Drivers and
# tools consume them unchanged; the file is the single source the sweep's
# common-mode discipline hangs on.

"""
    load_knobs(path=joinpath(pkgdir(CryspBrainSim), "config", "knobs.toml"))
        -> NamedTuple

The frozen knobs as typed values: `grid = (n, img_origin, voxsize)` (tuples,
Float32 where the projectors want them), `roi = (radius_mm, centre_mm)`,
`window = (z_lo, z_hi)`, `niter`, `n_sens`, `chunk`, `truth_selection`.
"""
function load_knobs(path::AbstractString=joinpath(pkgdir(CryspBrainSim),
                                                  "config", "knobs.toml"))
    k = TOML.parsefile(path)
    g = k["grid"]
    return (grid=(n=Tuple(Int.(g["n"])),
                  img_origin=Tuple(Float32.(g["img_origin_mm"])),
                  voxsize=Tuple(Float32.(g["voxsize_mm"]))),
            roi=(radius_mm=Float64(k["roi"]["radius_mm"]),
                 centre_mm=Tuple(Float64.(k["roi"]["centre_mm"]))),
            window=(Float64(k["window"]["z_lo_mm"]),
                    Float64(k["window"]["z_hi_mm"])),
            niter=Int(k["mlem"]["niter"]),
            n_sens=Int(k["sensitivity"]["n_sens"]),
            chunk=Int(k["sensitivity"]["chunk"]),
            truth_selection=String(k["selection"]["truth"]))
end
