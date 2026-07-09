# config.jl — the frozen common-mode run parameters (config/run_parameters.toml). Drivers and
# tools consume them unchanged; the file is the single source the sweep's
# common-mode discipline hangs on.

"""
    load_run_parameters(path=joinpath(pkgdir(CryspBrainSim), "config", "run_parameters.toml"))
        -> NamedTuple

The frozen run parameters as typed values: `grid = (n, img_origin, voxsize)` (tuples,
Float32 where the projectors want them), `roi = (radius_mm, centre_mm)`,
`window = (z_lo, z_hi)`, `niter`, `n_sens`, `chunk`, `truth_selection`, and
`config = (scenario, topology, scanner, crystal, leaf)` — the active arm of the
products tree that every driver and tool reads.
"""
function load_run_parameters(path::AbstractString=joinpath(pkgdir(CryspBrainSim),
                                                  "config", "run_parameters.toml"))
    k = TOML.parsefile(path)
    g = k["grid"]
    c = k["configuration"]
    return (config=(scenario=String(c["scenario"]), topology=String(c["topology"]),
                    scanner=String(c["scanner"]), crystal=String(c["crystal"]),
                    leaf=String(c["leaf"])),
            grid=(n=Tuple(Int.(g["n"])),
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
