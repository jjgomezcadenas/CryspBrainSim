# products.jl — navigation and IO for the PtCryspProds tree (the contract in
# dev/reference/PRODUCTS.md, mirrored at PtCryspProds/README.md). This file
# owns the self-describing interface: glob shards under a config leaf, verify
# the provenance attrs on every read, pool shards into one event list, and
# read the shared per-scenario inputs — scanner_geometry.json, the phantom/
# files, and the truth/ bundle (delivered upstream; see
# dev/upstream_response_truth_bundle.md).

# ---------------------------------------------------------------------------
# Small CSV table reader — the products CSVs are machine-written (no quoting,
# no embedded commas), so a dependency-free reader covers them. Each column
# parses to Float64 when every entry parses, otherwise stays String.
# ---------------------------------------------------------------------------
"""
    read_csv_table(path) -> Dict{String,Vector}

Read a simple comma-separated table with a header row into a column dict.
Columns whose every entry parses as a number come back as `Vector{Float64}`,
the rest as `Vector{String}`. Covers the machine-written products CSVs
(truth/, phantom/); it is not a general CSV parser.
"""
function read_csv_table(path::AbstractString)
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && error("read_csv_table: $path is empty")
    header = strip.(split(lines[1], ','))
    cells = [strip.(split(l, ',')) for l in lines[2:end]]
    all(length.(cells) .== length(header)) ||
        error("read_csv_table: ragged rows in $path")
    table = Dict{String,Vector}()
    for (j, name) in enumerate(header)
        raw = [row[j] for row in cells]
        nums = tryparse.(Float64, raw)
        table[name] = any(isnothing, nums) ? raw : Vector{Float64}(nums)
    end
    return table
end

# ---------------------------------------------------------------------------
# Config-leaf navigation and shard IO
# ---------------------------------------------------------------------------
"""
The provenance attrs every `lors_shardNNN.h5` must carry (dev/PLAN.md,
W1): source identity (scenario, budget, dose, master_seed,
realization = the shard index), detector identity (crystal, n_phi, n_z,
windows, response), normalization (nevents, nrows), and the quantization
scales. `shard_attrs` fails loudly when any is absent.
"""
const REQUIRED_SHARD_ATTRS = [
    "scenario", "crystal", "budget", "dose_Gy", "master_seed", "realization",
    "n_phi", "n_z", "tau_ns", "emin_keV", "eres", "sigma_xyz_mm",
    "nevents", "nrows", "xyz_scale_mm", "e_scale_keV",
]

# Attrs that identify the master a shard belongs to: equal across the shards
# of one config leaf (the realization/shard index and per-shard row counts are
# the ones that differ).
const COMMON_MODE_ATTRS = [
    "scenario", "crystal", "budget", "dose_Gy", "master_seed",
    "n_phi", "n_z", "tau_ns", "emin_keV", "eres", "sigma_xyz_mm",
    "xyz_scale_mm", "e_scale_keV",
]

"""
    leaf_dir(root; scenario, scanner, crystal=nothing, leaf) -> String

Path of a config leaf in the products tree:
`root/scenario/scanner/[crystal/]leaf` (e.g. `leaf = "fast_1Gy"`). Pass
`crystal = nothing` for a heterogeneous scanner, where the leaf hangs
directly under the scanner. Errors when the directory is absent.
"""
function leaf_dir(root::AbstractString; scenario::AbstractString,
                  scanner::AbstractString,
                  crystal::Union{Nothing,AbstractString}=nothing,
                  leaf::AbstractString)
    parts = crystal === nothing ? (scenario, scanner, leaf) :
                                  (scenario, scanner, crystal, leaf)
    dir = joinpath(root, parts...)
    isdir(dir) || error("leaf_dir: $dir does not exist")
    return dir
end

"""
    shard_files(leaf) -> Vector{String}

The sorted `lors_shard*.h5` paths under a config leaf — sorted by name, which
the zero-padded naming makes shard-index order. Errors when the leaf holds no
shards.
"""
function shard_files(leaf::AbstractString)
    files = sort(filter(f -> occursin(r"^lors_shard\d+\.h5$", basename(f)),
                        readdir(leaf; join=true)))
    isempty(files) && error("shard_files: no lors_shard*.h5 under $leaf")
    return files
end

"""
    shard_attrs(file) -> Dict{String,Any}

The HDF5 root attributes of a shard, verified against
[`REQUIRED_SHARD_ATTRS`](@ref) — a missing provenance attr is a broken
contract and errors with the file path and the missing keys.
"""
function shard_attrs(file::AbstractString)
    attrs = h5open(file, "r") do f
        Dict{String,Any}(k => read_attribute(f, k) for k in keys(attributes(f)))
    end
    missing_keys = filter(k -> !haskey(attrs, k), REQUIRED_SHARD_ATTRS)
    isempty(missing_keys) ||
        error("shard_attrs: $file lacks provenance attrs: $(join(missing_keys, ", "))")
    return attrs
end

"""
    read_shard(file) -> (coinc, attrs, n_dropped)

One shard as an `MCCoincidences` (via RecoCryspTools' reader, which rescales
the Int16 storage and drops degenerate LORs) plus its verified provenance
attrs and the degenerate-drop count `n_dropped = nrows - length(coinc)`,
logged per the contract.
"""
function read_shard(file::AbstractString)
    attrs = shard_attrs(file)
    coinc = read_coincidences(file)
    n_dropped = Int(attrs["nrows"]) - length(coinc)
    n_dropped > 0 &&
        @info "read_shard: $(basename(file)) dropped $n_dropped degenerate LOR(s)"
    return (coinc=coinc, attrs=attrs, n_dropped=n_dropped)
end

"""
    pool_shards(files) -> (coinc, shard_attrs, n_dropped)

Pool a config leaf's shards into one event list — the master the thinning
draws from. Verifies each shard's provenance attrs, checks the common-mode
attrs ([`COMMON_MODE_ATTRS`](@ref)) agree across shards and the
`realization` indices are all distinct (matching shards of one master),
and returns the concatenated `MCCoincidences`, the per-shard attr dicts,
and the total degenerate-drop count.
"""
function pool_shards(files::AbstractVector{<:AbstractString})
    isempty(files) && error("pool_shards: no files")
    shards = [read_shard(f) for f in files]
    ref = shards[1].attrs
    for s in shards[2:end]
        for k in COMMON_MODE_ATTRS
            s.attrs[k] == ref[k] || error(
                "pool_shards: attr $k differs across shards " *
                "($(s.attrs[k]) vs $(ref[k])) — not one master")
        end
    end
    # Refuse to mix products generations: the v2 shards moved the `t_decay_s`
    # zero to irradiation end and re-centred the phantom, so a v2/legacy pool is
    # meaningless (generation2_plan.md §5). Legacy shards carry no `generation`.
    gens = [shard_generation(s.attrs) for s in shards]
    allequal(gens) ||
        error("pool_shards: mixed products generations $(unique(gens)) — refuse to mix")
    reals = [Int(s.attrs["realization"]) for s in shards]
    allunique(reals) ||
        error("pool_shards: duplicate realization indices $reals — not one master")
    cs = [s.coinc for s in shards]
    coinc = MCCoincidences(
        hcat((c.xstart for c in cs)...), hcat((c.xend for c in cs)...),
        hcat((c.origin for c in cs)...), vcat((c.truth for c in cs)...),
        vcat((c.nscat for c in cs)...), hcat((c.elem1 for c in cs)...),
        hcat((c.elem2 for c in cs)...), hcat((c.energy for c in cs)...))
    return (coinc=coinc, shard_attrs=[s.attrs for s in shards],
            n_dropped=sum(s.n_dropped for s in shards))
end

# ---------------------------------------------------------------------------
# Shared per-scenario inputs: scanner geometry, phantom, truth bundle
# ---------------------------------------------------------------------------
"""
    shard_generation(attrs) -> String

The products generation of a shard, `"v2"` for generation-2 shards and
`"legacy"` for the off-centre masters that carry no `generation` attr. The v2
guard: consumers must refuse to mix generations, since v2 moved the `t_decay_s`
zero to irradiation end and re-centred the phantom (generation2_plan.md §5).
"""
shard_generation(attrs::AbstractDict) = String(get(attrs, "generation", "legacy"))

"""
    shard_t_decay(file) -> Vector{Float32}

The absolute decay time of every coincidence (`t_decay_s`, for randoms gamma 1's
decay). The zero is the `t_decay_zero` attr — `acquisition_start` in legacy
shards, `irradiation_end` in v2 (where each leaf is already band-cut to its
scenario window `[t_del, t_del+t_ac]`). Read raw from the file, so the vector
aligns with `read_shard`'s coincidence list only when no degenerate LORs were
dropped — assert `n_dropped == 0` before masking with it.
"""
function shard_t_decay(file::AbstractString)
    return h5open(f -> read(f, "t_decay_s"), file, "r")
end

"""
    shard_isotope(file) -> Vector{Int8}

The emitting isotope id of every coincidence (`isotope` column, v2 only:
0=O15, 1=C11, 2=N13, 3=C10, 4=O14; for randoms gamma 1's decay). Read raw, so it
aligns with `read_shard`'s coincidence list only when `n_dropped == 0`. Enables
the exact per-species σ_R and the exact per-species washout keep (`washout_g`),
replacing the label-free posterior surrogate.
"""
function shard_isotope(file::AbstractString)
    return h5open(f -> read(f, "isotope"), file, "r")
end

"""
    scanner_geometry(scanner_dir) -> NamedTuple

The ring geometry from `<scanner>/scanner_geometry.json`, in mm:
`(name, shape, r_inner_mm, wall_mm, half_length_mm, n_phi, n_z)`. `n_phi` is
the crystal count per wheel (azimuthal), `n_z` the wheel count (axial).
Accepts either the scanner directory or the JSON path itself.
"""
function scanner_geometry(scanner_dir::AbstractString)
    path = isdir(scanner_dir) ? joinpath(scanner_dir, "scanner_geometry.json") :
                                scanner_dir
    isfile(path) || error("scanner_geometry: $path does not exist")
    g = JSON3.read(read(path, String))
    s = g["scanner"]
    return (name=String(s["name"]), shape=String(s["shape"]),
            r_inner_mm=10.0 * s["r_inner_cm"],
            wall_mm=10.0 * s["wall_thickness_cm"],
            half_length_mm=10.0 * s["half_length_cm"],
            n_phi=Int(s["n_phi"]), n_z=Int(s["n_z"]))
end

"""
    phantom_region(scenario_dir; region=1) -> NamedTuple

One region row of `<scenario>/phantom/phantom_regions.csv`, in mm:
`(region, material, solid, semi_axes, centre)`. The uniform head is the
single region; multi-region phantoms index by row.
"""
function phantom_region(scenario_dir::AbstractString; region::Integer=1)
    t = read_csv_table(joinpath(scenario_dir, "phantom", "phantom_regions.csv"))
    n = length(t["region"])
    1 <= region <= n || error("phantom_region: region $region of $n")
    return (region=t["region"][region], material=t["material"][region],
            solid=t["solid"][region],
            semi_axes=(t["a_mm"][region], t["b_mm"][region], t["c_mm"][region]),
            centre=(t["cx_mm"][region], t["cy_mm"][region], t["cz_mm"][region]))
end

"""
    material_mu(scenario_dir, material) -> NamedTuple

μ at 511 keV for a phantom material, from
`<scenario>/phantom/phantom_material_<key>_meta.csv`:
`(material, mu_mm_inv, density_g_cm3, energy_keV)`. `material` is the
G4 name (e.g. `"G4_BRAIN_ICRP"`); the file key is its lower-case form.
"""
function material_mu(scenario_dir::AbstractString, material::AbstractString)
    key = lowercase(material)
    path = joinpath(scenario_dir, "phantom", "phantom_material_$(key)_meta.csv")
    isfile(path) || error("material_mu: $path does not exist")
    t = read_csv_table(path)
    return (material=t["material"][1], mu_mm_inv=t["mu_mm_inv"][1],
            density_g_cm3=t["density_g_cm3"][1], energy_keV=t["energy_keV"][1])
end

"""
    truth_dir(scenario_dir) -> String

The `<scenario>/truth/` bundle directory (detector-independent, shared
across every scanner and crystal). Errors when the bundle is absent — the
bundle is delivered for the current master, so an absence is a broken tree,
not a fallback case.
"""
function truth_dir(scenario_dir::AbstractString)
    dir = joinpath(scenario_dir, "truth")
    isdir(dir) || error("truth_dir: $dir does not exist — the truth/ bundle " *
                        "is part of the products contract")
    return dir
end
