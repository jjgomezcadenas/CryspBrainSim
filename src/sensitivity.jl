# sensitivity.jl — the MLEM normalization deliverable (dev/PLAN.md
# W5): the UNSCALED base = Aᵀ(a_geom) accumulated over chunked
# draws of geometric surface LORs from `ContinuousPET`, cached to NPZ with a
# provenance sidecar. `base` is realization-independent and reused across the
# whole sweep; the per-realization factor `n_events / n_sens` rides on each
# event list at reconstruction time (`scaled_sensitivity`).
#
# The project run parameter is n_sens = 10⁹ (dev/PLAN.md run-parameters table; the MC
# studies certified 5×10⁸ at 2.5 mm voxels, and our finer 1.5 mm grid holds
# ~4.6× fewer LOR crossings per voxel, so the sample is doubled — measured
# two-seed mottle 1.28% per image at the provisional corridor grid, ~37 s to
# build).
# Revalidate at the frozen grid before the sweep by comparing two seeds (the
# stage check in dev/PLAN.md).

"""
    sensitivity_base(; r_inner_mm, half_length_mm, n, img_origin, voxsize,
                     attenuation, n_sens, chunk=20_000_000, seed=1,
                     device=identity, progress=true) -> Array{Float32,3}

Accumulate the unscaled sensitivity `base = Aᵀ(a_geom)` for a closed ring:
draw `n_sens` uniform surface LORs from `ContinuousPET(diameter=2·r_inner,
afov=2·half_length)` in chunks, weight each chunk with
`attenuation(xs, xe) -> Vector{Float32}` (the [`attenuation_ellipsoid`](@ref)
or [`attenuation_mumap`](@ref) route), and backproject onto the
`(n, img_origin, voxsize)` grid.

LOR generation and attenuation run on CPU arrays; pass `device = MtlArray`
(with Metal loaded) to run each chunk's backprojection on the GPU. The RNG is
`MersenneTwister(seed)`, so a (seed, n_sens, chunk) triple reproduces the
draw; the backprojection accumulates atomically, so images reproduce to a
tolerance, never bit-for-bit.
"""
function sensitivity_base(; r_inner_mm::Real, half_length_mm::Real,
                          n::NTuple{3,<:Integer}, img_origin, voxsize,
                          attenuation, n_sens::Integer,
                          chunk::Integer=20_000_000, seed::Integer=1,
                          device=identity, progress::Bool=true)
    sc = ContinuousPET(diameter=2 * r_inner_mm, afov=2 * half_length_mm)
    rng = MersenneTwister(seed)
    base = device(zeros(Float32, n))
    done = 0
    while done < n_sens
        nb = Int(min(chunk, n_sens - done))
        gxs, gxe = sample_lors(sc, nb; rng=rng)
        ga = attenuation(gxs, gxe)
        base .+= sensitivity_image(device(gxs), device(gxe), n,
                                   img_origin, voxsize; weights=device(ga))
        done += nb
        progress && @info "sensitivity_base: $done / $n_sens LORs"
    end
    return Array(base)
end

"""
    sensitivity_cache_name(scanner_name, params; n_sens=params.n_sens) -> String

The cache basename (no directory, no extension) for a sensitivity base built on
`params`' frozen grid for `scanner_name`, e.g.
`crysp_ring_1m_grid64x64x96_1.5mm_orgm47.25_m47.25_m119.25_n1000000000`. The
grid origin is stamped in (`m` for a minus sign) so bases on the same-shape grid
at different origins never share a file. One rule, used by the builder
(`tools/make_sensitivity.jl`) and the drivers that load the cache.
"""
function sensitivity_cache_name(scanner_name::AbstractString, params;
                                n_sens::Integer=params.n_sens)
    n, vs, org = params.grid.n, params.grid.voxsize, params.grid.img_origin
    org_tag = join(replace.(string.(round.(Float64.(org); digits=2)), "-" => "m"), "_")
    return "$(scanner_name)_grid$(n[1])x$(n[2])x$(n[3])_$(vs[1])mm_org$(org_tag)_n$(n_sens)"
end

"""
    scaled_sensitivity(base, n_events, n_sens) -> Array{Float32,3}

The per-realization sensitivity `sens = base · (n_events / n_sens)` — the
scale convention of `sensitivity_image` for an independent geometric sample
(`RecoCryspUse.md` §5). Apply at reconstruction time; never bake it into the
cached base.
"""
scaled_sensitivity(base::AbstractArray{Float32,3}, n_events::Integer,
                   n_sens::Integer) = base .* Float32(n_events / n_sens)

"""
    recocrysp_sha() -> String

The git HEAD of the RecoCrysp checkout the session runs against (located via
`pathof(RecoCryspTools)`), for provenance records; `"unknown"` when git can't
answer. Every cached artifact carries the engine SHA that produced it.
"""
function recocrysp_sha()
    try
        root = dirname(dirname(dirname(pathof(RecoCryspTools))))
        return readchomp(`git -C $root rev-parse HEAD`)
    catch
        return "unknown"
    end
end

"""
    save_sensitivity(path, base; r_inner_mm, half_length_mm, n_sens, chunk,
                     seed, img_origin, voxsize, attenuation_meta) -> path

Cache the unscaled `base` to `<path>.npz` with the provenance sidecar
`<path>.toml`: the RecoCrysp SHA, the scanner and grid parameters, the draw
(n_sens, chunk, seed) and the attenuation description (a `Dict`/`NamedTuple`
of the μ route's parameters). `path` carries no extension.
"""
function save_sensitivity(path::AbstractString, base::Array{Float32,3};
                          r_inner_mm::Real, half_length_mm::Real,
                          n_sens::Integer, chunk::Integer, seed::Integer,
                          img_origin, voxsize, attenuation_meta)
    mkpath(dirname(abspath(path)))
    npzwrite(path * ".npz", Dict("base" => base))
    open(path * ".toml", "w") do io
        TOML.print(io, Dict(
            "recocrysp_sha" => recocrysp_sha(),
            "scanner" => Dict("r_inner_mm" => Float64(r_inner_mm),
                              "half_length_mm" => Float64(half_length_mm)),
            "draw" => Dict("n_sens" => n_sens, "chunk" => chunk, "seed" => seed),
            "grid" => Dict("n" => collect(size(base)),
                           "img_origin" => Float64.(collect(img_origin)),
                           "voxsize" => Float64.(collect(voxsize))),
            "attenuation" => Dict(String(k) => v isa Tuple ? Float64.(collect(v)) : v
                                  for (k, v) in pairs(attenuation_meta))))
    end
    return path
end

"""
    load_sensitivity(path) -> (base, meta)

Read a cached sensitivity base and its provenance sidecar (`<path>.npz` +
`<path>.toml`). Warns when the cache's RecoCrysp SHA differs from the
checkout in use — results then carry a different engine than the session.
"""
function load_sensitivity(path::AbstractString)
    base = npzread(path * ".npz")["base"]
    meta = TOML.parsefile(path * ".toml")
    sha = get(meta, "recocrysp_sha", "unknown")
    cur = recocrysp_sha()
    sha == cur ||
        @warn "load_sensitivity: cache built at RecoCrysp $sha, session runs $cur"
    return base, meta
end
