# qa.jl — validation ladder rung 2: a stored shard's statistics and sanity
# numbers as a struct the drivers assert on, turning "the shard looks fine"
# into an automatable gate. Runs on a shard at full statistics (distinct from
# a thinned realization); run it first on every new shard. The companion 3×3
# figure is tools/plot_shard.py.
#
# The QA reads the file's columns directly over HDF5 (not through
# read_coincidences): the gate scores the file as written — including the
# degenerate LORs the reconstruction reader drops — and needs the timing
# columns (`t1_ns`, `t2_ns`, `dt_ns`) that `MCCoincidences` does not carry.

"""
    ShardQA

One stored shard's statistics + sanity numbers (all positions mm, energies
keV). Fields:

- identity: `file`, `scenario`, `crystal`, `budget`, `shard` (the
  realization index).
- counting: `nrows`, `nevents`, `acceptance = nrows/nevents`, truth
  composition `n_true`/`n_scatter_single`/`n_scatter_multiple`/`n_random`
  with fractions `frac_true`/`frac_scatter`/`frac_random`, and
  `n_degenerate` (identical quantized endpoints — the LORs the
  reconstruction reader drops).
- per-hit sanity (both gammas pooled): `e_range` (min, max) keV,
  `r_range` (min, max) transverse hit radius, `doi_range = r − r_inner`.
- timing: `dt_frac_in_tau` — fraction of raw Δt = t1 − t2 within the
  coincidence window ±τ (1.0 by construction for an in-window list),
  `tau_ns`.
- source: `src_min`/`src_max` — the annihilation-point bounding box.
"""
struct ShardQA
    file::String
    scenario::String
    crystal::String
    budget::String
    shard::Int
    nrows::Int
    nevents::Int
    acceptance::Float64
    n_true::Int
    n_scatter_single::Int
    n_scatter_multiple::Int
    n_random::Int
    frac_true::Float64
    frac_scatter::Float64
    frac_random::Float64
    n_degenerate::Int
    e_range::Tuple{Float64,Float64}
    r_range::Tuple{Float64,Float64}
    doi_range::Tuple{Float64,Float64}
    dt_frac_in_tau::Float64
    tau_ns::Float64
    src_min::NTuple{3,Float64}
    src_max::NTuple{3,Float64}
end

"""
    shard_qa(file; r_inner_mm) -> ShardQA

Compute the QA statistics of one stored shard. `r_inner_mm` comes from
[`scanner_geometry`](@ref) and anchors the DOI range; the provenance attrs
are verified on read ([`shard_attrs`](@ref)).
"""
function shard_qa(file::AbstractString; r_inner_mm::Real)
    attrs = shard_attrs(file)
    xs = Float64(attrs["xyz_scale_mm"])
    es = Float64(attrs["e_scale_keV"])

    h5open(file, "r") do f
        col(n) = read(f[n])
        truth = col("truth")
        nrows = length(truth)
        Int(attrs["nrows"]) == nrows ||
            error("shard_qa: attr nrows=$(attrs["nrows"]) ≠ $(nrows) rows in $file")

        n_true = count(==(Int8(0)), truth)
        n_random = count(==(Int8(2)), truth)
        nscat = Int16.(col("nscat1")) .+ Int16.(col("nscat2"))
        scat = truth .== Int8(1)
        n_single = count(scat .& (nscat .== 1))
        n_multiple = count(scat .& (nscat .>= 2))

        # Quantized endpoints: degenerate = identical Int16 triples (the same
        # criterion RecoCryspTools drops at the file boundary, at quantized
        # precision).
        x1, y1, z1 = col("x1_mm"), col("y1_mm"), col("z1_mm")
        x2, y2, z2 = col("x2_mm"), col("y2_mm"), col("z2_mm")
        n_degenerate = count(i -> x1[i] == x2[i] && y1[i] == y2[i] && z1[i] == z2[i],
                             eachindex(x1))

        e1, e2 = Float64.(col("e1_keV")) .* es, Float64.(col("e2_keV")) .* es
        e_range = (min(minimum(e1), minimum(e2)), max(maximum(e1), maximum(e2)))

        r1 = hypot.(Float64.(x1) .* xs, Float64.(y1) .* xs)
        r2 = hypot.(Float64.(x2) .* xs, Float64.(y2) .* xs)
        r_range = (min(minimum(r1), minimum(r2)), max(maximum(r1), maximum(r2)))
        doi_range = r_range .- Float64(r_inner_mm)

        tau = Float64(attrs["tau_ns"])
        dtraw = Float64.(col("t1_ns")) .- Float64.(col("t2_ns"))
        dt_frac_in_tau = count(x -> abs(x) <= tau, dtraw) / nrows

        x0 = Float64.(col("x0_mm")) .* xs
        y0 = Float64.(col("y0_mm")) .* xs
        z0 = Float64.(col("z0_mm")) .* xs
        src_min = (minimum(x0), minimum(y0), minimum(z0))
        src_max = (maximum(x0), maximum(y0), maximum(z0))

        nevents = Int(attrs["nevents"])
        return ShardQA(
            String(file), String(attrs["scenario"]), String(attrs["crystal"]),
            String(attrs["budget"]), Int(attrs["realization"]),
            nrows, nevents, nrows / nevents,
            n_true, n_single, n_multiple, n_random,
            n_true / nrows, (n_single + n_multiple) / nrows, n_random / nrows,
            n_degenerate, e_range, r_range, doi_range,
            dt_frac_in_tau, tau, src_min, src_max)
    end
end

function Base.show(io::IO, ::MIME"text/plain", q::ShardQA)
    fmt(x) = string(round(x; sigdigits=4))
    println(io, "ShardQA  $(basename(q.file))")
    println(io, "  identity    $(q.scenario) / $(q.crystal) / $(q.budget) / shard $(q.shard)")
    println(io, "  counting    $(q.nrows) LORs from $(q.nevents) decays — acceptance $(fmt(100q.acceptance))%")
    println(io, "  composition true $(fmt(100q.frac_true))% · scatter $(fmt(100q.frac_scatter))% " *
                "(S $(q.n_scatter_single) / M $(q.n_scatter_multiple)) · random $(fmt(100q.frac_random))%")
    println(io, "  degenerate  $(q.n_degenerate) (dropped by the reconstruction reader)")
    println(io, "  energy      [$(fmt(q.e_range[1])), $(fmt(q.e_range[2]))] keV")
    println(io, "  hit radius  [$(fmt(q.r_range[1])), $(fmt(q.r_range[2]))] mm — DOI " *
                "[$(fmt(q.doi_range[1])), $(fmt(q.doi_range[2]))] mm")
    println(io, "  Δt window   $(fmt(100q.dt_frac_in_tau))% within ±$(q.tau_ns) ns")
    print(io,   "  source box  x $(fmt(q.src_min[1]))…$(fmt(q.src_max[1]))  " *
                "y $(fmt(q.src_min[2]))…$(fmt(q.src_max[2]))  " *
                "z $(fmt(q.src_min[3]))…$(fmt(q.src_max[3])) mm")
end
