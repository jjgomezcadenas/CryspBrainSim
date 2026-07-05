# tools/bench_sensitivity.jl — the step-4 scoping benchmark: measure the
# sensitivity pipeline's three stages (surface-LOR generation, analytic
# ellipsoid attenuation, backprojection on CPU and Metal) at the recipe's
# 20 M chunk, extrapolate to n_sens = 1×10⁸ and the locked 5×10⁸, and record
# the memory footprint of the pooled 174 M-event master. Writes
# out/sensitivity_scope/bench.toml and prints the table.
#
# Run:  julia -t auto --project=. tools/bench_sensitivity.jl [chunk_size]

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using TOML

const CHUNK = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20_000_000

# The crysp_ring_1m arm and the uniform head (the products-tree values).
const R_INNER = 387.0
const HALF_LENGTH = 512.0
const AXES = (72.0f0, 87.0f0, 102.0f0)
const CENTRE = (0.0f0, -30.0f0, 0.0f0)
const MU = 0.009913f0

# Provisional activity grid (frozen at the single-shard stage): the beam
# corridor, transverse voxel centres ±47.25 mm, axial z ∈ [−119.25, +23.25] mm,
# at 1.5 mm voxels. The offset origin covers the proximal activity a centered
# z-grid would clip.
const N = (64, 64, 96)
const VS = (1.5f0, 1.5f0, 1.5f0)
const ORG = (-47.25f0, -47.25f0, -119.25f0)

const N_MASTER = 174_296_897        # the pooled BGO master (10 shards)

function bench()
    println("threads: $(Threads.nthreads()); Metal functional: $(Metal.functional())")
    println("chunk = $CHUNK LORs; grid $N at $(VS[1]) mm\n")

    sc = ContinuousPET(diameter=2R_INNER, afov=2HALF_LENGTH)
    rng = MersenneTwister(1)

    # Warm-up at a small size compiles every kernel before timing.
    wxs, wxe = sample_lors(sc, 10_000; rng=rng)
    wa = attenuation_ellipsoid(wxs, wxe; semi_axes=AXES, centre=CENTRE, mu_mm_inv=MU)
    sensitivity_image(wxs, wxe, N, ORG, VS; weights=wa)
    if Metal.functional()
        Array(sensitivity_image(MtlArray(wxs), MtlArray(wxe), N, ORG, VS;
                                weights=MtlArray(wa)))
    end

    t_gen = @elapsed gxs, gxe = sample_lors(sc, CHUNK; rng=rng)
    t_att = @elapsed ga = attenuation_ellipsoid(gxs, gxe; semi_axes=AXES,
                                                centre=CENTRE, mu_mm_inv=MU)
    t_cpu = @elapsed sensitivity_image(gxs, gxe, N, ORG, VS; weights=ga)
    t_gpu = Metal.functional() ?
        (@elapsed Array(sensitivity_image(MtlArray(gxs), MtlArray(gxe), N, ORG, VS;
                                          weights=MtlArray(ga)))) : NaN

    rate(t) = CHUNK / t / 1e6
    @printf("sample_lors             %7.2f s   (%6.1f Mlors/s)\n", t_gen, rate(t_gen))
    @printf("attenuation (ellipsoid) %7.2f s   (%6.1f Mlors/s)\n", t_att, rate(t_att))
    @printf("backprojection CPU      %7.2f s   (%6.1f Mlors/s)\n", t_cpu, rate(t_cpu))
    isnan(t_gpu) ||
        @printf("backprojection Metal    %7.2f s   (%6.1f Mlors/s)\n", t_gpu, rate(t_gpu))

    per_chunk_gpu = t_gen + t_att + (isnan(t_gpu) ? t_cpu : t_gpu)
    per_chunk_cpu = t_gen + t_att + t_cpu
    println()
    for n_sens in (100_000_000, 500_000_000)
        f = n_sens / CHUNK
        @printf("n_sens = %.0e:  GPU chain %6.1f min   CPU chain %6.1f min\n",
                n_sens, f * per_chunk_gpu / 60, f * per_chunk_cpu / 60)
    end

    # Memory footprints (GiB).
    gib(b) = b / 2^30
    chunk_bytes = 2 * 3 * CHUNK * 4 + CHUNK * 4          # endpoints + weights
    master_recon = (2 * 3 + 1) * N_MASTER * 4            # xs+xe+mult on device
    master_struct = (3 * 3 * 4 + 2 * 4 + 2 * 2 * 2 + 1 + 2) * N_MASTER
    sysmem = Sys.total_memory()
    @printf("\nchunk footprint      %5.2f GiB\n", gib(chunk_bytes))
    @printf("master (recon xs/xe/mult, %d ev) %5.2f GiB\n", N_MASTER, gib(master_recon))
    @printf("master (full MCCoincidences)     %5.2f GiB\n", gib(master_struct))
    @printf("machine RAM          %5.1f GiB\n", gib(sysmem))

    out = joinpath(dirname(@__DIR__), "out", "sensitivity_scope")
    mkpath(out)
    open(joinpath(out, "bench.toml"), "w") do io
        TOML.print(io, Dict(
            "recocrysp_sha" => recocrysp_sha(),
            "threads" => Threads.nthreads(),
            "metal_functional" => Metal.functional(),
            "chunk" => CHUNK,
            "grid" => Dict("n" => collect(N), "voxsize" => Float64.(collect(VS))),
            "seconds" => Dict("sample_lors" => t_gen, "attenuation" => t_att,
                              "back_cpu" => t_cpu,
                              "back_metal" => isnan(t_gpu) ? -1.0 : t_gpu),
            "minutes_extrapolated" => Dict(
                "gpu_1e8" => 1e8 / CHUNK * per_chunk_gpu / 60,
                "gpu_5e8" => 5e8 / CHUNK * per_chunk_gpu / 60,
                "cpu_5e8" => 5e8 / CHUNK * per_chunk_cpu / 60),
            "memory_gib" => Dict("chunk" => gib(chunk_bytes),
                                 "master_recon" => gib(master_recon),
                                 "master_struct" => gib(master_struct),
                                 "machine" => gib(sysmem))))
    end
    println("\nwrote $(joinpath(out, "bench.toml"))")
end

bench()
