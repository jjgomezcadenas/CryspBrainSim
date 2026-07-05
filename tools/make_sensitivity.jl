# tools/make_sensitivity.jl — build and cache the unscaled sensitivity base
# for the crysp_ring_1m arm on the provisional activity grid, with the real
# phantom attenuation from the products tree. With `--check`, build a second
# base at an independent seed and report the Monte-Carlo mottle (the relative
# voxel-wise spread in the illuminated core) — the data behind the n_sens
# choice; rerun at the frozen grid before the sweep (the stage check in
# dev/PLAN.md).
#
# Run:  julia -t auto --project=. tools/make_sensitivity.jl [n_sens] [--check]
# Writes out/sensitivity/<name>.npz + .toml (and _seed2 with --check).

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Statistics: mean, std

const N_SENS = let a = filter(!startswith("--"), ARGS)
    isempty(a) ? 1_000_000_000 : parse(Int, a[1])       # the PLAN.md knob
end
const CHECK = "--check" in ARGS

const SCEN = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds",
                      "uniform_headep_sobp_1e8")
const N = (64, 64, 96)
const VS = (1.5f0, 1.5f0, 1.5f0)
const ORG = centered_grid(N, VS)
const CHUNK = 20_000_000

function main()
    geo = scanner_geometry(joinpath(SCEN, "crysp_ring_1m"))
    ph = phantom_attenuation(SCEN)
    atten(xs, xe) = attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes,
                                          centre=ph.centre, mu_mm_inv=ph.mu_mm_inv)
    dev = Metal.functional() ? MtlArray : identity
    println("scanner $(geo.name): r_inner $(geo.r_inner_mm) mm, " *
            "half_length $(geo.half_length_mm) mm; μ $(ph.mu_mm_inv) mm⁻¹; " *
            "n_sens $N_SENS; device $(dev === identity ? "CPU" : "Metal")")

    name = "$(geo.name)_grid$(N[1])x$(N[2])x$(N[3])_$(VS[1])mm_n$(N_SENS)"
    out = joinpath(dirname(@__DIR__), "out", "sensitivity", name)

    build(seed) = sensitivity_base(; r_inner_mm=geo.r_inner_mm,
                                   half_length_mm=geo.half_length_mm,
                                   n=N, img_origin=ORG, voxsize=VS,
                                   attenuation=atten, n_sens=N_SENS,
                                   chunk=CHUNK, seed=seed, device=dev,
                                   progress=false)

    t1 = @elapsed base = build(1)
    @printf("seed 1 built in %.1f s\n", t1)
    save_sensitivity(out, base;
                     r_inner_mm=geo.r_inner_mm, half_length_mm=geo.half_length_mm,
                     n_sens=N_SENS, chunk=CHUNK, seed=1, img_origin=ORG, voxsize=VS,
                     attenuation_meta=(route="ellipsoid", material=ph.material,
                                       mu_mm_inv=Float64(ph.mu_mm_inv),
                                       semi_axes=ph.semi_axes, centre=ph.centre))
    println("cached $out.npz (+ .toml)")

    if CHECK
        t2 = @elapsed base2 = build(2)
        @printf("seed 2 built in %.1f s\n", t2)
        core = base .> 0.25f0 * maximum(base)          # the illuminated core
        rel = (base[core] .- base2[core]) ./ (0.5f0 .* (base[core] .+ base2[core]))
        # Two independent draws differ by √2 × the per-image mottle.
        mottle = std(rel) / sqrt(2)
        @printf("MC mottle (illuminated core, %d voxels): per-image %.3f%%, pair spread %.3f%%, max pair diff %.2f%%\n",
                count(core), 100mottle, 100std(rel), 100maximum(abs.(rel)))
    end
end

main()
