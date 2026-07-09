# tools/recon_scatters.jl — the scatter-pedestal diagnostic: reconstruct one
# shard's SCATTER coincidences alone (truth flag 1) with the frozen chain
# (grid, attenuation, scaled sensitivity, 50 MLEM iterations). The resulting
# image shows where the reconstruction puts the mispositioned events — the
# shape that sits under the distal edge when all events are admitted
# uncorrected. tools/scatter_profile.py turns it into the profile figure and
# the numbers (pedestal level, slope through the fit window).
#
# Run:  julia -t auto --project=. tools/recon_scatters.jl [shard_index]
# Writes <config>/one_shard/recon_shardNNN_scatters.npz

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf

const SHARD = isempty(ARGS) ? 0 : parse(Int, ARGS[1])

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")

const PARAMS = load_run_parameters()
const CFG = PARAMS.config      # the active arm (run_parameters.toml)
const SCENARIO, TOPOLOGY, RING = CFG.scenario, CFG.topology, CFG.scanner

function main()
    cache = joinpath(sensitivity_out(SCENARIO, TOPOLOGY, RING),
                     sensitivity_cache_name(PARAMS))
    ctx = load_run_context(; products_root=ROOT, scenario=SCENARIO,
                           topology=TOPOLOGY, scanner=RING, crystal=CFG.crystal,
                           leaf=CFG.leaf, sens_cache=cache, params=PARAMS)
    ph = ctx.phantom

    r = read_shard(ctx.files[SHARD+1])
    mask = is_scatter(r.coinc)
    xs, xe = endpoints(r.coinc, mask)
    nev = size(xs, 2)
    @printf("shard %d: %d scatters of %d coincidences (%.1f%%)\n",
            SHARD, nev, length(r.coinc), 100nev / length(r.coinc))

    mult = attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes,
                                 centre=ph.centre, mu_mm_inv=ph.mu_mm_inv)
    sens = scaled_sensitivity(ctx.base, nev, PARAMS.n_sens)
    dev = Metal.functional() ? MtlArray : identity
    model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                 img_origin=PARAMS.grid.img_origin,
                                 voxsize=PARAMS.grid.voxsize, mult=dev(mult))
    x = dev(Float32.(sens .> 0))
    t = @elapsed x = mlem(model, x; niter=PARAMS.niter)

    out = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                   "one_shard")
    mkpath(out)
    name = @sprintf("recon_shard%03d_scatters.npz", SHARD)
    npzwrite(joinpath(out, name),
             Dict("image" => Array(x), "n_events" => nev))
    @printf("%d MLEM iters in %.1f s; wrote %s\n", PARAMS.niter, t,
            joinpath(out, name))
end

main()
