# drivers/ten_shards_dose.jl — the dose axis of the ten-shard stability
# study: thin every stored shard to sub-doses of the 1 Gy acquisition and
# reconstruct each realization with the frozen chain (trues selection,
# 50 MLEM iterations, whole-plane fits happen downstream in
# tools/ten_shards.py --dose-sweep).
#
# Each shard is one independent 1 Gy acquisition, so thinning shard i with
# keep-probability p = dose/1 Gy yields one independent realization of a
# `dose` acquisition — ten per dose from the ten shards, and a second thin
# seed doubles that at the low doses where the spread is largest (two thins
# at p ≤ 0.2 share almost no events). The thin applies to the FULL event
# list (the acquisition scales with dose), the trues selection follows.
#
# Seeding: realization_index = shard·1000 + 100·dose + seed inside the
# thinning namespace (THINNING_SEED_BASE) — unique per (shard, dose, seed),
# bit-for-bit reproducible.
#
# Run:  julia -t auto --project=. drivers/ten_shards_dose.jl
# One Julia session for all realizations: the sensitivity base loads once,
# each reconstruction reuses it. ~10 min on Metal.
# Writes <config>/ten_shards/recons/recon_shardNNN_d<dose>_s<seed>.npz

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf

# (dose_Gy, number of thin seeds); 1 Gy is the shards themselves.
const DOSES = [(0.5, 1), (0.2, 2), (0.1, 2)]
const TOP_DOSE = 1.0

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
    base, ph = ctx.base, ctx.phantom
    out = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                   "ten_shards", "recons")
    mkpath(out)

    dev = Metal.functional() ? MtlArray : identity
    n_done = 0
    t_all = @elapsed for (si, file) in enumerate(ctx.files)
        r = read_shard(file)
        tmask = is_true(r.coinc)
        M = length(r.coinc)
        for (dose, nseeds) in DOSES, seed in 0:nseeds-1
            target = round(Int, dose / TOP_DOSE * M)
            idx = si * 1000 + round(Int, 100 * dose) + seed
            keep = thin_lm(r.coinc, target, idx) .& tmask
            xs, xe = endpoints(r.coinc, keep)
            nev = size(xs, 2)
            mult = attenuation_ellipsoid(xs, xe; semi_axes=ph.semi_axes,
                                         centre=ph.centre, mu_mm_inv=ph.mu_mm_inv)
            sens = scaled_sensitivity(base, nev, PARAMS.n_sens)
            model = ListmodePoissonModel(dev(xs), dev(xe), dev(sens);
                                         img_origin=PARAMS.grid.img_origin,
                                         voxsize=PARAMS.grid.voxsize,
                                         mult=dev(mult))
            x = dev(Float32.(sens .> 0))
            t = @elapsed x = mlem(model, x; niter=PARAMS.niter)
            name = @sprintf("recon_shard%03d_d%g_s%d.npz", si - 1, dose, seed)
            npzwrite(joinpath(out, name),
                     Dict("image" => Array(x), "dose_Gy" => dose,
                          "seed" => seed, "n_events" => nev,
                          "realization_index" => idx))
            n_done += 1
            @printf("shard %d  d=%.1f Gy  s%d: %8d trues, %d iters in %5.1f s  (%s)\n",
                    si - 1, dose, seed, nev, PARAMS.niter, t, name)
        end
    end
    @printf("\n%d reconstructions in %.1f min\n", n_done, t_all / 60)
end

main()
