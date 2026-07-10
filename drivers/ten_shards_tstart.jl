# drivers/ten_shards_tstart.jl — the acquisition-start-time axis of the
# ten-shard study: apply a delayed start as the pure event selection
# t_decay >= t_start on every stored shard of the active scanner and
# reconstruct each with the frozen chain, in both selections (trues and all
# events). The fits happen downstream in tools/ten_shards.py --t-start.
#
# The stored acquisition covers [0, t_meas_s]; a delayed start keeps the
# window's end, so the selection is exact — Poisson counts, per-isotope
# decay laws, and the randoms' real time structure all included. The cut
# vector is read raw from the file (shard_t_decay) and aligned with the
# coincidence list under the assertion that no degenerate LORs were dropped.
#
# Run:  julia -t auto --project=. drivers/ten_shards_tstart.jl [t_start_s]
# Writes <config>/ten_shards/recons/recon_shardNNN_t<T>[_all].npz

using CryspBrainSim
using RecoCryspTools
using Metal
using NPZ: npzwrite
using Printf

const T_START = isempty(ARGS) ? 60.0 : parse(Float64, ARGS[1])

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
    t_all = @elapsed for (si, file) in enumerate(ctx.files)
        r = read_shard(file)
        r.n_dropped == 0 ||
            error("t_decay alignment needs zero dropped LORs; $(basename(file)) dropped $(r.n_dropped)")
        keep_t = shard_t_decay(file) .>= Float32(T_START)
        for (sel, tag_sel) in ((is_true(r.coinc), ""), (trues(length(r.coinc)), "_all"))
            keep = keep_t .& sel
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
            name = @sprintf("recon_shard%03d_t%g%s.npz", si - 1, T_START, tag_sel)
            npzwrite(joinpath(out, name),
                     Dict("image" => Array(x), "t_start_s" => T_START,
                          "n_events" => nev))
            frac = nev / count(sel)
            @printf("shard %d  t>=%gs %-5s %9d events (%.1f%% kept), %d iters in %4.1f s\n",
                    si - 1, T_START, tag_sel == "" ? "trues" : "all", nev,
                    100frac, PARAMS.niter, t)
        end
    end
    @printf("\n20 reconstructions in %.1f min\n", t_all / 60)
end

main()
