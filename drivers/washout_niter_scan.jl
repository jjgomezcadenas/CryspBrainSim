# drivers/washout_niter_scan.jl — σ_R as a function of MLEM iterations, for the
# washed thinned realizations, to find the niter that minimises σ_R per geometry.
# Motivation: at a fixed niter each geometry sits at a different point on its own
# MLEM noise–resolution curve, so σ_R(niter) is expected to be U-shaped — high at
# low niter (the distal edge is not yet formed, the erfc fit is unstable) and high
# at large niter (MLEM noise dominates), with a minimum in between whose location
# can differ by geometry. This maps that curve.
#
# Efficiency: MLEM is iterative, so one reconstruction to NITER_MAX records R50 at
# every checkpoint (like drivers/one_shard.jl's r50_vs_iter) — the whole curve at
# the cost of one reconstruction pass per realization, not a scan.
#
# Washed selection only (the anomaly), 1 Gy default, over the pooled master, same
# thinning as drivers/washout_sigma_r.jl --thinned. Run per active arm:
#   julia -t auto --project=. drivers/washout_niter_scan.jl [--tstart 120,180]
#        [--dose 1.0] [--realizations 50] [--niter-max 150] [--check-every 10] [--seed N]
# Writes <config>/washout/sigma_r_vs_niter.toml

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using Statistics: mean, std
using TOML
using DelimitedFiles: readdlm

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const REPO = dirname(@__DIR__)
const PARAMS = load_run_parameters()
const CFG = PARAMS.config
const SCENARIO, TOPOLOGY, RING = CFG.scenario, CFG.topology, CFG.scanner
const LN2 = log(2.0)
const DEV = Metal.functional() ? MtlArray : identity

const SEED_BASE = let i = findfirst(==("--seed"), ARGS)
    i === nothing ? 2_000_000 : parse(Int, ARGS[i+1])
end
const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 50 : parse(Int, ARGS[i+1])
end
const T_STARTS = let i = findfirst(==("--tstart"), ARGS)
    i !== nothing ? parse.(Float64, split(ARGS[i+1], ",")) : [120.0, 180.0]
end
const THIN_DOSE = let i = findfirst(==("--dose"), ARGS)
    i === nothing ? 1.0 : parse(Float64, ARGS[i+1])
end
const NITER_MAX = let i = findfirst(==("--niter-max"), ARGS)
    i === nothing ? 150 : parse(Int, ARGS[i+1])
end
const CHECK_EVERY = let i = findfirst(==("--check-every"), ARGS)
    i === nothing ? 10 : parse(Int, ARGS[i+1])
end

# --- washout machinery (mirrors drivers/washout_sigma_r.jl) ---
phi(a, tirr, t1, t2) = (expm1(a * tirr) * (exp(-a * t1) - exp(-a * t2))) / a^2
gfac(λ, M, μ, tirr, t1, t2) =
    sum(Mk * phi(λ + μk, tirr, t1, t2) for (Mk, μk) in zip(M, μ)) / phi(λ, tirr, t1, t2)

function interp_clamped(xg, yg, xq)
    out = similar(xq, Float64); n = length(xg)
    @inbounds for k in eachindex(xq)
        x = xq[k]
        if x <= xg[1]; out[k] = yg[1]
        elseif x >= xg[n]; out[k] = yg[n]
        else
            j = searchsortedlast(xg, x); t = (x - xg[j]) / (xg[j+1] - xg[j])
            out[k] = (1 - t) * yg[j] + t * yg[j+1]
        end
    end
    out
end

function washout_config()
    wp = TOML.parsefile(joinpath(REPO, "config", "washout_brain.toml"))
    isos = String.(wp["physical"]["isotopes"])
    λ = Dict(iso => LN2 / T for (iso, T) in zip(isos, wp["physical"]["T_half_s"]))
    M = Float64.(wp["model"]["M"]); μ = [LN2 / T for T in wp["model"]["T_s"]]
    tirr, t1, t2 = wp["timing"]["t_irr_s"], wp["timing"]["t1_s"], wp["timing"]["t2_s"]
    csv = readdlm(joinpath(ROOT, SCENARIO, "truth", "activity_profile_fast.csv"), ',', header=true)
    data, hdr = csv[1], vec(csv[2])
    col(name) = Float64.(data[:, findfirst(==(name), hdr)])
    return (isos=isos, λ=λ, M=M, μ=μ, tirr=tirr, t1=t1, t2=t2, tmeas=t2 - t1,
            zprof=col("z_mm"), P=Dict(iso => col(iso) for iso in isos))
end

function washout_weight(cfg, g, z0, td)
    num = zeros(Float64, length(z0)); den = zeros(Float64, length(z0))
    for iso in cfg.isos
        λi = cfg.λ[iso]; Pz = interp_clamped(cfg.zprof, cfg.P[iso], z0)
        norm = 1.0 - exp(-λi * cfg.tmeas)
        @inbounds for k in eachindex(z0)
            rate = Pz[k] * λi * exp(-λi * td[k]) / norm
            den[k] += rate; num[k] += rate * g[iso]
        end
    end
    @inbounds for k in eachindex(num)
        num[k] = den[k] > 0 ? num[k] / den[k] : 0.0
    end
    num
end

# checkpointed reconstruction: R50 at each niter checkpoint from one MLEM pass
function r50_vs_iter(ctx, xs, xe, mult; device=DEV)
    p = ctx.params
    sens = scaled_sensitivity(ctx.base, size(xs, 2), p.n_sens)
    model = ListmodePoissonModel(device(xs), device(xe), device(sens);
                                 img_origin=p.grid.img_origin, voxsize=p.grid.voxsize,
                                 mult=device(mult))
    x = device(Float32.(sens .> 0))
    iters = Int[]; r50s = Float64[]; it = 0
    while it < NITER_MAX
        x = mlem(model, x; niter=CHECK_EVERY); it += CHECK_EVERY
        z, prof = depth_profile(Array(x); voxel_size_mm=p.grid.voxsize, beam_axis=3,
                                roi_radius_mm=p.roi.radius_mm, z_origin_mm=p.grid.img_origin[3])
        fit = fit_endpoint(z, prof; window=ctx.ref.window, weighted=true)
        push!(iters, it); push!(r50s, fit.z0)
    end
    iters, r50s
end

function main()
    cache = joinpath(sensitivity_out(SCENARIO, TOPOLOGY, RING), sensitivity_cache_name(PARAMS))
    ctx = load_run_context(; products_root=ROOT, scenario=SCENARIO, topology=TOPOLOGY,
                           scanner=RING, crystal=CFG.crystal, leaf=CFG.leaf,
                           sens_cache=cache, params=PARAMS)
    cfg = washout_config()
    cfgdir = config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal)

    println("pooling $(length(ctx.files)) shards…")
    pool = pool_shards(ctx.files)
    pool.n_dropped == 0 || error("pool dropped $(pool.n_dropped) LORs")
    z0 = Float64.(pool.coinc.origin[3, :])
    td = Float64.(vcat([shard_t_decay(f) for f in ctx.files]...))
    a_all = lor_attenuation(ctx, pool.coinc.xstart, pool.coinc.xend)
    n_shards = length(ctx.files); M = length(pool.coinc)
    p_dose = THIN_DOSE / n_shards
    checkpoints = collect(CHECK_EVERY:CHECK_EVERY:NITER_MAX)
    @printf("pooled %d coincidences; dose %.2g Gy; %d t_start × %d washed realizations; niter %d..%d step %d\n",
            M, THIN_DOSE, length(T_STARTS), REALIZATIONS, CHECK_EVERY, NITER_MAX, CHECK_EVERY)

    points = []
    for ts in T_STARTS
        g = Dict(iso => gfac(cfg.λ[iso], cfg.M, cfg.μ, cfg.tirr, cfg.t1 + ts, cfg.t2) for iso in cfg.isos)
        w = washout_weight(cfg, g, z0, td)
        tcut = td .>= ts
        R = fill(NaN, REALIZATIONS, length(checkpoints))   # R50[realization, checkpoint]
        t = @elapsed for zi in 1:REALIZATIONS
            rw = MersenneTwister(SEED_BASE + 2_000_000 + Int(ts) * 1000 + zi)
            kw = tcut .& (rand(rw, M) .< p_dose .* w)
            xs, xe = endpoints(pool.coinc, kw)
            _, r50s = r50_vs_iter(ctx, xs, xe, a_all[kw])
            R[zi, :] = r50s
        end
        σ = [ (v = filter(isfinite, R[:, c]); length(v) > 1 ? std(v) : NaN) for c in eachindex(checkpoints) ]
        μR = [ (v = filter(isfinite, R[:, c]); isempty(v) ? NaN : mean(v)) for c in eachindex(checkpoints) ]
        cmin = argmin(map(x -> isnan(x) ? Inf : x, σ))
        push!(points, (t_start=ts, sigma=σ, meanR=μR))
        @printf("t_start %3.0f s: σ_R(niter) min %.3f at niter %d | σ_R@%d = %.3f | %.0fs\n",
                ts, σ[cmin], checkpoints[cmin], PARAMS.niter,
                σ[findfirst(==(PARAMS.niter), checkpoints)], t)
    end

    out = joinpath(cfgdir, "washout"); mkpath(out)
    open(joinpath(out, "sigma_r_vs_niter.toml"), "w") do io
        TOML.print(io, Dict(
            "scanner" => RING, "crystal" => CFG.crystal, "dose_Gy" => THIN_DOSE,
            "realizations" => REALIZATIONS, "frozen_niter" => PARAMS.niter,
            "niter" => checkpoints,
            "point" => [Dict("t_start_s" => p.t_start,
                             "washed_sigma_R_mm" => p.sigma,
                             "washed_mean_R50_mm" => p.meanR) for p in points]))
    end
    println("wrote $(joinpath(out, "sigma_r_vs_niter.toml"))")
end

main()
