# drivers/washout_sigma_r.jl — the detected-level isotope-washout study, over
# the acquisition-start axis. Washout is applied as a per-event thinning of each
# stored shard, reconstructed with the frozen chain — the in-Geant4 washout
# realised downstream and event-exact (md/washout-g4-formulation.md). Each
# recorded decay is kept with probability
#
#   w(z0, t_decay) = Σ_i P(i | z0, t_decay) · g_i(t_start),
#   P(i | z0, t_decay) ∝ P_i(z0) · λ_i e^{−λ_i t_decay}/(1 − e^{−λ_i t_meas}),
#
# where P_i(z) is the truth per-isotope depth profile and (z0, t_decay) the
# event's true origin depth and decay time — both in the shard, so no isotope
# tag and no upstream data. The survival factor g_i is recomputed for the
# EFFECTIVE window [t1 + t_start, t2]: a delayed start keeps the later-decaying,
# more-washed events, so washout and start-time delay compound on ¹⁵O.
#
# Working protocol (all events). For each t_start: cut t_decay ≥ t_start, thin
# by w, reconstruct the ten shards, and read σ_R and the Δ_R50 shift against the
# matched nominal (delayed-start, no washout) from ten_shards/tstart_<T>.toml.
#
# Run:  julia -t auto --project=. drivers/washout_sigma_r.jl [t_start ...]
#       (default t_start = 60 120 180 300 s)
# Writes <config>/washout/sigma_r_washout_tstart.toml

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
const WASHOUT_SEED_BASE = 2_000_000
const LN2 = log(2.0)
const DEV = Metal.functional() ? MtlArray : identity
const THINNED = "--thinned" in ARGS
const WASHED_ONLY = "--washed-only" in ARGS   # skip the nominal fit (low-count probe)
const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 50 : parse(Int, ARGS[i+1])
end
const T_STARTS = let i = findfirst(==("--tstart"), ARGS)
    i !== nothing ? parse.(Float64, split(ARGS[i+1], ",")) :
        (THINNED ? [0.0, 60.0, 120.0, 180.0, 300.0] : [60.0, 120.0, 180.0, 300.0])
end
# The dose fraction for thinned realizations (default 1 Gy = the shards' dose).
const THIN_DOSE = let i = findfirst(==("--dose"), ARGS)
    i === nothing ? 1.0 : parse(Float64, ARGS[i+1])
end

# The elementary window/build-up integral and the survival factor (note Eqs.
# 6–7), recomputed here so the effective window can shift with t_start.
phi(a, tirr, t1, t2) = (expm1(a * tirr) * (exp(-a * t1) - exp(-a * t2))) / a^2
gfac(λ, M, μ, tirr, t1, t2) =
    sum(Mk * phi(λ + μk, tirr, t1, t2) for (Mk, μk) in zip(M, μ)) / phi(λ, tirr, t1, t2)

function interp_clamped(xg, yg, xq)
    out = similar(xq, Float64)
    n = length(xg)
    @inbounds for k in eachindex(xq)
        x = xq[k]
        if x <= xg[1]
            out[k] = yg[1]
        elseif x >= xg[n]
            out[k] = yg[n]
        else
            j = searchsortedlast(xg, x)
            t = (x - xg[j]) / (xg[j+1] - xg[j])
            out[k] = (1 - t) * yg[j] + t * yg[j+1]
        end
    end
    out
end

function washout_config()
    wp = TOML.parsefile(joinpath(REPO, "config", "washout_brain.toml"))
    isos = String.(wp["physical"]["isotopes"])
    λ = Dict(iso => LN2 / T for (iso, T) in zip(isos, wp["physical"]["T_half_s"]))
    M = Float64.(wp["model"]["M"])
    μ = [LN2 / T for T in wp["model"]["T_s"]]
    tirr, t1, t2 = wp["timing"]["t_irr_s"], wp["timing"]["t1_s"], wp["timing"]["t2_s"]
    csv = readdlm(joinpath(ROOT, SCENARIO, "truth", "activity_profile_fast.csv"),
                  ',', header=true)
    data, hdr = csv[1], vec(csv[2])
    col(name) = Float64.(data[:, findfirst(==(name), hdr)])
    return (isos=isos, λ=λ, M=M, μ=μ, tirr=tirr, t1=t1, t2=t2, tmeas=t2 - t1,
            zprof=col("z_mm"), P=Dict(iso => col(iso) for iso in isos))
end

# Per-event keep-probability w(z0, t_decay) for survival factors g.
function washout_weight(cfg, g, z0, td)
    num = zeros(Float64, length(z0))
    den = zeros(Float64, length(z0))
    for iso in cfg.isos
        λi = cfg.λ[iso]
        Pz = interp_clamped(cfg.zprof, cfg.P[iso], z0)
        norm = 1.0 - exp(-λi * cfg.tmeas)
        @inbounds for k in eachindex(z0)
            rate = Pz[k] * λi * exp(-λi * td[k]) / norm
            den[k] += rate
            num[k] += rate * g[iso]
        end
    end
    @inbounds for k in eachindex(num)
        num[k] = den[k] > 0 ? num[k] / den[k] : 0.0
    end
    num
end

# Matched nominal (delayed-start, no washout) per-shard R50, from the
# acquisition-start study's outputs.
function nominal_r50(cfgdir, t_start)
    if t_start == 0
        t = TOML.parsefile(joinpath(cfgdir, "ten_shards", "results.toml"))
        return Float64.(t["all_ev"]["erfc"]["R50_mm"])
    end
    t = TOML.parsefile(joinpath(cfgdir, "ten_shards", "tstart_$(Int(t_start)).toml"))
    return Float64.(t["all_ev"]["R50_mm"])
end

# The thinned firm-up: pool the shards and, at each t_start, draw REALIZATIONS
# nominal-thinned and washout-thinned 1 Gy acquisitions (per-event keep-prob
# 1/n_shards, times w for the washed set), so σ_R and its inflation are pinned
# to ±1/√(2(N−1)) with matched statistics on both selections.
function thinned_sweep(ctx, cfg, cfgdir)
    println("pooling $(length(ctx.files)) shards…")
    pool = pool_shards(ctx.files)
    pool.n_dropped == 0 ||
        error("thinned washout needs zero dropped LORs; pool dropped $(pool.n_dropped)")
    z0 = Float64.(pool.coinc.origin[3, :])
    td = Float64.(vcat([shard_t_decay(f) for f in ctx.files]...))
    length(td) == length(pool.coinc) || error("t_decay/pool length mismatch")
    a_all = lor_attenuation(ctx, pool.coinc.xstart, pool.coinc.xend)
    n_shards = length(ctx.files)
    M = length(pool.coinc)
    p_dose = THIN_DOSE / n_shards          # keep-prob for a THIN_DOSE-Gy realization
    band = 1 / sqrt(2 * (REALIZATIONS - 1))
    @printf("dose %.2g Gy per realization\n", THIN_DOSE)
    @printf("pooled %d coincidences; %d t_start × %d realizations × 2 selections\n",
            M, length(T_STARTS), REALIZATIONS)

    points = []
    for ts in T_STARTS
        g = Dict(iso => gfac(cfg.λ[iso], cfg.M, cfg.μ, cfg.tirr, cfg.t1 + ts, cfg.t2)
                 for iso in cfg.isos)
        w = washout_weight(cfg, g, z0, td)
        tcut = td .>= ts
        f_wo = sum(w[tcut]) / count(tcut)          # washout survival among cut events
        nom, wo = Float64[], Float64[]
        wnev = Int[]
        t = @elapsed for z in 1:REALIZATIONS
            if !WASHED_ONLY
                rn = MersenneTwister(WASHOUT_SEED_BASE + 1_000_000 + Int(ts) * 1000 + z)
                kn = tcut .& (rand(rn, M) .< p_dose)
                xs, xe = endpoints(pool.coinc, kn)
                push!(nom, reconstruct_endpoint(ctx, xs, xe, a_all[kn]; device=DEV).r50_fit)
            end
            rw = MersenneTwister(WASHOUT_SEED_BASE + 2_000_000 + Int(ts) * 1000 + z)
            kw = tcut .& (rand(rw, M) .< p_dose .* w)
            xs, xe = endpoints(pool.coinc, kw)
            res = reconstruct_endpoint(ctx, xs, xe, a_all[kw]; device=DEV)
            push!(wo, res.r50_fit); push!(wnev, res.nev)
        end
        # Failure guard: an unstable fit returns a non-finite R50; drop it from
        # σ_R and report the count so a low-dose corner cannot pass silently.
        nom_ok, wo_ok = filter(isfinite, nom), filter(isfinite, wo)
        nfail_n, nfail_w = length(nom) - length(nom_ok), REALIZATIONS - length(wo_ok)
        σn = length(nom_ok) > 1 ? std(nom_ok) : NaN
        σw = length(wo_ok) > 1 ? std(wo_ok) : NaN
        shift = isempty(nom_ok) ? NaN : mean(wo_ok) - mean(nom_ok)
        push!(points, (t_start=ts, f_wo=f_wo, σ_nom=σn, σ_wo=σw,
                       nfail_n=nfail_n, nfail_w=nfail_w, wnev_min=minimum(wnev),
                       shift=shift, nom=copy(nom), wo=copy(wo)))
        @printf("t_start %3.0f s: nom σ_R %.3f | wash σ_R %.3f | infl %.2f (±%.0f%%) | ΔR50 %+.3f mm | fails n/w %d/%d | min wash nev %d | %.0fs\n",
                ts, σn, σw, σw / σn, 100band, shift,
                nfail_n, nfail_w, minimum(wnev), t)
    end

    out = joinpath(cfgdir, "washout")
    mkpath(out)
    open(joinpath(out, "sigma_r_washout_thinned.toml"), "w") do io
        TOML.print(io, Dict(
            "method" => "thinned realizations of the pooled master",
            "dose_Gy" => THIN_DOSE,
            "realizations" => REALIZATIONS, "sigma_band" => band,
            "seed_base" => WASHOUT_SEED_BASE,
            "point" => [Dict(
                "t_start_s" => p.t_start, "washout_survival" => p.f_wo,
                "nominal_sigma_R_mm" => p.σ_nom, "washed_sigma_R_mm" => p.σ_wo,
                "inflation" => p.σ_wo / p.σ_nom,
                "n_fail_nominal" => p.nfail_n, "n_fail_washed" => p.nfail_w,
                "min_washed_events" => p.wnev_min,
                "delta_R50_washout_mm" => p.shift,
                "nominal_R50_mm" => p.nom, "washed_R50_mm" => p.wo)
                for p in points]))
    end
    println("wrote $(joinpath(out, "sigma_r_washout_thinned.toml"))")
end

function main()
    cache = joinpath(sensitivity_out(SCENARIO, TOPOLOGY, RING),
                     sensitivity_cache_name(PARAMS))
    ctx = load_run_context(; products_root=ROOT, scenario=SCENARIO,
                           topology=TOPOLOGY, scanner=RING, crystal=CFG.crystal,
                           leaf=CFG.leaf, sens_cache=cache, params=PARAMS)
    cfg = washout_config()
    cfgdir = config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal)

    # Cross-check the Julia g_i against the Python washout.toml (t_start = 0).
    let wt = TOML.parsefile(joinpath(scenario_out(SCENARIO), "washout", "washout.toml"))
        for iso in cfg.isos
            gj = gfac(cfg.λ[iso], cfg.M, cfg.μ, cfg.tirr, cfg.t1, cfg.t2)
            abs(gj - wt["g_factor"][iso]) < 1e-6 ||
                error("g_i mismatch $iso: julia $gj vs python $(wt["g_factor"][iso])")
        end
        println("g_i cross-check vs washout.toml (t_start=0): OK")
    end

    if THINNED
        thinned_sweep(ctx, cfg, cfgdir)
        return
    end

    # Precompute each shard's origins, decay times, and per-isotope P_i(z0).
    println("loading shards…")
    shards = map(ctx.files) do file
        r = read_shard(file)
        r.n_dropped == 0 ||
            error("washout thinning needs zero dropped LORs; $(basename(file)) dropped $(r.n_dropped)")
        (coinc=r.coinc, real=Int(r.attrs["realization"]),
         z0=Float64.(r.coinc.origin[3, :]), td=Float64.(shard_t_decay(file)))
    end

    points = []
    for ts in T_STARTS
        g = Dict(iso => gfac(cfg.λ[iso], cfg.M, cfg.μ, cfg.tirr, cfg.t1 + ts, cfg.t2)
                 for iso in cfg.isos)
        r50 = Float64[]
        fracs = Float64[]
        t = @elapsed for sh in shards
            w = washout_weight(cfg, g, sh.z0, sh.td)
            rng = MersenneTwister(WASHOUT_SEED_BASE + Int(ts) * 100 + sh.real)
            keep = (sh.td .>= ts) .& (rand(rng, length(w)) .< w)
            xs, xe = endpoints(sh.coinc, keep)
            res = reconstruct_endpoint(ctx, xs, xe, lor_attenuation(ctx, xs, xe); device=DEV)
            push!(r50, res.r50_fit)
            push!(fracs, count(keep) / length(w))
        end
        nom = nominal_r50(cfgdir, ts)
        σ_wo, σ_nom = std(r50), std(nom)
        push!(points, (t_start=ts, g=g, r50=r50, nom=nom, kept=mean(fracs),
                       σ_wo=σ_wo, σ_nom=σ_nom, shift=mean(r50) - mean(nom)))
        @printf("t_start %3.0f s: g(O15) %.3f | kept %.3f | washed R50 %8.3f σ_R %.3f | nominal σ_R %.3f | ΔR50 %+.3f mm | %.0fs\n",
                ts, g["O15"], mean(fracs), mean(r50), σ_wo, σ_nom, mean(r50) - mean(nom), t)
    end

    out = joinpath(cfgdir, "washout")
    mkpath(out)
    open(joinpath(out, "sigma_r_washout_tstart.toml"), "w") do io
        TOML.print(io, Dict(
            "selection" => "all-events (working protocol)",
            "n_shards" => length(shards), "seed_base" => WASHOUT_SEED_BASE,
            "point" => [Dict(
                "t_start_s" => p.t_start, "kept_fraction" => p.kept,
                "g_O15" => p.g["O15"],
                "washed_R50_mean_mm" => mean(p.r50), "washed_sigma_R_mm" => p.σ_wo,
                "washed_R50_per_shard_mm" => p.r50,
                "nominal_R50_mean_mm" => mean(p.nom), "nominal_sigma_R_mm" => p.σ_nom,
                "delta_R50_washout_mm" => p.shift) for p in points]))
    end
    println("wrote $(joinpath(out, "sigma_r_washout_tstart.toml"))")
end

main()
