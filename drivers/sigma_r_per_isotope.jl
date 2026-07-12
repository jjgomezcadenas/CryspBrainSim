# drivers/sigma_r_per_isotope.jl — σ_R measured for one isotope at a time, to
# test whether positron range (longest for ¹⁵O, shortest for ¹¹C) inflates the
# per-event position noise and hence σ_R. The data are isotope-blind, so we
# select an "isotope-i-like" sample by thinning the pooled master with that
# isotope's posterior
#
#   P(i | z0, t_decay) = r_i / Σ_j r_j,   r_i = P_i(z0) λ_i e^{-λ_i t_decay}/(1-e^{-λ_i t_meas})
#
# — the same per-event rates the washout weight marginalises, here kept per
# isotope. Each isotope enters with its OWN natural statistics (× p_dose, no
# count matching): ¹⁵O ~3× the counts of ¹¹C, so counting alone predicts
# σ_R(¹⁵O) ≈ 0.6 σ_R(¹¹C); range would erode that advantage. The posterior is a
# soft selector (leakage between isotopes dilutes the contrast → any measured
# difference is a lower bound).
#
# Uses only (z0, t_decay) + the truth profiles — no isotope tag, no upstream.
# Optional washout thinning (× w) and a t_decay ≥ t_start cut for phase 2.
#
# Run (per active arm, e.g. ring 1 m CsI at t=0):
#   julia -t auto --project=. drivers/sigma_r_per_isotope.jl [--tstart 0]
#        [--dose 1.0] [--realizations 100] [--isotopes O15,C11,N13,O14] [--washed]
# Writes <config>/washout/sigma_r_per_isotope.toml
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
const SEED_BASE = 5_000_000

const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 100 : parse(Int, ARGS[i+1])
end
const T_STARTS = let i = findfirst(==("--tstart"), ARGS)
    i !== nothing ? parse.(Float64, split(ARGS[i+1], ",")) : [0.0]
end
const THIN_DOSE = let i = findfirst(==("--dose"), ARGS)
    i === nothing ? 1.0 : parse(Float64, ARGS[i+1])
end
const ISOTOPES = let i = findfirst(==("--isotopes"), ARGS)
    i !== nothing ? String.(split(ARGS[i+1], ",")) : ["O15", "C11", "N13", "O14"]
end
const WASHED = "--washed" in ARGS

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

# per-event rates r_i and their sum; posteriors are r_i ./ den
function isotope_rates(cfg, z0, td)
    r = Dict(iso => zeros(Float64, length(z0)) for iso in cfg.isos)
    den = zeros(Float64, length(z0))
    for iso in cfg.isos
        λi = cfg.λ[iso]; Pz = interp_clamped(cfg.zprof, cfg.P[iso], z0)
        norm = 1.0 - exp(-λi * cfg.tmeas)
        ri = r[iso]
        @inbounds for k in eachindex(z0)
            ri[k] = Pz[k] * λi * exp(-λi * td[k]) / norm
            den[k] += ri[k]
        end
    end
    r, den
end

# washout marginalised weight w(z0,t_decay) for the effective window [t1+ts, t2]
function washout_weight(cfg, ts, z0, td, rates, den)
    g = Dict(iso => gfac(cfg.λ[iso], cfg.M, cfg.μ, cfg.tirr, cfg.t1 + ts, cfg.t2) for iso in cfg.isos)
    w = zeros(Float64, length(z0))
    for iso in cfg.isos
        gi = g[iso]; ri = rates[iso]
        @inbounds for k in eachindex(w); w[k] += ri[k] * gi; end
    end
    @inbounds for k in eachindex(w); w[k] = den[k] > 0 ? w[k] / den[k] : 0.0; end
    w
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
    rates, den = isotope_rates(cfg, z0, td)
    @printf("pooled %d coincidences; dose %.2g Gy; %s; isotopes %s; N=%d\n",
            M, THIN_DOSE, WASHED ? "WASHED" : "nominal", join(ISOTOPES, ","), REALIZATIONS)

    points = []
    for ts in T_STARTS
        tcut = td .>= ts
        w = WASHED ? washout_weight(cfg, ts, z0, td, rates, den) : ones(Float64, M)
        for (ii, iso) in enumerate(ISOTOPES)
            post = iso == "ALL" ? ones(Float64, M) : rates[iso] ./ den   # ALL = combined mix
            keepp = p_dose .* post .* w              # natural isotope-i statistics (× w if washed)
            r50 = Float64[]; nevs = Int[]
            t = @elapsed for z in 1:REALIZATIONS
                rng = MersenneTwister(SEED_BASE + ii * 100_000 + Int(ts) * 1000 + z)
                keep = tcut .& (rand(rng, M) .< keepp)
                xs, xe = endpoints(pool.coinc, keep)
                res = reconstruct_endpoint(ctx, xs, xe, a_all[keep]; device=DEV)
                push!(r50, res.r50_fit); push!(nevs, res.nev)
            end
            ok = filter(isfinite, r50)
            σ = length(ok) > 1 ? std(ok) : NaN
            push!(points, (t_start=ts, iso=iso, σ=σ, meanR=mean(ok), nev=mean(nevs),
                           nfail=REALIZATIONS - length(ok)))
            @printf("t_start %3.0f %s %-4s: nev %8.0f | σ_R %.3f | mean R50 %8.3f | fails %d | %.0fs\n",
                    ts, WASHED ? "wash" : "nom ", iso, mean(nevs), σ, mean(ok),
                    REALIZATIONS - length(ok), t)
        end
    end

    out = joinpath(cfgdir, "washout"); mkpath(out)
    # one file per (t_start set, selection) so runs accumulate instead of clobbering
    tag = "t" * join(string.(Int.(round.(T_STARTS))), "_")
    fname = WASHED ? "sigma_r_per_isotope_washed_$(tag).toml" : "sigma_r_per_isotope_$(tag).toml"
    open(joinpath(out, fname), "w") do io
        TOML.print(io, Dict(
            "scanner" => RING, "crystal" => CFG.crystal, "dose_Gy" => THIN_DOSE,
            "washed" => WASHED, "realizations" => REALIZATIONS,
            "point" => [Dict("t_start_s" => p.t_start, "isotope" => p.iso,
                             "mean_events" => p.nev, "sigma_R_mm" => p.σ,
                             "mean_R50_mm" => p.meanR, "n_fail" => p.nfail) for p in points]))
    end
    println("wrote $(joinpath(out, fname))")
end

main()
