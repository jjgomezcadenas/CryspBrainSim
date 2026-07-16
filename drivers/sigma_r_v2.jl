# drivers/sigma_r_v2.jl — the generation-2 σ_R study for one scanner arm, over
# the three acquisition-scenario leaves (del120/180/300, each a fixed 300 s
# window on the irradiation-end clock — the delay axis is now the leaf axis, not
# a t_decay cut). Two deliverables from one pooling pass per leaf:
#
#   1. Washout σ_R (working protocol, all-events): nominal vs EXACT washed. v2
#      carries the emitting isotope per LOR and the stamped per-isotope survival
#      washout_g, so washout is the exact per-species Bernoulli keep of
#      washout_brain.tex §5 — keep event e with probability p_dose · g_i[iso(e)],
#      not the label-free marginalised posterior. Reports σ_R nominal/washed, the
#      inflation, and the calibratable edge shift ΔR50 = mean(washed) − mean(nom).
#
#   2. Per-isotope σ_R: select a PURE single-species sub-sample by the isotope
#      column (exact, no posterior leakage) at its natural abundance (× p_dose),
#      reconstruct, take σ_R — the positron-range test, now exact. O15 and C11.
#
# Each leaf's stamped washout_g is cross-checked against a recompute from the
# stamped Mizuno params + window (note Eq. 8), so a self-contradicting shard fails.
#
# Run (active arm = ring 1 m CsI v2, config/run_parameters_csi_v2.toml):
#   julia -t auto --project=. drivers/sigma_r_v2.jl [--realizations 100]
#        [--dose 1.0] [--isotopes O15,C11] [--leaves del120s_ac300s_1Gy,...] [--tend 300]
# `--tend T` models a SHORTER scan by sub-cutting a leaf window [t1,t2] to [t1,T]
# (T ≤ t2): keep only t_decay ≤ T and recompute g_i for [t1,T] — reuses the shards,
# no new products. Output filenames gain a `_t<t1>_<T>` tag so they don't clobber
# the full-window results.
# Writes <config>/washout_v2/{sigma_r_washout_v2, sigma_r_per_isotope_v2}[_t<t1>_<T>].toml

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using Statistics: mean, std
using TOML

const ROOT = joinpath(dirname(dirname(@__DIR__)), "PtCryspProds")
const PARAMS = load_run_parameters()
const CFG = PARAMS.config
const SCENARIO, TOPOLOGY, RING = CFG.scenario, CFG.topology, CFG.scanner
const DEV = Metal.functional() ? MtlArray : identity
const SEED_BASE = 7_000_000
const LN2 = log(2.0)

const REALIZATIONS = let i = findfirst(==("--realizations"), ARGS)
    i === nothing ? 100 : parse(Int, ARGS[i+1])
end
const THIN_DOSE = let i = findfirst(==("--dose"), ARGS)
    i === nothing ? 1.0 : parse(Float64, ARGS[i+1])
end
# --tend T: sub-cut each leaf's window to a shorter scan ending at T seconds.
const TEND = let i = findfirst(==("--tend"), ARGS)
    i === nothing ? nothing : parse(Float64, ARGS[i+1])
end
# Isotope ids follow the shard's `isotope_names` order (0=O15,1=C11,2=N13,3=C10,4=O14).
const ISO_ID = Dict("O15" => 0, "C11" => 1, "N13" => 2, "C10" => 3, "O14" => 4)
const ISOTOPES = let i = findfirst(==("--isotopes"), ARGS)
    i !== nothing ? String.(split(ARGS[i+1], ",")) : ["O15", "C11"]
end
const LEAVES = let i = findfirst(==("--leaves"), ARGS)
    i !== nothing ? String.(split(ARGS[i+1], ",")) :
        ["del120s_ac300s_1Gy", "del180s_ac300s_1Gy", "del300s_ac300s_1Gy"]
end

# note Eq. 8: window/build-up integral and the closed-form per-isotope survival.
phi(a, tirr, t1, t2) = (expm1(a * tirr) * (exp(-a * t1) - exp(-a * t2))) / a^2
gfac(λ, M, μ, tirr, t1, t2) =
    sum(Mk * phi(λ + μk, tirr, t1, t2) for (Mk, μk) in zip(M, μ)) / phi(λ, tirr, t1, t2)

σ_of(v) = (ok = filter(isfinite, v); length(ok) > 1 ? std(ok) : NaN)
nfail(v) = count(!isfinite, v)

# One leaf: pool, verify washout_g, run nominal/washed/per-isotope realizations.
function run_leaf(leaf)
    cache = joinpath(sensitivity_out(SCENARIO, TOPOLOGY, RING), sensitivity_cache_name(PARAMS))
    ctx = load_run_context(; products_root=ROOT, scenario=SCENARIO, topology=TOPOLOGY,
                           scanner=RING, crystal=CFG.crystal, leaf=leaf,
                           sens_cache=cache, params=PARAMS)
    a = shard_attrs(ctx.files[1])
    t_del = Float64(a["t_del_s"])
    t1, t2, tirr = Float64(a["t1_s"]), Float64(a["t2_s"]), Float64(a["t_irr_s"])
    g_stamped = Float64.(a["washout_g"])         # per-isotope, id order, for the stamped [t1,t2]
    Thalf = Float64.(a["isotope_half_lives"])
    M_k = Float64.(a["washout_fractions"]); T_k = Float64.(a["washout_Thalf_s"])
    μ = LN2 ./ T_k
    # integrity: stamped washout_g == recompute from stamped Mizuno + stamped window
    for i in eachindex(g_stamped)
        gj = gfac(LN2 / Thalf[i], M_k, μ, tirr, t1, t2)
        abs(gj - g_stamped[i]) < 1e-6 ||
            error("$leaf: washout_g[$i] $(g_stamped[i]) ≠ recomputed $gj")
    end
    # a shorter scan (--tend) sub-cuts the window to [t1, t2_eff]; recompute g_i there.
    t2_eff = TEND === nothing ? t2 : TEND
    t2_eff <= t2 || error("$leaf: --tend $t2_eff exceeds the leaf window end $t2")
    g = TEND === nothing ? g_stamped :
        [gfac(LN2 / Thalf[i], M_k, μ, tirr, t1, t2_eff) for i in eachindex(Thalf)]
    t_ac = t2_eff - t1

    println("── $leaf  window [$t1,$t2_eff]s  g(O15,C11)=($(round(g[1],digits=3)),$(round(g[2],digits=3)))")
    pool = pool_shards(ctx.files)
    pool.n_dropped == 0 || error("$leaf: pool dropped $(pool.n_dropped) LORs")
    iso = vcat([shard_isotope(f) for f in ctx.files]...)
    length(iso) == length(pool.coinc) || error("$leaf: isotope/pool length mismatch")
    td = Float64.(vcat([shard_t_decay(f) for f in ctx.files]...))
    length(td) == length(pool.coinc) || error("$leaf: t_decay/pool length mismatch")
    inwin = TEND === nothing ? trues(length(iso)) : (td .<= t2_eff)
    a_all = lor_attenuation(ctx, pool.coinc.xstart, pool.coinc.xend)
    n_shards = length(ctx.files); Mtot = length(pool.coinc)
    p_dose = THIN_DOSE / n_shards
    gvec = [g[Int(i) + 1] for i in iso]          # per-event survival = its isotope's g_i
    @printf("   pooled %d coincs (%d in [%.0f,%.0f]s); %d shards; p_dose %.3f; N=%d\n",
            Mtot, count(inwin), t1, t2_eff, n_shards, p_dose, REALIZATIONS)

    recon(keep) = reconstruct_endpoint(ctx, endpoints(pool.coinc, keep)...,
                                       a_all[keep]; device=DEV).r50_fit
    reconv(keep) = (r = reconstruct_endpoint(ctx, endpoints(pool.coinc, keep)...,
                                             a_all[keep]; device=DEV); (r.r50_fit, r.nev))

    # (1) washout: nominal vs exact washed
    nom, wsh = Float64[], Float64[]
    t = @elapsed for z in 1:REALIZATIONS
        rn = MersenneTwister(SEED_BASE + Int(t_del) * 1000 + z)
        push!(nom, recon(inwin .& (rand(rn, Mtot) .< p_dose)))
        rw = MersenneTwister(SEED_BASE + 500_000 + Int(t_del) * 1000 + z)
        push!(wsh, recon(inwin .& (rand(rw, Mtot) .< p_dose .* gvec)))
    end
    σn, σw = σ_of(nom), σ_of(wsh)
    survival = sum(gvec[inwin]) / count(inwin)
    shift = mean(filter(isfinite, wsh)) - mean(filter(isfinite, nom))
    @printf("   washout: nom σ_R %.3f | washed σ_R %.3f | infl %.2f | ΔR50 %+.3f | surv %.3f | fails n/w %d/%d | %.0fs\n",
            σn, σw, σw / σn, shift, survival, nfail(nom), nfail(wsh), t)

    # (2) per-isotope: pure single-species sub-sample at natural abundance
    iso_pts = []
    for name in ISOTOPES
        id = ISO_ID[name]; sel = (iso .== id) .& inwin
        r50 = Float64[]; nevs = Int[]
        ti = @elapsed for z in 1:REALIZATIONS
            ri = MersenneTwister(SEED_BASE + 900_000 + id * 100_000 + Int(t_del) * 1000 + z)
            r, nv = reconv(sel .& (rand(ri, Mtot) .< p_dose))
            push!(r50, r); push!(nevs, nv)
        end
        push!(iso_pts, (iso=name, σ=σ_of(r50), meanR=mean(filter(isfinite, r50)),
                        nev=mean(nevs), nfail=nfail(r50)))
        @printf("   iso %-4s: nev %8.0f | σ_R %.3f | mean R50 %8.3f | fails %d | %.0fs\n",
                name, mean(nevs), σ_of(r50), mean(filter(isfinite, r50)), nfail(r50), ti)
    end

    return (leaf=leaf, t_del=t_del, t_ac=t_ac, t2=t2_eff, g=g, survival=survival,
            crystal_label=ctx.crystal, σ_nom=σn, σ_wsh=σw, shift=shift,
            R_nom=mean(filter(isfinite, nom)), R_wsh=mean(filter(isfinite, wsh)),
            nfail_n=nfail(nom), nfail_w=nfail(wsh), iso=iso_pts)
end

function main()
    @printf("v2 σ_R study — %s / %s, dose %.2g Gy, N=%d\n", RING, CFG.crystal, THIN_DOSE, REALIZATIONS)
    results = [run_leaf(leaf) for leaf in LEAVES]
    # the crystal LABEL dir (e.g. csi_2X0) — the tree convention shared with the tools
    cfgdir = config_out(SCENARIO, TOPOLOGY, RING, results[1].crystal_label)
    out = joinpath(cfgdir, "washout_v2"); mkpath(out)
    # tag a shorter-scan run so it doesn't clobber the full-window results
    tag = TEND === nothing ? "" :
        "_t$(Int(round(results[1].t_del)))_$(Int(round(TEND)))"

    open(joinpath(out, "sigma_r_washout_v2$(tag).toml"), "w") do io
        TOML.print(io, Dict(
            "generation" => "v2", "scanner" => RING, "crystal" => CFG.crystal,
            "dose_Gy" => THIN_DOSE, "realizations" => REALIZATIONS,
            "sigma_band" => 1 / sqrt(2 * (REALIZATIONS - 1)),
            "method" => "exact per-species g_i keep on the pooled leaf",
            "point" => [Dict(
                "leaf" => r.leaf, "t_del_s" => r.t_del, "t_ac_s" => r.t_ac,
                "t2_s" => r.t2, "washout_survival" => r.survival,
                "nominal_sigma_R_mm" => r.σ_nom, "washed_sigma_R_mm" => r.σ_wsh,
                "inflation" => r.σ_wsh / r.σ_nom,
                "delta_R50_washout_mm" => r.shift,
                "nominal_R50_mean_mm" => r.R_nom, "washed_R50_mean_mm" => r.R_wsh,
                "n_fail_nominal" => r.nfail_n, "n_fail_washed" => r.nfail_w)
                for r in results]))
    end
    open(joinpath(out, "sigma_r_per_isotope_v2$(tag).toml"), "w") do io
        TOML.print(io, Dict(
            "generation" => "v2", "scanner" => RING, "crystal" => CFG.crystal,
            "dose_Gy" => THIN_DOSE, "realizations" => REALIZATIONS,
            "method" => "pure per-species selection by the isotope column",
            "point" => [Dict(
                "leaf" => r.leaf, "t_del_s" => r.t_del, "t_ac_s" => r.t_ac,
                "isotope" => p.iso, "sigma_R_mm" => p.σ, "mean_R50_mm" => p.meanR,
                "mean_events" => p.nev, "n_fail" => p.nfail)
                for r in results for p in r.iso]))
    end
    println("wrote $(joinpath(out, "sigma_r_washout_v2$(tag).toml"))")
    println("wrote $(joinpath(out, "sigma_r_per_isotope_v2$(tag).toml"))")
end

main()
