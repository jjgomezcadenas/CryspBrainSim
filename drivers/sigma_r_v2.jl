# drivers/sigma_r_v2.jl — the generation-2 σ_R study for one scanner arm, over
# the three acquisition-scenario leaves (del120/180/300, each a fixed 300 s
# window on the irradiation-end clock — the delay axis is now the leaf axis, not
# a t_decay cut). Two deliverables from one pooling pass per leaf:
#
#   1. Washout σ_R (working protocol, all-events): nominal vs EXACT washed. v2
#      carries the emitting isotope per LOR and the stamped per-isotope survival
#      washout_g, so washout uses a per-species, window-integrated Bernoulli keep
#      washout_brain.tex §5 — keep event e with probability p_dose · g_i[iso(e)],
#      not the label-free marginalised posterior. Reports σ_R nominal/washed, the
#      inflation, and the calibratable edge shift ΔR50 = mean(washed) − mean(nom).
#
#   2. Per-isotope σ_R: select a PURE single-species sub-sample by the isotope
#      column (exact, no posterior leakage) at its natural abundance (× p_dose),
#      reconstruct, take σ_R — the positron-range test, now exact. O15 and C11.
#      `--isotopes none` skips this pass.
#
#   3. Grogg-estimator σ_R (Grogg et al., IEEE TNS 60 (2013) 3290): every
#      realization also carries the linear-x-intercept endpoint on the same
#      profile — Poisson-weighted primary + paper-literal unweighted — so the
#      estimator comparison to the erfc R50 is paired, realization by
#      realization. Reported nominal and washed, with the in-window O15
#      fraction (the paper's "~80% ¹⁵O" premise, measured).
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
# Writes <config>/washout_v2/{sigma_r_washout_v2, sigma_r_per_isotope_v2,
# sigma_r_grogg_v2}[_t<t1>_<T>].toml

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using Statistics: mean, median, std
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
const MAX_FAILURE_FRACTION = let i = findfirst(==("--max-fit-failure-fraction"), ARGS)
    i === nothing ? 0.05 : parse(Float64, ARGS[i+1])
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
    i === nothing ? ["O15", "C11"] :
        ARGS[i+1] == "none" ? String[] : String.(split(ARGS[i+1], ","))
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

0.0 <= MAX_FAILURE_FRACTION < 1.0 ||
    error("--max-fit-failure-fraction must lie in [0,1)")

function require_stable_fits(label, endpoints)
    failures = nfail(endpoints)
    nfits = length(endpoints)
    max_failures = floor(Int, MAX_FAILURE_FRACTION * nfits)
    failures <= max_failures ||
        error("$label: $failures/$nfits endpoint fits failed; at most " *
              "$max_failures failures are allowed for N=$nfits " *
              "($(100MAX_FAILURE_FRACTION)% limit)")
    return failures
end

function finite_summary(values)
    valid = sort(filter(isfinite, values))
    isempty(valid) && return (median=NaN, minimum=NaN, maximum=NaN)
    return (median=median(valid), minimum=first(valid), maximum=last(valid))
end

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
                                       a_all[keep]; device=DEV)
    reconv(keep) = (r = recon(keep); (r.r50_fit, r.nev))

    # (1+3) washout nominal vs exact washed — every realization carries the
    # erfc R50 and the Grogg x-intercept (weighted + unweighted), paired.
    est = Dict(k => Float64[] for k in
               (:nom, :wsh, :nom_chi2, :wsh_chi2,
                :nomg, :wshg, :nomgu, :wshgu, :nomgs, :wshgs,
                :nomg_zf, :nomg_zl, :nomg_sl, :wshg_zf, :wshg_zl, :wshg_sl,
                :nomgs_zf, :nomgs_zl, :nomgs_sl, :wshgs_zf, :wshgs_zl, :wshgs_sl))
    t = @elapsed for z in 1:REALIZATIONS
        rng = MersenneTwister(SEED_BASE + Int(t_del) * 1000 + z)
        uniform = rand(rng, Mtot)
        nominal_keep = inwin .& (uniform .< p_dose)
        washed_keep = inwin .& (uniform .< p_dose .* gvec)
        all(washed_keep .<= nominal_keep) ||
            error("$leaf: washed selection is not a subset of nominal")

        r = recon(nominal_keep)
        push!(est[:nom], r.r50_fit); push!(est[:nomg], r.rx_grogg)
        push!(est[:nom_chi2], r.erfc_chi2_dof)
        push!(est[:nomgu], r.rx_grogg_unw); push!(est[:nomgs], r.rx_grogg_sm)
        push!(est[:nomg_zf], r.grogg_z_first); push!(est[:nomg_zl], r.grogg_z_last)
        push!(est[:nomg_sl], r.grogg_slope)
        push!(est[:nomgs_zf], r.grogg_sm_z_first); push!(est[:nomgs_zl], r.grogg_sm_z_last)
        push!(est[:nomgs_sl], r.grogg_sm_slope)
        r = recon(washed_keep)
        push!(est[:wsh], r.r50_fit); push!(est[:wshg], r.rx_grogg)
        push!(est[:wsh_chi2], r.erfc_chi2_dof)
        push!(est[:wshgu], r.rx_grogg_unw); push!(est[:wshgs], r.rx_grogg_sm)
        push!(est[:wshg_zf], r.grogg_z_first); push!(est[:wshg_zl], r.grogg_z_last)
        push!(est[:wshg_sl], r.grogg_slope)
        push!(est[:wshgs_zf], r.grogg_sm_z_first); push!(est[:wshgs_zl], r.grogg_sm_z_last)
        push!(est[:wshgs_sl], r.grogg_sm_slope)
    end
    nom, wsh = est[:nom], est[:wsh]
    nfail_nominal = require_stable_fits("$leaf nominal", nom)
    nfail_washout = require_stable_fits("$leaf washout", wsh)
    correction_nominal = finite_pool_correction(p_dose)
    correction_washout = finite_pool_correction(p_dose .* gvec[inwin])
    σn_raw, σw_raw = σ_of(nom), σ_of(wsh)
    σn = correction_nominal * σn_raw
    σw = correction_washout * σw_raw
    survival = sum(gvec[inwin]) / count(inwin)
    shift = mean(filter(isfinite, wsh)) - mean(filter(isfinite, nom))
    # the measured isotope mix behind Grogg's ¹⁵O-dominance premise
    o15win = (iso .== ISO_ID["O15"]) .& inwin
    fO15 = count(o15win) / count(inwin)
    fO15w = sum(gvec[o15win]) / sum(gvec[inwin])
    @printf("   washout: nom σ_R %.3f→%.3f | washed σ_R %.3f→%.3f | infl %.2f | ΔR50 %+.3f | surv %.3f | fails n/w %d/%d | %.0fs\n",
            σn_raw, σn, σw_raw, σw, σw / σn, shift, survival,
            nfail(nom), nfail(wsh), t)
    @printf("   grogg-w: nom σ_R %.3f | washed σ_R %.3f | mean nom %8.3f | fails n/w %d/%d | O15 frac %.3f (washed %.3f)\n",
            correction_nominal * σ_of(est[:nomg]),
            correction_washout * σ_of(est[:wshg]), mean(filter(isfinite, est[:nomg])),
            nfail(est[:nomg]), nfail(est[:wshg]), fO15, fO15w)
    @printf("   grogg-u: nom σ_R %.3f | washed σ_R %.3f | mean nom %8.3f | fails n/w %d/%d\n",
            correction_nominal * σ_of(est[:nomgu]),
            correction_washout * σ_of(est[:wshgu]), mean(filter(isfinite, est[:nomgu])),
            nfail(est[:nomgu]), nfail(est[:wshgu]))
    @printf("   grogg-7mm: nom σ_R %.3f | washed σ_R %.3f | mean nom %8.3f | fails n/w %d/%d\n",
            correction_nominal * σ_of(est[:nomgs]),
            correction_washout * σ_of(est[:wshgs]), mean(filter(isfinite, est[:nomgs])),
            nfail(est[:nomgs]), nfail(est[:wshgs]))

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
        isotope_failures = require_stable_fits("$leaf isotope $name", r50)
        σ_raw = σ_of(r50)
        correction = finite_pool_correction(p_dose)
        push!(iso_pts, (iso=name, σ_raw=σ_raw, σ=correction * σ_raw,
                        correction=correction, meanR=mean(filter(isfinite, r50)),
                        nev=mean(nevs), nfail=isotope_failures))
        @printf("   iso %-4s: nev %8.0f | σ_R %.3f | mean R50 %8.3f | fails %d | %.0fs\n",
                name, mean(nevs), σ_of(r50), mean(filter(isfinite, r50)), nfail(r50), ti)
    end

    # Grogg summary per variant: σ, mean, fails, nominal and washed.
    gsum(kn, kw) = (σn_raw=σ_of(est[kn]), σw_raw=σ_of(est[kw]),
                    σn=correction_nominal * σ_of(est[kn]),
                    σw=correction_washout * σ_of(est[kw]),
                    Rn=mean(filter(isfinite, est[kn])),
                    Rw=mean(filter(isfinite, est[kw])),
                    fn=nfail(est[kn]), fw=nfail(est[kw]))
    return (leaf=leaf, t_del=t_del, t_ac=t_ac, t2=t2_eff, g=g, survival=survival,
            crystal_label=ctx.crystal, σ_nom_raw=σn_raw, σ_wsh_raw=σw_raw,
            σ_nom=σn, σ_wsh=σw,
            correction_nominal=correction_nominal,
            correction_washout=correction_washout, shift=shift,
            R_nom=mean(filter(isfinite, nom)), R_wsh=mean(filter(isfinite, wsh)),
            nfail_n=nfail_nominal, nfail_w=nfail_washout,
            chi2_nominal=finite_summary(est[:nom_chi2]),
            chi2_washout=finite_summary(est[:wsh_chi2]), iso=iso_pts,
            grogg_w=gsum(:nomg, :wshg), grogg_u=gsum(:nomgu, :wshgu),
            grogg_s=gsum(:nomgs, :wshgs),
            fO15=fO15, fO15_washed=fO15w, est=est)
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
            "method" => "paired nominal and isotope-indexed window-integrated survival thinning with common uniforms",
            "point" => [Dict(
                "leaf" => r.leaf, "t_del_s" => r.t_del, "t_ac_s" => r.t_ac,
                "t2_s" => r.t2, "washout_survival" => r.survival,
                "sigma_R_convention" => "finite-pool corrected",
                "nominal_sigma_R_raw_mm" => r.σ_nom_raw,
                "washed_sigma_R_raw_mm" => r.σ_wsh_raw,
                "nominal_sigma_R_mm" => r.σ_nom, "washed_sigma_R_mm" => r.σ_wsh,
                "nominal_finite_pool_correction" => r.correction_nominal,
                "washed_finite_pool_correction" => r.correction_washout,
                "inflation" => r.σ_wsh / r.σ_nom,
                "delta_R50_washout_mm" => r.shift,
                "nominal_R50_mean_mm" => r.R_nom, "washed_R50_mean_mm" => r.R_wsh,
                "max_fit_failure_fraction" => MAX_FAILURE_FRACTION,
                "n_fail_nominal" => r.nfail_n, "n_fail_washed" => r.nfail_w,
                "nominal_erfc_chi2_dof_median" => r.chi2_nominal.median,
                "nominal_erfc_chi2_dof_min" => r.chi2_nominal.minimum,
                "nominal_erfc_chi2_dof_max" => r.chi2_nominal.maximum,
                "washed_erfc_chi2_dof_median" => r.chi2_washout.median,
                "washed_erfc_chi2_dof_min" => r.chi2_washout.minimum,
                "washed_erfc_chi2_dof_max" => r.chi2_washout.maximum)
                for r in results]))
    end
    if !isempty(ISOTOPES)
        open(joinpath(out, "sigma_r_per_isotope_v2$(tag).toml"), "w") do io
            TOML.print(io, Dict(
                "generation" => "v2", "scanner" => RING, "crystal" => CFG.crystal,
                "dose_Gy" => THIN_DOSE, "realizations" => REALIZATIONS,
                "method" => "pure per-species selection by the isotope column",
                "point" => [Dict(
                    "leaf" => r.leaf, "t_del_s" => r.t_del, "t_ac_s" => r.t_ac,
                    "isotope" => p.iso, "sigma_R_convention" => "finite-pool corrected",
                    "sigma_R_raw_mm" => p.σ_raw, "sigma_R_mm" => p.σ,
                    "finite_pool_correction" => p.correction,
                    "mean_R50_mm" => p.meanR,
                    "mean_events" => p.nev, "n_fail" => p.nfail)
                    for r in results for p in r.iso]))
        end
        println("wrote $(joinpath(out, "sigma_r_per_isotope_v2$(tag).toml"))")
    end
    # The Grogg-estimator comparison: erfc R50 vs linear x-intercept, same
    # realizations (paired seeds), weighted primary + unweighted variant.
    open(joinpath(out, "sigma_r_grogg_v2$(tag).toml"), "w") do io
        TOML.print(io, Dict(
            "generation" => "v2", "scanner" => RING, "crystal" => CFG.crystal,
            "dose_Gy" => THIN_DOSE, "realizations" => REALIZATIONS,
            "sigma_band" => 1 / sqrt(2 * (REALIZATIONS - 1)),
            "method" => "Grogg linear x-intercept vs erfc R50, paired realizations",
            "point" => [Dict(
                "leaf" => r.leaf, "t_del_s" => r.t_del, "t_ac_s" => r.t_ac,
                "t2_s" => r.t2, "o15_fraction" => r.fO15,
                "o15_fraction_washed" => r.fO15_washed,
                "sigma_R_convention" => "finite-pool corrected",
                "nominal_finite_pool_correction" => r.correction_nominal,
                "washed_finite_pool_correction" => r.correction_washout,
                "erfc_nominal_sigma_R_raw_mm" => r.σ_nom_raw,
                "erfc_washed_sigma_R_raw_mm" => r.σ_wsh_raw,
                "erfc_nominal_sigma_R_mm" => r.σ_nom,
                "erfc_washed_sigma_R_mm" => r.σ_wsh,
                "erfc_nominal_R50_mean_mm" => r.R_nom,
                "grogg_nominal_sigma_R_mm" => r.grogg_w.σn,
                "grogg_washed_sigma_R_mm" => r.grogg_w.σw,
                "grogg_nominal_sigma_R_raw_mm" => r.grogg_w.σn_raw,
                "grogg_washed_sigma_R_raw_mm" => r.grogg_w.σw_raw,
                "grogg_nominal_Rx_mean_mm" => r.grogg_w.Rn,
                "grogg_washed_Rx_mean_mm" => r.grogg_w.Rw,
                "grogg_n_fail_nominal" => r.grogg_w.fn,
                "grogg_n_fail_washed" => r.grogg_w.fw,
                "grogg_unw_nominal_sigma_R_mm" => r.grogg_u.σn,
                "grogg_unw_washed_sigma_R_mm" => r.grogg_u.σw,
                "grogg_unw_nominal_sigma_R_raw_mm" => r.grogg_u.σn_raw,
                "grogg_unw_washed_sigma_R_raw_mm" => r.grogg_u.σw_raw,
                "grogg_unw_nominal_Rx_mean_mm" => r.grogg_u.Rn,
                "grogg_unw_washed_Rx_mean_mm" => r.grogg_u.Rw,
                "grogg_unw_n_fail_nominal" => r.grogg_u.fn,
                "grogg_unw_n_fail_washed" => r.grogg_u.fw,
                # Grogg's full pipeline: 7 mm FWHM smoothing before the fit
                "grogg_sm7_nominal_sigma_R_mm" => r.grogg_s.σn,
                "grogg_sm7_washed_sigma_R_mm" => r.grogg_s.σw,
                "grogg_sm7_nominal_sigma_R_raw_mm" => r.grogg_s.σn_raw,
                "grogg_sm7_washed_sigma_R_raw_mm" => r.grogg_s.σw_raw,
                "grogg_sm7_nominal_Rx_mean_mm" => r.grogg_s.Rn,
                "grogg_sm7_washed_Rx_mean_mm" => r.grogg_s.Rw,
                "grogg_sm7_n_fail_nominal" => r.grogg_s.fn,
                "grogg_sm7_n_fail_washed" => r.grogg_s.fw,
                # per-realization endpoints + chosen Grogg range, paired by
                # index — the raw material for paired/mechanism analysis
                "realizations_erfc_nominal_mm" => r.est[:nom],
                "realizations_erfc_washed_mm" => r.est[:wsh],
                "realizations_erfc_nominal_chi2_dof" => r.est[:nom_chi2],
                "realizations_erfc_washed_chi2_dof" => r.est[:wsh_chi2],
                "realizations_grogg_nominal_mm" => r.est[:nomg],
                "realizations_grogg_washed_mm" => r.est[:wshg],
                "realizations_grogg_unw_nominal_mm" => r.est[:nomgu],
                "realizations_grogg_unw_washed_mm" => r.est[:wshgu],
                "realizations_grogg_nominal_z_first_mm" => r.est[:nomg_zf],
                "realizations_grogg_nominal_z_last_mm" => r.est[:nomg_zl],
                "realizations_grogg_nominal_slope" => r.est[:nomg_sl],
                "realizations_grogg_washed_z_first_mm" => r.est[:wshg_zf],
                "realizations_grogg_washed_z_last_mm" => r.est[:wshg_zl],
                "realizations_grogg_washed_slope" => r.est[:wshg_sl],
                "realizations_grogg_sm7_nominal_mm" => r.est[:nomgs],
                "realizations_grogg_sm7_washed_mm" => r.est[:wshgs],
                "realizations_grogg_sm7_nominal_z_first_mm" => r.est[:nomgs_zf],
                "realizations_grogg_sm7_nominal_z_last_mm" => r.est[:nomgs_zl],
                "realizations_grogg_sm7_nominal_slope" => r.est[:nomgs_sl],
                "realizations_grogg_sm7_washed_z_first_mm" => r.est[:wshgs_zf],
                "realizations_grogg_sm7_washed_z_last_mm" => r.est[:wshgs_zl],
                "realizations_grogg_sm7_washed_slope" => r.est[:wshgs_sl])
                for r in results]))
    end
    println("wrote $(joinpath(out, "sigma_r_washout_v2$(tag).toml"))")
    println("wrote $(joinpath(out, "sigma_r_grogg_v2$(tag).toml"))")
end

main()
