#!/usr/bin/env julia
# Restartable GPU/CPU jobs for the statistical-procedure reference case.
# Each expensive reconstruction writes one immutable TOML result immediately.
#
# Examples (run from a normal Terminal for Metal access):
#   julia -t 1 --project=. drivers/statistical_procedure_jobs.jl --stage shard --index 0
#   julia -t 1 --project=. drivers/statistical_procedure_jobs.jl --stage ensemble --mode washed --first 1 --last 10
#   julia --project=. drivers/statistical_procedure_jobs.jl --stage combine --mode washed

using CryspBrainSim
using RecoCryspTools
using Metal
using Printf
using Random: MersenneTwister
using Statistics: mean, std
using TOML

argvalue(flag, default, parsefn=identity) = let i = findfirst(==(flag), ARGS)
    i === nothing ? default : parsefn(ARGS[i + 1])
end
hasflag(flag) = any(==(flag), ARGS)

const REPO = dirname(@__DIR__)
const PRODUCTS_ROOT = joinpath(dirname(REPO), "PtCryspProds")
const CONFIG_PATH = abspath(argvalue("--config", joinpath(REPO, "config", "run_parameters_ring_bgo_v2.toml"), String))
const PARAMS = load_run_parameters(CONFIG_PATH)
const CFG = PARAMS.config
const LEAF = argvalue("--leaf", "del120s_ac300s_1Gy", String)
const DOSE_GY = argvalue("--dose", 1.0, x -> parse(Float64, x))
const STAGE = argvalue("--stage", "", String)
const MODE = argvalue("--mode", "washed", String)
const INDEX = argvalue("--index", 0, x -> parse(Int, x))
const FIRST = argvalue("--first", 1, x -> parse(Int, x))
const LAST = argvalue("--last", FIRST, x -> parse(Int, x))
const SEED_BASE = argvalue("--seed", 9_000_000, x -> parse(Int, x))
const DEV = Metal.functional() ? MtlArray : identity

function context()
    cache = joinpath(sensitivity_out(CFG.scenario, CFG.topology, CFG.scanner), sensitivity_cache_name(PARAMS))
    load_run_context(products_root=PRODUCTS_ROOT, scenario=CFG.scenario, topology=CFG.topology,
                     scanner=CFG.scanner, crystal=CFG.crystal, leaf=LEAF,
                     sens_cache=cache, params=PARAMS)
end

case_dir(ctx) = joinpath(config_out(ctx.scenario, ctx.topology, ctx.ring, ctx.crystal),
                          "statistical_procedure", "$(LEAF)_D$(replace(string(DOSE_GY), "." => "p"))Gy")

function write_toml(path, data)
    mkpath(dirname(path))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        TOML.print(io, data)
    end
    mv(temporary, path; force=true)
end

function metadata(ctx)
    Dict("schema" => "statistical-procedure-v1", "config_path" => CONFIG_PATH,
         "scenario" => ctx.scenario, "topology" => ctx.topology, "scanner" => ctx.ring,
         "crystal" => ctx.crystal, "leaf" => LEAF, "dose_Gy" => DOSE_GY,
         "niter" => PARAMS.niter, "backend" => Metal.functional() ? "metal" : "cpu")
end

function shard_job()
    ctx = context()
    0 <= INDEX < length(ctx.files) || error("--index must be in 0:$(length(ctx.files)-1)")
    path = joinpath(case_dir(ctx), "shards", @sprintf("shard%03d.toml", INDEX))
    isfile(path) && (println("exists: $path"); return)
    shard = read_shard(ctx.files[INDEX + 1])
    shard.n_dropped == 0 || error("shard $INDEX dropped $(shard.n_dropped) LORs")
    xs, xe = endpoints(shard.coinc)
    result = reconstruct_endpoint(ctx, xs, xe, lor_attenuation(ctx, xs, xe);
                                  device=DEV, return_profile=true)
    data = metadata(ctx)
    merge!(data, Dict("stage" => "shard", "status" => "complete", "index" => INDEX,
                      "input_file" => ctx.files[INDEX + 1], "events" => result.nev,
                      "R50_mm" => result.r50_fit, "R50_fit_error_mm" => result.z0_err,
                      "erfc_chi2_dof" => result.erfc_chi2_dof, "erfc_popt" => result.erfc_popt,
                      "fit_window_mm" => collect(result.fit_window_mm),
                      "profile_z_mm" => result.profile_z_mm, "profile" => result.profile))
    write_toml(path, data)
    run(`$(Sys.which("python")) tools/plot_statistical_procedure.py --shard $path`)
    println("wrote $path")
end

function ensemble_job()
    MODE in ("nominal", "washed") || error("--mode must be nominal or washed")
    FIRST >= 1 && LAST >= FIRST || error("invalid --first/--last range")
    ctx = context(); root = joinpath(case_dir(ctx), MODE)
    pool = pool_shards(ctx.files); pool.n_dropped == 0 || error("pooled sample dropped LORs")
    isotope = vcat([shard_isotope(file) for file in ctx.files]...)
    survival = Float64.(shard_attrs(ctx.files[1])["washout_g"])
    p_dose = DOSE_GY / length(ctx.files)
    q = MODE == "nominal" ? fill(p_dose, length(isotope)) :
        p_dose .* [survival[Int(id)+1] for id in isotope]
    expected, conditional = sum(q), sum(q .* (1 .- q))
    xs, xe = endpoints(pool.coinc); attenuation = lor_attenuation(ctx, xs, xe)
    for index in FIRST:LAST
        path = joinpath(root, @sprintf("realization%04d.toml", index))
        isfile(path) && (println("exists: $path"); continue)
        keep = rand(MersenneTwister(SEED_BASE + index), length(q)) .< q
        result = reconstruct_endpoint(ctx, xs[:, keep], xe[:, keep], attenuation[keep]; device=DEV)
        data = metadata(ctx)
        merge!(data, Dict("stage" => "ensemble", "mode" => MODE, "status" => "complete",
                          "index" => index, "seed" => SEED_BASE + index, "events" => result.nev,
                          "R50_mm" => result.r50_fit, "R50_fit_error_mm" => result.z0_err,
                          "erfc_chi2_dof" => result.erfc_chi2_dof,
                          "expected_events" => expected, "conditional_count_variance" => conditional))
        write_toml(path, data)
        println("wrote $path")
    end
end

function combine_job()
    ctx = context(); root = joinpath(case_dir(ctx), MODE)
    files = sort(filter(path -> endswith(path, ".toml"), readdir(root; join=true)))
    isempty(files) && error("no completed $MODE realisations in $root")
    rows = TOML.parsefile.(files)
    all(row -> row["status"] == "complete" && row["mode"] == MODE, rows) || error("incompatible result files")
    indices = [Int(row["index"]) for row in rows]
    length(unique(indices)) == length(indices) || error("duplicate realization indices")
    endpoints = Float64[row["R50_mm"] for row in rows]
    all(isfinite, endpoints) || error("non-finite endpoint in completed results")
    expected = Float64(rows[1]["expected_events"]); conditional = Float64(rows[1]["conditional_count_variance"])
    correction = sqrt(expected / conditional)
    output = Dict("stage" => "combined", "mode" => MODE, "realizations" => length(rows),
                  "indices" => indices, "R50_mm" => endpoints, "raw_sigma_R_mm" => std(endpoints),
                  "finite_pool_correction" => correction,
                  "corrected_sigma_R_mm" => correction * std(endpoints),
                  "mean_R50_mm" => mean(endpoints), "n_fail" => 0,
                  "relative_sigma_uncertainty" => 1 / sqrt(2 * (length(rows)-1)))
    merge!(output, metadata(ctx))
    path = joinpath(case_dir(ctx), "combined", "$(MODE)_N$(length(rows)).toml")
    write_toml(path, output); println("wrote $path")
end

STAGE == "shard" ? shard_job() : STAGE == "ensemble" ? ensemble_job() :
STAGE == "combine" ? combine_job() : error("--stage must be shard, ensemble, or combine")
