"""
    CryspBrainSim

Reconstruction and range analysis for the proton-therapy PET study: read the list-mode LOR shards
a PET scanner records for a proton SOBP field, reconstruct the β⁺ activity, extract the distal
range endpoint R, and report its statistical precision σ_R vs dose — one curve per scanner geometry.

The build plan is `dev/PLAN.md`; the upstream contracts are vendored in `dev/reference/`. The
reconstruction engine comes from `RecoCryspTools` (re-exports the RecoCrysp core).
"""
module CryspBrainSim

using RecoCryspTools
using LsqFit: curve_fit, coef, vcov, PrecisionWeights
using SpecialFunctions: erfc, erfcinv
using Statistics: mean, median, std
using LinearAlgebra: diag
using HDF5: h5open, read_attribute, attributes
using JSON3
using TOML
using NPZ: npzread, npzwrite
using Random: MersenneTwister

export depth_profile, distal_window, fit_endpoint, sigma_R
export read_csv_table, leaf_dir, shard_files, shard_attrs, read_shard,
       shard_t_decay, shard_isotope, shard_generation, pool_shards,
       scanner_geometry, phantom_region,
       material_mu, truth_dir
export ShardStats, shard_stats
export TruthReference, characterize, distal_crossing, windowed_crossing,
       read_depth_dose, read_activity_profile, write_reference
export phantom_attenuation, attenuation_ellipsoid, centered_grid, build_mumap,
       attenuation_mumap
export sensitivity_base, scaled_sensitivity, save_sensitivity,
       load_sensitivity, recocrysp_sha, sensitivity_cache_name
export load_run_parameters
export THINNING_SEED_BASE, dose_to_counts, thin_mask, thin_lm
export dose_tag, load_run_context, lor_attenuation, reconstruct_endpoint
export crystal_label, CRYSTAL_X0_MM, out_root, scenario_out, truth_out,
       mumap_out, ring_out, sensitivity_out, config_out, validation_out
export scanner_spec, write_ring_geometry, write_crystal_spec, write_descriptors

include("config.jl")
include("output.jl")
include("profile.jl")
include("endpoint.jl")
include("products.jl")
include("shard_stats.jl")
include("characterize.jl")
include("mumap.jl")
include("sensitivity.jl")
include("thinning.jl")
include("reconstruct.jl")

# The remaining analysis surface arrives over dev/PLAN.md build steps 5–6:
#   thinning.jl · dualhead_sampler.jl

end # module CryspBrainSim
