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

# The analysis surface is filled in over dev/PLAN.md build steps 2–6:
#   endpoint.jl · profile.jl · products.jl · qa.jl · characterize.jl · mumap.jl
#   sensitivity.jl · thinning.jl · dualhead_sampler.jl

end # module CryspBrainSim
