# thinning.jl — W6: downstream realizations by seeded Bernoulli thinning of
# the pooled master (dev/reference/data_generation_strategy.md §4). The three
# non-negotiable properties, each fixing a failure mode:
#
#   1. POOL across all shards — p = target/M_total over the union; thinning
#      per shard file is the "p = 1 per file" degenerate trap.
#   2. BERNOULLI, never exact-count — each event kept independently, so a
#      realization's count fluctuates as Binomial(M_total, p) ≈ Poisson(target),
#      the same fluctuation a real acquisition has; that fluctuation is part
#      of σ_R.
#   3. OWN seed namespace — the RNG is seeded by the downstream realization
#      index, never entangled with the upstream (master_seed, shard_index);
#      thinning operates purely on the produced LORs and reproduces
#      bit-for-bit.

"""
The thinning RNG seed base: `MersenneTwister(THINNING_SEED_BASE +
realization_index)`. Its own namespace, distinct from every upstream seed
(master_seed 1, transport/time 1234) and the sensitivity draw seeds (1, 2) —
realization k is fully determined by k alone.
"""
const THINNING_SEED_BASE = 1_000_000

"""
    dose_to_counts(dose_Gy, top_dose_Gy, M_total, n_shards) -> Int

The physically anchored target count for a realization at `dose_Gy`:
`round(dose/top_dose · M_total/n_shards)` — ONE acquisition's pooled count
scaled by the dose ratio (β⁺ activity is linear in dose).

This resolves the PLAN W6 confirmation thread: the recipe's yield model
(activation yield/Gy × dose × washout-window fraction × geometry sensitivity)
collapses to the dose ratio, because each shard IS one simulated acquisition
at the top dose — its count already carries the yield, the window and the
sensitivity. The dose ratio scales one acquisition's statistics
(`M_total/n_shards`), and the keep-probability over the pooled union is
`p = target/M_total = (dose/top_dose)/n_shards` — ≤ 1/n_shards at the top
dose, the strategy's "p ≤ 0.1" rule. Thin down, never up: `dose_Gy` must not
exceed `top_dose_Gy`.
"""
function dose_to_counts(dose_Gy::Real, top_dose_Gy::Real, M_total::Integer,
                        n_shards::Integer)
    0 < dose_Gy <= top_dose_Gy ||
        throw(ArgumentError("dose $dose_Gy Gy outside (0, $top_dose_Gy] — thin down, never up"))
    return round(Int, dose_Gy / top_dose_Gy * M_total / n_shards)
end

"""
    thin_mask(M_total, target_counts, realization_index) -> BitVector

The seeded Bernoulli stream filter: keep each of the `M_total` pooled events
independently with `p = target_counts/M_total`, RNG =
`MersenneTwister(THINNING_SEED_BASE + realization_index)`. The kept count
fluctuates as Binomial(M_total, p); the mask reproduces bit-for-bit for the
same `(M_total, target_counts, realization_index)`.
"""
function thin_mask(M_total::Integer, target_counts::Integer,
                   realization_index::Integer)
    0 <= target_counts <= M_total ||
        throw(ArgumentError("target_counts $target_counts outside [0, $M_total]"))
    p = target_counts / M_total
    rng = MersenneTwister(THINNING_SEED_BASE + realization_index)
    mask = BitVector(undef, M_total)
    @inbounds for i in 1:M_total
        mask[i] = rand(rng) < p
    end
    return mask
end

"""
    thin_lm(coinc::MCCoincidences, target_counts, realization_index) -> BitVector

One realization over a pooled event list, as a mask aligned with `coinc` —
compose with the class masks and column-select:
`keep = thin_lm(c, n, k) .& is_true(c); xs, xe = endpoints(c, keep)`.
"""
thin_lm(coinc::MCCoincidences, target_counts::Integer, realization_index::Integer) =
    thin_mask(length(coinc), target_counts, realization_index)

"""
    finite_pool_correction(probabilities) -> Float64

Count-level correction from the conditional variance produced by repeated
Bernoulli thinning of a fixed event pool to the variance of an independently
generated Poisson acquisition.

For event keep probabilities `q_e`,

    Poisson variance      = sum(q_e)
    conditional variance  = sum(q_e * (1 - q_e))

and the standard-deviation correction is

    sqrt(sum(q_e) / sum(q_e * (1 - q_e))).

For a uniform probability `q`, this reduces to `1/sqrt(1-q)`.
"""
function finite_pool_correction(probabilities)
    q = collect(Float64, probabilities)
    isempty(q) && throw(ArgumentError("probabilities must not be empty"))
    all(isfinite, q) ||
        throw(ArgumentError("probabilities must be finite"))
    all(value -> 0.0 <= value <= 1.0, q) ||
        throw(ArgumentError("probabilities must lie in [0,1]"))

    poisson_variance = sum(q)
    conditional_variance = sum(value * (1.0 - value) for value in q)

    poisson_variance > 0 ||
        throw(ArgumentError("at least one probability must be positive"))
    conditional_variance > 0 ||
        throw(ArgumentError(
            "conditional variance is zero; the pool cannot estimate a spread"
        ))

    return sqrt(poisson_variance / conditional_variance)
end

finite_pool_correction(q::Real) = finite_pool_correction((q,))
