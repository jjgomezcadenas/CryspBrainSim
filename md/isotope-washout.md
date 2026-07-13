# Isotope washout (IW) — loss study DONE (2026-07-10)

Productions model **physical decay only**; in tissue perfusion/metabolism clears the β⁺ emitters
(¹⁵O, ¹¹C, ¹³N…) before they decay. Washout **as loss** is now measured, fully downstream, both
1 m arms. Numbers + method: [`md/results.md`](results.md) ("Isotope washout"); the note carries it
as §7 of `latex/endpoint_precision.tex`; the G4-vs-downstream exchange in
[`washout-g4-formulation.md`](washout-g4-formulation.md).

**`latex/washout_brain.tex` is the complete, self-contained methodology spec — written to send
upstream to implement washout in the production chain.** It carries: the Mizuno brain model and the
window-integrated per-isotope survival g_i (closed form); the Poisson-thinning exactness proof; the
mix reweighting; and the **implementation recipe** — the recommended upstream form is a single
**per-decay Bernoulli keep with probability g_i(species)** (exact for uniform clearance, no per-event
age, no re-simulation), with the age-resolved W(τ) variant, the label-free marginalised weight
w(z₀,t_decay) = Σ_i P(i|z₀,t_decay) g_i we use downstream (posterior derived), the per-species
posterior selector for σ_R^(i), and the spatial-non-uniformity route to the one genuine bias.

**Realized in generation-2 (2026-07-13):** upstream now stamps the per-LOR isotope label and the
per-window `washout_g`, so `drivers/sigma_r_v2.jl` applies the **exact** per-species Bernoulli keep
`p_dose·g_i[isotope]` (the recommended path above) — the marginalised w is demoted to the label-free
cross-check. The loss study was re-run this way across six v2 scanners (both crystals, three bores):
same conclusions — **no bias** (ΔR₅₀ ±0.08 mm), **~1.5× σ_R cost** tracking counts not geometry —
and the exact per-isotope selection **definitively refuted** the positron-range σ_R floor (¹⁵O more
precise per count than ¹¹C; BGO gives the first clean ¹¹C point). Numbers: [`results.md`](results.md)
"Generation-2 σ_R study"; note §8.

## The result (uniform Mizuno brain clearance)

The two questions, answered:

- **(a) Precision — yes, a real cost.** Washout removes ~57% of counts **nearly uniformly** (all
  per-isotope survivals g_i ≈ 0.4–0.5), so it is the ordinary counting penalty **σ_R × ~1.5×**
  (0.11 → ~0.16 mm at 1 Gy), roughly flat vs start time and scanner, rising to ~2× only at the
  count-starved 300 s. **Not** the ¹⁵O variance-drain that protects the delayed *start* (those were
  conflated in the first pass; the "washout free at t=0" reading was an n=10 outlier artifact,
  corrected).
- **(b) Bias — none, for uniform clearance.** The edge shift (ΔR₅₀^wo ≈ +0.22 mm at t=0, shrinking
  with delay) is a deterministic constant the per-scanner reference calibrates away; its
  washout-**parameter** systematic is ±0.02 mm ≪ σ_R.

## How it was done (the key methods)

- **Truth level** (`tools/washout.py`): the per-isotope survival g_i is a closed-form scalar (note
  Eq. 7), cross-checked vs direct integration; reweight the five truth columns, refit Δ_R50; MC over
  the Mizuno uncertainties gives the ±0.02 mm band.
- **Detected level** (`drivers/washout_sigma_r.jl --thinned`): per-event thinning of the pooled
  master by w(z₀,t_decay) = Σ_i P(i|z₀,t_decay) g_i — the isotope-**marginalised** survival, built
  from the truth profiles and decay laws. Uses only (z₀, t_decay), **already in the shards** — the
  earlier belief that detected level "needs creation-time + isotope columns" was wrong; detection
  is isotope-blind, so no new columns, no upstream, no Geant4. σ_R measured with the thinned method
  (dose-adaptive 0.2–1 Gy to keep the count-starved washed corner in the stable-fit regime; 0.1 Gy
  was verified to fail).

## Open — spatial non-uniformity (the one genuine bias route)

Uniform clearance is first-order. A real perfusion field is heterogeneous; a **depth-dependent**
clearance would make the edge shift depth-dependent — a genuine bias no single calibration constant
absorbs. Capturing it needs a **downstream perfusion/compartment transport model** on the production
points (not Geant4, not the range estimator). Not pursued here; the uniform result bounds the loss
effect. Redistribution (cleared atoms decaying elsewhere in the FOV) is the related deferred item —
also downstream perfusion physics, not Geant4 (see washout-g4-formulation.md).

## Related

- **Composite-erfc edge model** (per-isotope components) = the same per-isotope decomposition as the
  IW reweighting; still worth building (listed in [`pending.md`](pending.md)).
- Per-arm ready configs `config/run_parameters_{bgo,csi}.toml` (`cp` to activate) were added this
  round; the thinned σ_R machinery + failure guard in `washout_sigma_r.jl` are reusable.
