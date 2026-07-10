# Isotope washout — problem, upstream responses, and resolution

**Resolution up front:** the washout-as-loss study is **fully downstream and needs nothing from
upstream** — not a new column, not a clock revision, not a master regeneration. This note records
the problem, the two upstream responses (they disagree), the assessment, and exactly how the whole
loss study is carried on the frozen source we already have.

## The physical problem

The β⁺ emitters produced along the beam path (¹⁵O 122 s, ¹¹C 1223 s, ¹³N 598 s, ¹⁰C 19 s, ¹⁴O 71 s)
do not sit where they were created: perfusion and metabolism clear the atoms before they decay —
the Mizuno three-exponential survival W(t) (fast/medium/slow) in `latex/washout_brain.tex`. Each
produced nucleus carries two competing clocks, physical decay τ_d and clearance τ_c, and the fate of
a nucleus with τ_c < τ_d splits the problem in two:

- **Loss.** The cleared atom leaves the field of view; its annihilation is never recorded near the
  beam. Washout is a pure removal process on the local activity.
- **Redistribution.** The cleared atom is transported by perfusion and may still decay **inside** the
  FOV, contributing a **displaced** annihilation — removal plus transport.

The observable is the beam-corridor depth profile's distal-edge midpoint R₅₀ and its run-to-run
spread σ_R. Loss removes counts (σ_R up) and reshapes the isotope mix (R₅₀ may move); redistribution
can additionally bias the edge shape.

## Upstream responses (two voices, they disagree)

### G4 (the simulation expert) — the correct analysis

- **The equivalence.** Washout-as-loss is a per-nucleus independent survival, so imposing it inside
  Geant4 (keep if τ_c > τ_d, probability W(τ_d)) and imposing it downstream (Bernoulli-thin by
  W(Δt), Δt = age at decay) select the **same in-window set** with the same mix shift, σ_R, and
  normalization. A removal process cannot move a decay time, so the expensive in-MC route adds
  nothing.
- **No time column is needed** (the key point). The per-event age integrates out analytically:
  the window-integrated survival for isotope j at position x,

    g_j(x) = n_j(x) / m_j(x),   n_j = ∫∫ ρ_j·s(t′)·λ_j e^{−λ_jτ}·W(τ)·[t′+τ ∈ window] dτ dt′,

  with m_j the same integral without W, is a **closed-form scalar** over the known creation shape
  s(t′) and decay law. Because production is Poisson and clearance is independent per nucleus, the
  true washout-recorded process and the constant-g_j-thinned process are **both Poisson with the
  identical rate n_j** — so g_j reproduces the correct noise (σ_R), not just the mean, with **no
  per-event timing**. If W is spatially uniform, g_j is one scalar per isotope and the whole
  depth/R₅₀ effect comes through the isotope mix varying with depth.
- **Redistribution is not Geant4 either.** It is transport of the *nucleus* (perfusion), a
  downstream geometry/physiology model on the production points — not proton transport. Stage A
  stays frozen in both regimes.
- **Redistribution deferred.** Matches the forward-MC literature; physically the acquisition delay
  suppresses the redistributing (fast, diffusible) component. If ever revived it is compartmental
  kinetics on time-activity curves — a different tool, still not Geant4.

### PTCrysp (the pipeline) — the offer, and where it reverts

The follow-up offer proposed applying loss "by thinning on **t_prod_s** / isotope / position,"
supplying per-event **t′ + isotope columns** and a **clock revision**, and asked which decisions to
settle. This **contradicts G4's own analysis one step earlier**: it re-introduces exactly the
per-event time column and clock question that the analytic g_j retired, and it would break the
pipeline's stated invariant that `emitters.csv` is timeless.

**The data model settles the tie in G4's favor.** `sobp_layers.csv` carries **energy + weight only —
no delivery-time or painting-order column**. The production in this scenario is time-integrated with
no depth–time structure. The only thing that could justify a per-event t_prod over the analytic
scalar is a **depth–time coupling** (a layer-painting sequence correlating creation time with
depth), and it **does not exist in the current model**. A per-event t_prod would be *synthesized*
from the assumed uniform R(t′) — carrying zero information beyond the uniform model, hence strictly
equivalent to g_j and more fragile. Emitting it would cost the invariant and buy nothing.

## Assessment

G4 is right; PTCrysp's offer is a regression to a superseded design. Of PTCrysp's three "open
decisions," two are artifacts of that regression:

- **Production profile for t′** — uniform over the 60 s irradiation is not merely a default; it is
  the only structure the model contains, and under it the analytic g_j needs **no column**.
- **Clock / delay convention** — **moot** under the analytic route: there is no per-event clock. The
  delay-era losses are already inside g_j, which integrates survival over the window [120, 1320] s
  with creation in [0, 60] s.
- **Regeneration scope** — would be the only real question *if* a regeneration were needed; it is
  not (see below).

## What we do with nothing from upstream

### Truth level — runnable now

Reweight the five per-isotope depth columns of `truth/activity_profile_fast.csv` (`O15,C11,N13,C10,
O14`) by their g_j, re-sum to a washed total, refit Δ_R50, compare to nominal. g_j is computed here
from the Mizuno W (we own it), the acquisition timing we already own, and the uniform creation
model. **Zero upstream.** This is the first-order answer to "does washout move the edge."

### Detected level — from the transfer functions we already measured

The reconstructed signal is **isotope-blind**: an annihilation gives a 511 keV pair regardless of
parent isotope, so attenuation, acceptance, and MLEM act on the annihilation distribution, not on
which isotope produced it (the small per-isotope positron-range difference is already folded into
the stored annihilation positions). Writing the detector transfer T as linear and isotope-blind, the
washed detected image is T(Σ_j g_j A_j) = Σ_j g_j T(A_j). Two consequences let us reach the detected
observables without tagging any event by isotope:

- **Bias.** Δ_R50(detected, washout) = Δ_R50(detected, nominal — already measured) + [Δ_R50(truth
  washed) − Δ_R50(truth nominal)]. The bracket is the truth-level shift above; the transfer holds
  because the truth→detected offset was measured to be a stable, mix-independent constant
  (≈ −0.4 mm, shard-independent, in the endpoint study).
- **Precision.** σ_R(washout, 1 Gy) = σ_R(nominal, 1 Gy) / √f, with the overall kept fraction
  f = Σ_j g_j N_j / Σ_j N_j from the truth abundances N_j (the column integrals) and g_j, applied
  through the **exact 1/√dose law** the dose sweep already established. Washout-as-loss is, for
  precision, an isotope-weighted dose reduction.

Both detected-level numbers come from truth-level g_j composed with detector transfer functions the
endpoint study already delivered. **Zero upstream.**

### The one second-order term — bounded here, not asked of upstream

The only thing a per-LOR isotope tag would add is the isotope-*differential* re-thinning of the
actual detected events — the change in edge **shape** beyond the linear mix-shift, i.e. T's blur
interacting with the mix change and the offset's weak dependence on it. This is second order. We
**bound it downstream** (e.g. from the spread of per-isotope edge positions and widths in the truth
profiles) rather than request a column. Per the standing guidance — do not ask upstream for
verifications or cosmetics — we leave the pipeline untouched.

## Conclusion

- **The loss study needs nothing from upstream.** Truth level: reweight the isotope columns by an
  analytic g_j. Detected level: compose that mix-shift with the already-measured truth→detected
  offset (bias) and 1/√dose law (precision). No column, no clock change, no regeneration, no G4.
- **Redistribution stays deferred**; if revived it is a downstream perfusion model on the production
  points, still not Geant4.
- **Send upstream:** we accept the redistribution deferral and adopt loss as downstream, and we do
  **not** need the t′/isotope columns or the clock revision — G4's analytic result covers it. Leave
  the pipeline as is.
- **Next action here:** run the truth-level g_j reweight; that first-order number gates whether the
  second-order shape term is even worth bounding.

## Related

- Loss-route scoping and the isotope-mix framing: [`isotope-washout.md`](isotope-washout.md).
- Washout model and parameters: `latex/washout_brain.tex`.
