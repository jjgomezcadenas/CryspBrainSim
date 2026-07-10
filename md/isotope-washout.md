# Next problem: isotope washout (IW)

Every production to date models **physical decay only** — the β⁺ emitters sit where they were
created until they decay. In tissue they do not: dissolved gases and metabolites (¹⁵O, ¹¹C, ¹³N)
are cleared by perfusion/diffusion before decaying. This is **isotope washout (IW)** — the missing
physics between our simulation and a real patient.

Literature model: a sum of exponentials (fast/medium/slow biological components), per isotope and
chemical form, folding with the physical half-life into an effective one.

## The two questions (mirroring the endpoint study)

- **(a) Does IW degrade precision?** It removes counts, worst at late times, compounding the
  delayed-start loss already measured.
- **(b) Does IW bias the edge?** Only if clearance is depth-/tissue-dependent. In the uniform-brain
  phantom it is roughly uniform, so first-order the edge should hold; the test is whether Δ_R50
  moves under IW.

## Where IW enters — the scoping decision for the next session

- **(i) Truth level — computable now, the natural first step.**
  `truth/activity_profile_fast.csv` has per-isotope depth profiles (columns
  `O15,C11,N13,C10,O14,total`). **Spatially-uniform IW is a per-isotope multiplicative reweighting**
  of those columns — each isotope's P(z) scaled by its window-integrated survival fraction under its
  clearance model — then refit for Δ_R50. Needs no new detected-level data. Tests question (a) at
  truth level and, since the reweighting shifts the isotope mix, whether the edge moves (b).
- **(ii) Detected level — needs more data.**
  - *Upstream* — a new production with IW kinetics applied (cleanest, but needs a PTCryspMC feature
    + run).
  - *Downstream reweighting* on the stored runs — an event decaying at `t_decay_s` had to survive
    clearance from its creation to its decay, so it could be Bernoulli-thinned by a survival
    probability. That needs the per-event **creation time and isotope** — the isotope id was
    explicitly *dropped* from the decay-time request, and creation time is not stored either — so
    downstream is **not currently possible** without more columns.

Decide the axis (start with (i)), then whether an upstream request (or per-event isotope +
creation-time columns) is warranted.

## Related

- **G4 vs downstream — resolved:** [`washout-g4-formulation.md`](washout-g4-formulation.md).
  The loss study is **fully downstream, zero upstream**: an analytic per-isotope survival scalar
  g_j (from the Mizuno W we own) reweights the truth isotope columns; the detected-level bias and
  σ_R follow from the already-measured truth→detected offset and 1/√dose law, because detection is
  isotope-blind. No column, no clock change, no regeneration, no Geant4. Redistribution deferred
  (and if revived, a downstream perfusion model, still not Geant4).
- The **composite-erfc edge model** (per-isotope components) is the same object as the IW
  reweighting — the per-isotope decomposition. Worth building together.
- **Washout is spatially uniform only to first order** — a real perfusion field is not; if the
  truth-level study shows the edge is sensitive to the mix, spatially-dependent clearance becomes
  the reason to go upstream.
