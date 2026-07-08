# Request to PTCryspMC — add the absolute decay time to the coincidence files

**From:** the CryspBrainSim analysis repo (reconstruction + range endpoint).
**Goal:** enable the acquisition-start-time study downstream — emulate any delayed acquisition
window as a pure event selection on the existing `fast` master, with no new productions and no
reweighting.

## Why

The endpoint study quantifies, rung by rung, what degrades the β⁺ range prediction: the
reconstruction, the event mix, and the **acquisition timing** — a later start loses the short-lived
isotopes (C¹⁰ T½ 19 s, O¹⁴ 71 s, O¹⁵ 122 s) and reshapes the distal edge. The truth side of that
study is computable today from the per-isotope columns of the `truth/` bundle. The detected side is
not: the coincidence files carry only photon times **relative to the decay** (`t1_ns`, `t2_ns`),
so a decay's position inside the acquisition window is unrecoverable downstream, and each timing
scenario would need its own production.

One column closes this. With the absolute decay time per coincidence, a delayed start is the cut
`t_decay ≥ t_start` on the stored list:

- **statistically exact** — the surviving events are precisely the acquisition that a scanner
  starting at `t_start` records: Poisson counts, the correct isotope mix (each isotope's decay
  times were already drawn from its own decay law via `time_seed`), no analytic decay factors;
- **randoms come out right for free** — randoms follow the instantaneous activity², and the stored
  randoms carry the real time structure, so the time cut selects exactly the randoms of the
  sub-window;
- **composes with the thinning** — the downstream dose axis (Bernoulli `thin_lm`) and the timing
  axis (this cut) apply independently to the same pooled master.

The quantity exists upstream already — the per-isotope decay times are drawn when the source is
materialized — and is dropped when the coincidences are written. This is plumbing, not physics.

We considered also requesting a per-event isotope id and dropped it: with the decay time in hand
it adds nothing actionable (the isotope is unobservable in a real measurement, and the per-isotope
truth profiles already answer the mix questions at truth level).

## The ask

One new column in `lors_shardNNN.h5` (all three modes: `lors_truth`, `lors_det`, `randoms`),
appended additively to the existing schema:

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `t_decay_s` | Float32 | s | absolute decay time of the annihilation within the acquisition window (gamma 1's decay for randoms, matching the `x0` convention) |

Conventions:

- **Zero point:** the acquisition start of the leaf's budget (`fast`: the window `[0, t_meas_s]`),
  consistent with the existing `t_meas_s` / `t0_s` clock attrs.
- **Randoms:** the two gammas come from different decays; store gamma 1's decay time — the same
  convention the schema already uses for `x0/y0/z0`.
- **Precision:** Float32 seconds resolves ~0.1 ms at 1200 s, ample for window cuts (the sub-ns
  physics stays in `t1_ns`/`t2_ns`).
- **Attrs:** stamp `t_decay_zero = "acquisition_start"` (or equivalent) so the file remains
  self-describing; regenerate `docs/SCHEMA.md` as usual.

Cost: 4 bytes × ~17.4 M rows ≈ 70 MB per shard (~14% of the current 515 MB).

## Backfill

Backfill the ten existing BGO shards of `uniform_headep_sobp_1e8/crysp_ring_1m/bgo/fast_1Gy/`
(and any future arm by default). The decay-time draw is seeded (`time_seed`), so regeneration
reproduces the same events with the column attached.

## Scope note

The cut emulates any window **inside** the produced one: a delayed start up to `t_meas_s`, with
acquisition ending at the produced window's end. Windows extending beyond it (e.g. a true
`offline` budget) remain separate productions — this request does not replace the budget axis,
it makes the start-time direction of it free.

## What we do with it (first analyses)

- Δ_R50(t_start) and Δ_R_p(t_start): the activity−dose endpoint differences as the acquisition
  start is delayed from 0 to ~600 s, at full statistics over the ten shards — the timing rung of
  the endpoint-degradation ladder.
- The same cut composed with `thin_lm` gives (dose × start-time) maps when part (b) needs them.
