# CsI-ring washout σ_R — N=200 firm-up (del120, del180)

Provenance note for `sigma_r_washout_v2_firmup_N200.toml` in this directory.
Read this before trusting or regenerating those numbers.

## TL;DR

The CsI reference-ring (TBP: `crysp_ring_1m_csi_2x0`, r_inner 387 mm, AFOV
1024 mm) **full-window** washed σ_R at N=100 was flaky. del180 (0.199) had
dipped *below* del120 (0.232), which is non-physical (a later, count-poorer
window must have σ_R ≥ the earlier one). Re-running del120 and del180 at
**N=200** fixed it. These two leaves are now quoted at N=200; **del300 is still
N=100**. The N=200 values live in the separate `_firmup_N200.toml` file, **not**
in the main `sigma_r_washout_v2.toml` (which is deliberately left at N=100 — see
"Why a separate file").

## The numbers (1 Gy, all events, whole-plane erfc R50)

Counting expectation for washout is inflation = 1/√survival ≈ 1/√0.48 ≈ **1.44**.

| leaf | nominal | washed | inflation | N | note |
|------|---------|--------|-----------|---|------|
| del120 | 0.113 → **0.118** | 0.232 → **0.226** | 2.05 → **1.91** | 100→200 | washed held; inflation genuinely high |
| del180 | 0.174 → **0.164** | 0.199 → **0.221** | 1.15 → **1.35** | 100→200 | washed rose; N=100 was a low fluctuation |
| del300 | 0.228 | 0.309 | 1.35 | 100 | not rerun (looked fine) |

Interpretation:

- **del180 was the flaky one.** Its N=100 washed (0.199) was a genuine ~2σ low
  excursion; at N=200 it rises to 0.221 and monotonicity is restored
  (del120 ≈ del180 within the ±5% N=200 bar, both < del300).
- **del120 is NOT an outlier.** Its washed value barely moved (0.232 → 0.226)
  and the per-realization distribution is clean (MAD-std ≈ sample-std, no single
  blown-up fit). Its washout inflation of ~1.9 — well above the counting 1.44 —
  is a **real, stable feature** of the ¹⁵O-richest window, not noise. Why it is
  that high (isotope-mix reweighting by the non-uniform g_i? something else?) is
  **not yet explained** and would be worth understanding.
- The earlier "CsI is flat vs AFOV while BGO rises" reading (see the short-scan
  figure `washed_shortscan_afov.png`) leaned partly on the inflated N=100 ring
  point; with N=200 the crystal-dependent trends are within errors. Do not build
  a physics story on them without more statistics. Both crystals lose the same
  count fraction as AFOV shrinks (verified on the shards: ~0.47 trues at CAFOV
  for both) and are counting-limited.

## Why a separate file (not edited into the main TOML)

The main `sigma_r_washout_v2.toml` (N=100, all three leaves) feeds the §8
figures of `latex/endpoint_precision.tex` via `tools/plot_sigma_r_v2.py` and
`plot_washed_v2_scanners.py`, which read the top-level `realizations` field for
the error band. Splicing mixed-N leaves into it would give those figures a wrong
band for del300 and silently perturb a published note. So the main file is left
untouched at N=100, and the firm-up is quoted from `_firmup_N200.toml`. The
paper table `tab:sigmaR_combined` in `latex/cbs.tex` uses the N=200 values for
del120/del180 and the N=100 value for del300.

## How each number is backed (important — not fully symmetric)

- **del180**: recomputed **exactly** from its retained per-realization dump
  `firmup_N200_del180_dump.toml` (the `sigma_r_grogg_v2` output of the N=200 run,
  preserved under that stable name). σ_R = std(realizations_erfc_{nominal,washed}_mm).
- **del120**: taken from the **recorded N=200 run summary** (nominal 0.118,
  washed 0.226). Its per-realization dump was **overwritten** by the del180 run
  (both write the default `sigma_r_grogg_v2.toml` path), so there is no retained
  dump for del120. The value is still reproducible because the driver uses fixed
  seeds (see below).

`tools/make_ringcsi_firmup.py` writes `_firmup_N200.toml` from these two
sources; the `source` field on each point records which.

## How to reproduce / finish cleanly

Config: `config/run_parameters_csi_v2.toml` active (scanner
`crysp_ring_1m_csi_2x0`, crystal `csi`). Fixed seeds ⇒ deterministic output.

```
cp config/run_parameters_csi_v2.toml config/run_parameters.toml
julia -t auto --project=. drivers/sigma_r_v2.jl \
    --realizations 200 --leaves del120s_ac300s_1Gy --isotopes none    # → 0.118 / 0.226
julia -t auto --project=. drivers/sigma_r_v2.jl \
    --realizations 200 --leaves del180s_ac300s_1Gy --isotopes none    # → 0.164 / 0.221
```

Each single-leaf run **overwrites** `sigma_r_washout_v2.toml` (and
`sigma_r_grogg_v2.toml`) with just that leaf — back up the main file first and
restore it after, or run all three leaves in one invocation.

**To make it fully uniform-N** (recommended if this matters for publication):
run all three leaves at N=200 in one pass (`--realizations 200` with the default
leaves), which regenerates a clean N=200 `sigma_r_washout_v2.toml`. Then the
separate firm-up file + this note become redundant and can be retired, and the
§8 figures in `endpoint_precision.tex` should be regenerated and its CsI-ring
numbers moved to N=200 as well (they are still N=100 as of this writing).

## Open items

- del300 CsI-ring still N=100 (0.228 / 0.309).
- del120 N=200 per-realization dump not retained (value reproducible by command).
- The ~1.9 washout inflation at del120 is real and unexplained.
- `endpoint_precision.tex` §8 still carries the old N=100 CsI-ring numbers.
