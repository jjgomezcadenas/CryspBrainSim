#!/bin/bash
# collect_note_figures.sh — copy the figures latex/endpoint_precision.tex
# includes from the (git-ignored) output tree into latex/figs/, so the note
# compiles from a fresh clone. Regenerate the sources first when needed:
#   tools/fit_activity_profile.py --no-pulls          (fit figure)
#   tools/plot_recon_projections.py                   (projections)
#   tools/ten_shards.py                               (ladder)
#   tools/ten_shards.py --dose-sweep                  (dose sweep)
#   tools/recon_scatters.jl + tools/scatter_profile.py (scatter profile)
set -euo pipefail
cd "$(dirname "$0")/.."

CFG=out/uniform_headep_sobp_1e8/closed/crysp_ring_1m/bgo_3X0
FIGS=latex/figs

cp "$CFG/one_shard/fits/figures/recon_activity.png"     "$FIGS/fit_recon_activity.png"
cp "$CFG/one_shard/fits/figures/recon_projections.png"  "$FIGS/recon_projections.png"
cp "$CFG/one_shard/fits/figures/scatters_activity.png"  "$FIGS/scatters_profile.png"
cp "$CFG/ten_shards/figures/delta_r50.png"              "$FIGS/ladder_delta_r50.png"
cp "$CFG/ten_shards/figures/delta_r50_vs_dose.png"      "$FIGS/dose_sweep_r50.png"
echo "collected 5 figures into $FIGS/"
