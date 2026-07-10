#!/bin/bash
# collect_note_figures.sh — copy the figures latex/endpoint_precision.tex
# includes from the (git-ignored) output tree into latex/figs/, so the note
# compiles from a fresh clone. Both scanner arms are collected regardless of
# the active [configuration]. Regenerate the sources first when needed (per
# arm, with that arm active):
#   tools/fit_activity_profile.py --no-pulls          (fit figure)
#   tools/plot_recon_projections.py                   (projections)
#   tools/ten_shards.py  /  --dose-sweep              (ladder, dose sweep)
#   tools/recon_scatters.jl + tools/scatter_profile.py (scatter profile)
set -euo pipefail
cd "$(dirname "$0")/.."

FIGS=latex/figs
BGO=out/uniform_headep_sobp_1e8/closed/crysp_ring_1m_bgo_2x0/bgo_195k_2X0
CSI=out/uniform_headep_sobp_1e8/closed/crysp_ring_1m_csi_2x0/csi_2X0

cp "$BGO/one_shard/fits/figures/recon_activity.png"     "$FIGS/fit_recon_activity_bgo.png"
cp "$BGO/one_shard/fits/figures/recon_projections.png"  "$FIGS/recon_projections_bgo.png"
cp "$BGO/one_shard/fits/figures/scatters_activity.png"  "$FIGS/scatters_profile_bgo.png"
cp "$BGO/ten_shards/figures/delta_r50.png"              "$FIGS/ladder_delta_r50_bgo.png"
cp "$BGO/ten_shards/figures/delta_r50_vs_dose.png"      "$FIGS/dose_sweep_r50_bgo.png"
cp "$CSI/one_shard/fits/figures/scatters_activity.png"  "$FIGS/scatters_profile_csi.png"
cp "$CSI/ten_shards/figures/delta_r50.png"              "$FIGS/ladder_delta_r50_csi.png"
cp "$CSI/ten_shards/figures/delta_r50_vs_dose.png"      "$FIGS/dose_sweep_r50_csi.png"
echo "collected 8 figures into $FIGS/"
