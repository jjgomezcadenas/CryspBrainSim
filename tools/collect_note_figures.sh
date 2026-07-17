#!/bin/bash
# collect_note_figures.sh — copy the figures latex/endpoint_precision.tex
# includes from the (git-ignored) output tree into latex/figs/, so the note
# compiles from a fresh clone. Both scanner arms are collected regardless of
# the active [configuration]. Regenerate the sources first when needed (per
# arm, with that arm active):
#   tools/plot_activity_time.py                       (induced-activity vs time)
#   tools/fit_activity_profile.py --no-pulls          (fit figure)
#   tools/plot_recon_projections.py                   (projections)
#   tools/ten_shards.py  /  --dose-sweep              (ladder, dose sweep)
#   tools/recon_scatters.jl + tools/scatter_profile.py (scatter profile)
#   drivers/sigma_r_v2.jl + tools/plot_sigma_r_v2.py  (v2 washout + per-isotope,
#                                                      ring CsI v2 arm active)
set -euo pipefail
cd "$(dirname "$0")/.."

FIGS=latex/figs
TRUTH=out/uniform_headep_sobp_1e8/truth
BGO=out/uniform_headep_sobp_1e8/closed/crysp_ring_1m_bgo_2x0/bgo_195k_2X0
CSI=out/uniform_headep_sobp_1e8/closed/crysp_ring_1m_csi_2x0/csi_2X0

cp "$TRUTH/figures/activity_time.png"                   "$FIGS/activity_time.png"

cp "$BGO/one_shard/fits/figures/recon_all_events_activity.png" "$FIGS/fit_all_events_bgo.png"
cp "$CSI/one_shard/fits/figures/recon_all_events_activity.png" "$FIGS/fit_all_events_csi.png"
cp "$BGO/one_shard/fits/figures/recon_projections.png"  "$FIGS/recon_projections_bgo.png"
cp "$BGO/one_shard/fits/figures/scatters_activity.png"  "$FIGS/scatters_profile_bgo.png"
cp "$BGO/ten_shards/figures/delta_r50.png"              "$FIGS/ladder_delta_r50_bgo.png"
cp "$BGO/ten_shards/figures/delta_r50_vs_dose.png"      "$FIGS/dose_sweep_r50_bgo.png"
cp "$CSI/one_shard/fits/figures/scatters_activity.png"  "$FIGS/scatters_profile_csi.png"
cp "$CSI/ten_shards/figures/delta_r50.png"              "$FIGS/ladder_delta_r50_csi.png"
cp "$CSI/ten_shards/figures/delta_r50_vs_dose.png"      "$FIGS/dose_sweep_r50_csi.png"
cp out/uniform_headep_sobp_1e8/closed/comparison/figures/tstart_r50.png "$FIGS/tstart_r50.png"
CMP=out/uniform_headep_sobp_1e8/closed/comparison/figures
cp "$CSI/washout_v2/figures/sigma_r_v2.png"     "$FIGS/sigma_r_v2_csi.png"
cp "$BGO/washout_v2/figures/sigma_r_v2.png"     "$FIGS/sigma_r_v2_bgo.png"
cp "$CMP/washed_sigma_r_scanners_v2_csi.png"    "$FIGS/washed_sigma_r_scanners_v2_csi.png"
cp "$CMP/washed_sigma_r_scanners_v2_bgo.png"    "$FIGS/washed_sigma_r_scanners_v2_bgo.png"
cp "$CMP/washed_bgo_vs_csi_v2.png"              "$FIGS/washed_bgo_vs_csi_v2.png"
cp "$CMP/grogg_v2.png"                          "$FIGS/grogg_v2.png"
echo "collected 17 figures into $FIGS/"
