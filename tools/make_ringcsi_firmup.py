#!/usr/bin/env python3
"""
make_ringcsi_firmup.py — persist the N=200 firm-up of the CsI-ring (TBP)
full-window washout σ_R for del120 and del180, which supersede the N=100 values
in sigma_r_washout_v2.toml (del120 washed 0.232→0.226, del180 washed
0.199→0.221; the N=100 del180 was a low fluctuation, del120 held).

The driver uses fixed seeds, so these are deterministic and reproducible with
  julia -t auto --project=. drivers/sigma_r_v2.jl \\
      --realizations 200 --leaves del<t>s_ac300s_1Gy --isotopes none
(config config/run_parameters_csi_v2.toml active). This script writes the
firm-up result file from the retained run products:

- del180: σ_R recomputed EXACTLY from its retained per-realization dump
  (firmup_N200_del180_dump.toml, the sigma_r_grogg_v2 output of that run).
- del120: taken from the recorded N=200 run summary (nominal 0.118, washed
  0.226); its per-realization dump was overwritten, but the value is
  reproducible with the command above (fixed seeds).

del300 is unchanged (N=100, in the main file). Writes
sigma_r_washout_v2_firmup_N200.toml next to the main file. Run:
  python3 tools/make_ringcsi_firmup.py
"""
import os
import tomllib

import numpy as np

BASE = ("out/uniform_headep_sobp_1e8/closed/crysp_ring_1m_csi_2x0/"
        "csi_2X0/washout_v2")

# del120: recorded N=200 run summary (dump not retained; reproducible by command)
DEL120 = {"t_del_s": 120.0, "realizations": 200,
          "nominal_sigma_R_mm": 0.118, "washed_sigma_R_mm": 0.226,
          "source": "N=200 run summary (per-realization dump not retained; "
                    "reproduce with drivers/sigma_r_v2.jl --realizations 200 "
                    "--leaves del120s_ac300s_1Gy --isotopes none)"}


def del180_from_dump():
    with open(os.path.join(BASE, "firmup_N200_del180_dump.toml"), "rb") as f:
        p = tomllib.load(f)["point"][0]
    nom = np.array(p["realizations_erfc_nominal_mm"], float)
    wsh = np.array(p["realizations_erfc_washed_mm"], float)
    return {"t_del_s": 180.0, "realizations": len(nom),
            "nominal_sigma_R_mm": round(float(nom.std(ddof=1)), 4),
            "washed_sigma_R_mm": round(float(wsh.std(ddof=1)), 4),
            "source": "recomputed exactly from firmup_N200_del180_dump.toml "
                      "(retained per-realization dump of the N=200 run)"}


def esc(s):
    return s.replace("\\", "\\\\").replace('"', '\\"')


def main():
    pts = [DEL120, del180_from_dump()]
    lines = [
        'generation = "v2"',
        'scanner = "crysp_ring_1m_csi_2x0"',
        'crystal = "csi"',
        "dose_Gy = 1.0",
        'note = "N=200 firm-up of the CsI-ring full-window washout for del120 '
        "and del180; supersedes those two leaves in sigma_r_washout_v2.toml "
        '(N=100). del300 unchanged (N=100)."',
    ]
    for p in pts:
        lines += [
            "",
            "[[point]]",
            f"t_del_s = {p['t_del_s']}",
            f"realizations = {p['realizations']}",
            f"nominal_sigma_R_mm = {p['nominal_sigma_R_mm']}",
            f"washed_sigma_R_mm = {p['washed_sigma_R_mm']}",
            f'source = "{esc(p["source"])}"',
        ]
    out = os.path.join(BASE, "sigma_r_washout_v2_firmup_N200.toml")
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"wrote {out}")
    for p in pts:
        print(f"  del{int(p['t_del_s'])}: N={p['realizations']} "
              f"nominal {p['nominal_sigma_R_mm']} washed {p['washed_sigma_R_mm']}")


if __name__ == "__main__":
    main()
