#!/usr/bin/env python3
"""
shard_summary.py — attrs + truth breakdown of stored LOR shards. Run this
first on every new shard (validation ladder rung 2, the text half; the 3×3
figure is plot_shard.py, the assertable struct src/qa.jl).

Run:  python3 tools/shard_summary.py <shard.h5> [<shard.h5> ...]
      python3 tools/shard_summary.py <leaf_dir>          # all shards in a leaf
"""
import glob
import os
import sys

import h5py
import numpy as np


def dec(v):
    return v.decode() if isinstance(v, bytes) else v


def summarize(path):
    with h5py.File(path, "r") as f:
        a = {k: dec(v) for k, v in f.attrs.items()}
        truth = f["truth"][:]
        nscat = f["nscat1"][:].astype(np.int64) + f["nscat2"][:].astype(np.int64)

    n = truth.size
    n_t = int((truth == 0).sum())
    n_s = int((truth == 1).sum())
    n_r = int((truth == 2).sum())
    n_single = int(((truth == 1) & (nscat == 1)).sum())
    n_multi = int(((truth == 1) & (nscat >= 2)).sum())
    nev = int(a.get("nevents", 0))

    print(f"=== {os.path.basename(path)}")
    print(f"  {a.get('scenario', '?')} / {a.get('crystal', '?')} / "
          f"{a.get('budget', '?')} — shard {a.get('realization', '?')} "
          f"(master_seed {a.get('master_seed', '?')})")
    print(f"  nrows {n:,} / nevents {nev:,} — acceptance {n / max(nev, 1):.2%}")
    print(f"  true {n_t / n:.2%} | scatter {n_s / n:.2%} "
          f"(single {n_single / n:.2%}, multiple {n_multi / n:.2%}) | "
          f"random {n_r / n:.3%}")
    print(f"  detector: eres {a.get('eres', '?')}, Emin {a.get('emin_keV', '?')} keV, "
          f"σ_xyz {a.get('sigma_xyz_mm', '?')} mm, τ {a.get('tau_ns', '?')} ns, "
          f"blocks {a.get('n_phi', '?')}×{a.get('n_z', '?')}")
    print(f"  dose {a.get('dose_Gy', '?')} Gy, t_meas {a.get('t_meas_s', '?')} s, "
          f"escaped dropped {a.get('n_escaped_dropped', '?')}")
    return n, n_t, n_s, n_r, nev


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    paths = []
    for arg in sys.argv[1:]:
        if os.path.isdir(arg):
            paths += sorted(glob.glob(os.path.join(arg, "lors_shard*.h5")))
        else:
            paths.append(arg)
    if not paths:
        sys.exit("no shards found")

    tot = np.zeros(5, dtype=np.int64)
    for p in paths:
        tot += np.array(summarize(p), dtype=np.int64)
    if len(paths) > 1:
        n, n_t, n_s, n_r, nev = tot
        print(f"=== pooled: {len(paths)} shards")
        print(f"  nrows {n:,} / nevents {nev:,} — acceptance {n / max(nev, 1):.2%}")
        print(f"  true {n_t / n:.2%} | scatter {n_s / n:.2%} | random {n_r / n:.3%}")


if __name__ == "__main__":
    main()
