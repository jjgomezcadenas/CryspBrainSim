"""
crysp_paths.py — the output directory layout, mirroring src/output.jl so the
Python tools write into the same tree the Julia drivers do (dev/PLAN.md →
"Output layout"):

    out/<scenario>/truth/                              scenario tier
    out/<scenario>/<topology>/<ring>/sensitivity/      scenario × ring tier
    out/<scenario>/<topology>/<ring>/<crystal>/        config tier
    out/validation/                                    package cross-checks

The active scanner arm comes from config/run_parameters.toml [configuration]
(the same block the Julia drivers read): `active_config()` resolves it into
the products paths and the output config directory, so switching arms is an
edit of that one block.

Keep this in step with src/output.jl.
"""
import json
import os
import tomllib
from types import SimpleNamespace

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODUCTS = os.path.join(os.path.dirname(REPO), "PtCryspProds")

# Radiation length X0 (mm) per crystal material — matches src/output.jl.
CRYSTAL_X0_MM = {"BGO": 11.18, "CSI": 18.6}


def crystal_label(crystal, wall_mm):
    """Crystal directory name folding the crystal arm and thickness, e.g.
    crystal_label("BGO", 37.0) == "bgo_3X0", crystal_label("bgo_195k", 22.36)
    == "bgo_195k_2X0". The X0 lookup uses the arm name's first token, so
    operating-point variants keep distinct labels."""
    key = crystal.upper().split("_")[0]
    if key not in CRYSTAL_X0_MM:
        raise ValueError(f"no radiation length on file for material {key} ({crystal})")
    return f"{crystal.lower()}_{round(wall_mm / CRYSTAL_X0_MM[key])}X0"


def out_root():
    return os.path.join(REPO, "out")


def scenario_out(scenario):
    return os.path.join(out_root(), scenario)


def truth_out(scenario):
    return os.path.join(scenario_out(scenario), "truth")


def sensitivity_out(scenario, topology, ring):
    return os.path.join(scenario_out(scenario), topology, ring, "sensitivity")


def config_out(scenario, topology, ring, crystal):
    return os.path.join(scenario_out(scenario), topology, ring, crystal)


def validation_out():
    return os.path.join(out_root(), "validation")


def active_config():
    """The active scanner arm (run_parameters.toml [configuration]) resolved
    into paths: scenario/topology/scanner/crystal/leaf as named upstream, the
    out-tree crystal label (material + thickness from the arm's
    scanner_geometry.json), the output config dir, and the products leaf
    holding the shard files."""
    with open(os.path.join(REPO, "config", "run_parameters.toml"), "rb") as f:
        c = tomllib.load(f)["configuration"]
    scanner_dir = os.path.join(PRODUCTS, c["scenario"], c["scanner"])
    with open(os.path.join(scanner_dir, "scanner_geometry.json")) as f:
        s = json.load(f)["scanner"]
    label = crystal_label(c["crystal"], 10.0 * s["wall_thickness_cm"])
    return SimpleNamespace(
        scenario=c["scenario"], topology=c["topology"], scanner=c["scanner"],
        crystal=c["crystal"], leaf=c["leaf"], label=label,
        cfg_dir=config_out(c["scenario"], c["topology"], c["scanner"], label),
        products_leaf=os.path.join(scanner_dir, c["crystal"], c["leaf"]),
        truth_dir=os.path.join(PRODUCTS, c["scenario"], "truth"),
        scenario_dir=os.path.join(PRODUCTS, c["scenario"]))
