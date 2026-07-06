"""
crysp_paths.py — the output directory layout, mirroring src/output.jl so the
Python tools write into the same tree the Julia drivers do (dev/PLAN.md →
"Output layout"):

    out/<scenario>/truth/                              scenario tier
    out/<scenario>/<topology>/<ring>/sensitivity/      scenario × ring tier
    out/<scenario>/<topology>/<ring>/<crystal>/        config tier
    out/validation/                                    package cross-checks

Keep this in step with src/output.jl.
"""
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Radiation length X0 (mm) per crystal material — matches src/output.jl.
CRYSTAL_X0_MM = {"BGO": 11.18, "CSI": 18.6}


def crystal_label(material, wall_mm):
    """Crystal directory name folding material and thickness, e.g.
    crystal_label("BGO", 37.0) == "bgo_3X0"."""
    key = material.upper()
    if key not in CRYSTAL_X0_MM:
        raise ValueError(f"no radiation length on file for material {material}")
    return f"{material.lower()}_{round(wall_mm / CRYSTAL_X0_MM[key])}X0"


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
