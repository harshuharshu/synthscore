# -*- coding: utf-8 -*-
"""
medchem.py - property filters and desirability scoring.

Filters run cheap->expensive. Objectives are combined with a Derringer
desirability function (geometric mean of per-property desirabilities) rather
than a weighted sum, so a single dealbreaker property cannot be averaged away.

Dependencies: RDKit only (including its bundled Ertl SA scorer and
FilterCatalogs), plus code/synthesizability/standardize.py.

Note on thresholds: the property windows in PROFILES are convention-derived
defaults, not calibrated values. In particular `sa_max` is an arbitrary cutoff;
using it as a hard gate makes any downstream statement about SAscore's behaviour
on the resulting library circular (AUDIT.md, Defect 1). Calibrated, target-class
specific thresholds live in code/synthesizability/.
"""

import os
import sys
import math

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import FilterCatalog
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.Chem import AllChem

# RDKit ships the Ertl synthetic-accessibility scorer as a contrib module.
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402


# ---------------------------------------------------------------------------
# Structural alert catalogs (PAINS + Brenk + NIH). Built once, reused.
# ---------------------------------------------------------------------------
def _build_alert_catalog():
    params = FilterCatalog.FilterCatalogParams()
    for cat in (
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
        FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH,
    ):
        params.AddCatalog(cat)
    return FilterCatalog.FilterCatalog(params)


_ALERT_CATALOG = _build_alert_catalog()


def structural_alerts(mol):
    """Return list of structural-alert names the molecule triggers (PAINS/Brenk/NIH)."""
    matches = _ALERT_CATALOG.GetMatches(mol)
    return [m.GetDescription() for m in matches]


# ---------------------------------------------------------------------------
# Property panel - the descriptors a med chemist actually looks at.
# ---------------------------------------------------------------------------
_STD_PATH = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "..", "synthesizability"))
if _STD_PATH not in sys.path:
    sys.path.insert(0, _STD_PATH)
from standardize import Standardizer  # noqa: E402

_STANDARDIZER = Standardizer()


def compute_properties(mol, standardize=True):
    """
    Descriptor panel. By default the structure is standardized first (salt
    stripped, charges neutralized, canonical tautomer, stereochemistry retained),
    so that a salt and its parent give identical descriptors. Pass
    standardize=False only to reproduce pre-audit numbers.
    """
    if standardize:
        std, _info = _STANDARDIZER.standardize(mol)
        if std is not None:
            mol = std
    return {
        "mw": Descriptors.MolWt(mol),
        "clogp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "n_rings": mol.GetRingInfo().NumRings(),
        # CalcNumAtomStereoCenters counts only ASSIGNED centres, so it reads 0 on a
        # flat structure even where stereochemistry is possible. The unassigned
        # count is reported alongside so the pair is interpretable.
        "n_stereo": rdMolDescriptors.CalcNumAtomStereoCenters(mol),
        "n_stereo_unassigned": rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol),
        "qed": QED.qed(mol),
        "sa": sascorer.calculateScore(mol),
        "n_heavy": mol.GetNumHeavyAtoms(),
    }


# ---------------------------------------------------------------------------
# Hard property windows. "balanced" = generous lead/drug-like windows;
# anything outside is a real liability, not a stylistic preference.
# ---------------------------------------------------------------------------
PROFILES = {
    "strict": dict(mw=(250, 450), clogp=(-0.5, 4.0), tpsa=(40, 110), hbd=(0, 4),
                   hba=(1, 8), rotb=(0, 7), fsp3_min=0.30, arom_max=3,
                   max_alerts=0, sa_max=4.0),
    "balanced": dict(mw=(200, 500), clogp=(-1.0, 5.0), tpsa=(20, 130), hbd=(0, 5),
                     hba=(0, 10), rotb=(0, 10), fsp3_min=0.20, arom_max=4,
                     max_alerts=1, sa_max=5.0),
    "permissive": dict(mw=(150, 650), clogp=(-3.0, 6.5), tpsa=(0, 180), hbd=(0, 7),
                       hba=(0, 12), rotb=(0, 15), fsp3_min=0.0, arom_max=6,
                       max_alerts=99, sa_max=7.0),
}


def passes_windows(props, alerts, profile):
    """Return (bool_pass, list_of_reasons). In permissive mode, reasons are advisory."""
    P = PROFILES[profile]
    reasons = []
    if not (P["mw"][0] <= props["mw"] <= P["mw"][1]):
        reasons.append("MW %.0f outside %s" % (props["mw"], P["mw"]))
    if not (P["clogp"][0] <= props["clogp"] <= P["clogp"][1]):
        reasons.append("cLogP %.1f outside %s" % (props["clogp"], P["clogp"]))
    if not (P["tpsa"][0] <= props["tpsa"] <= P["tpsa"][1]):
        reasons.append("TPSA %.0f outside %s" % (props["tpsa"], P["tpsa"]))
    if props["hbd"] > P["hbd"][1]:
        reasons.append("HBD %s > %s" % (props["hbd"], P["hbd"][1]))
    if props["hba"] > P["hba"][1]:
        reasons.append("HBA %s > %s" % (props["hba"], P["hba"][1]))
    if props["rotb"] > P["rotb"][1]:
        reasons.append("RotB %s > %s" % (props["rotb"], P["rotb"][1]))
    if props["fsp3"] < P["fsp3_min"]:
        reasons.append("Fsp3 %.2f < %s (too flat)" % (props["fsp3"], P["fsp3_min"]))
    if props["n_aromatic_rings"] > P["arom_max"]:
        reasons.append("AromRings %s > %s" % (props["n_aromatic_rings"], P["arom_max"]))
    if props["sa"] > P["sa_max"]:
        reasons.append("SA %.1f > %s (hard to make)" % (props["sa"], P["sa_max"]))
    if len(alerts) > P["max_alerts"]:
        reasons.append("%d structural alerts > %s" % (len(alerts), P["max_alerts"]))
    return (len(reasons) == 0, reasons)


# ---------------------------------------------------------------------------
# Derringer desirability - the right way to fuse objectives.
# Each property maps to a desirability d in [0,1]; overall D is their geometric
# mean, so any single d=0 zeroes the molecule (a true dealbreaker), which a
# weighted sum would wrongly let slide.
# ---------------------------------------------------------------------------
def _d_target(x, lo, best_lo, best_hi, hi):
    """Trapezoidal desirability: 1.0 in [best_lo,best_hi], ramps to 0 at lo/hi."""
    if x <= lo or x >= hi:
        return 0.0
    if best_lo <= x <= best_hi:
        return 1.0
    if x < best_lo:
        return (x - lo) / (best_lo - lo)
    return (hi - x) / (hi - best_hi)


def _d_max(x, lo, hi):
    """Higher is better, saturating at hi."""
    if x <= lo:
        return 0.0
    if x >= hi:
        return 1.0
    return (x - lo) / (hi - lo)


def _d_min(x, lo, hi):
    """Lower is better, saturating at lo."""
    if x <= lo:
        return 1.0
    if x >= hi:
        return 0.0
    return (hi - x) / (hi - lo)


def desirability(props):
    """Composite candidate-quality score in [0,1] (geometric mean of components)."""
    d = {
        "mw":    _d_target(props["mw"], 180, 280, 420, 560),
        "clogp": _d_target(props["clogp"], -1.5, 1.0, 3.5, 5.5),
        "tpsa":  _d_target(props["tpsa"], 15, 50, 100, 150),
        "fsp3":  _d_max(props["fsp3"], 0.10, 0.45),       # reward 3D character
        "qed":   _d_max(props["qed"], 0.30, 0.80),
        "sa":    _d_min(props["sa"], 2.5, 6.0),           # reward synthesizability
        "rotb":  _d_min(props["rotb"], 6, 12),            # penalize floppiness
    }
    vals = [max(v, 1e-6) for v in d.values()]  # avoid log(0); keep dealbreaker signal via tiny value
    logmean = sum(math.log(v) for v in vals) / len(vals)
    D = math.exp(logmean)
    # if any true zero, force D toward zero (geometric-mean dealbreaker behavior)
    if min(d.values()) == 0.0:
        D *= 0.05
    return round(D, 4), {k: round(v, 3) for k, v in d.items()}


# ---------------------------------------------------------------------------
# Novelty utilities
# ---------------------------------------------------------------------------
def murcko_scaffold(mol):
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaf)
    except Exception:
        return ""


def morgan_fp(mol, radius=2, nbits=2048):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def max_tanimoto(mol, ref_fps):
    """Max Tanimoto similarity of mol to a list of reference fingerprints."""
    if not ref_fps:
        return 0.0
    fp = morgan_fp(mol)
    sims = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
    return max(sims) if sims else 0.0
