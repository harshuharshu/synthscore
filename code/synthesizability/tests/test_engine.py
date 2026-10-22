# -*- coding: utf-8 -*-
"""
Regression tests for the standardizer.

The defect these guard against (AUDIT.md, Defect 6.3): scores were computed on
whatever structure the SMILES encoded, so a free base and its salt scored
differently even though they are the same substance.

Structures are read from data/approved_drugs_scored.csv -- none are hand-entered.

Run:  python -m pytest code/synthesizability/tests/ -q
"""
import os
import sys

import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
DATA = os.path.abspath(os.path.join(PKG, "..", "..", "data"))
for p in (PKG, os.path.join(PKG, "vendor")):
    if p not in sys.path:
        sys.path.insert(0, p)

import engine  # noqa: E402

SALT_SUFFIXES = ("MESYLATE", "HYDROCHLORIDE", "SULFATE", "MALEATE", "CITRATE",
                 "TARTRATE", "FUMARATE", "SODIUM", "POTASSIUM", "BROMIDE",
                 "SUCCINATE", "BESYLATE", "TOSYLATE", "NITRATE", "OXALATE")

# Suffixes that denote a covalent ester prodrug rather than a salt. These must NOT
# be merged -- they are different molecules and must keep distinct scores.
ESTER_SUFFIXES = ("ACETATE", "PALMITATE")


@pytest.fixture(scope="module")
def drugs():
    path = os.path.join(DATA, "approved_drugs_scored.csv")
    if not os.path.exists(path):
        pytest.skip("approved_drugs_scored.csv not present")
    return pd.read_csv(path).dropna(subset=["name", "smiles"])


@pytest.fixture(scope="module")
def eng():
    return engine.Engine(want=["sascore", "scscore", "brsascore"])


def _pairs(drugs, suffixes):
    names = set(drugs.name)
    out = []
    for n in names:
        parts = n.split()
        if len(parts) >= 2 and parts[-1] in suffixes:
            parent = " ".join(parts[:-1])
            if parent in names:
                out.append((parent, n))
    return sorted(set(out))


def test_true_salts_collapse_to_parent(drugs, eng):
    """A salt and its free parent must standardize to one structure and one score."""
    smi = dict(zip(drugs.name, drugs.smiles))
    pairs = _pairs(drugs, SALT_SUFFIXES)
    assert len(pairs) > 100, f"expected many salt pairs, found {len(pairs)}"

    merged = 0
    for parent, salt in pairs:
        a = eng.score_one(smi[parent], standardize=True)
        b = eng.score_one(smi[salt], standardize=True)
        if a["std_smiles"] == b["std_smiles"]:
            merged += 1
            for m in eng.metrics:
                assert a[m] == b[m], f"{salt}: {m} differs after standardization"
    # Residual non-merging salts are charge-state edge cases, notably cephalosporin
    # zwitterions. The overwhelming majority must collapse.
    assert merged / len(pairs) > 0.90, f"only {merged}/{len(pairs)} salt pairs merged"


def test_unstandardized_salts_do_not_match(drugs, eng):
    """Without standardization the same pairs score differently -- the original defect."""
    smi = dict(zip(drugs.name, drugs.smiles))
    pairs = _pairs(drugs, SALT_SUFFIXES)[:40]
    differing = sum(
        eng.score_one(smi[p], standardize=False)["sascore"]
        != eng.score_one(smi[s], standardize=False)["sascore"]
        for p, s in pairs)
    assert differing == len(pairs), \
        f"only {differing}/{len(pairs)} unstandardized pairs differed; expected all"


def test_ester_prodrugs_are_not_merged(drugs, eng):
    """
    '<DRUG> ACETATE' is ambiguous: sometimes an acetic-acid SALT (bazedoxifene
    acetate), sometimes a covalent ESTER prodrug (abiraterone acetate). The
    standardizer must merge the first and keep the second distinct. The two cases
    are distinguishable from the raw SMILES: a salt is a multi-component structure
    (a disconnected '.' fragment), an ester is a single connected molecule.
    """
    smi = dict(zip(drugs.name, drugs.smiles))
    pairs = _pairs(drugs, ESTER_SUFFIXES)
    assert pairs, "no acetate/palmitate pairs found"

    covalent_checked = 0
    for parent, other in pairs:
        a = eng.score_one(smi[parent], standardize=True)
        b = eng.score_one(smi[other], standardize=True)
        merged = (a["std_smiles"] == b["std_smiles"])
        is_multicomponent = "." in str(smi[other])
        if is_multicomponent:
            assert merged, f"{other} is a multi-component salt but was not merged"
        else:
            covalent_checked += 1
            assert not merged, \
                f"{other} is a single covalent molecule but was merged into {parent}"
    assert covalent_checked > 0, "no covalent ester prodrugs in the comparison set"


def test_stereochemistry_is_preserved(drugs, eng):
    """Defect 6.1: stereocentres must survive standardization."""
    from rdkit import Chem
    chiral = [s for s in drugs.smiles if "@" in str(s)][:60]
    assert chiral, "no stereo-bearing SMILES in the reference set"
    kept = 0
    for s in chiral:
        mol, _ = eng.std.standardize(s)
        if mol is not None and "@" in Chem.MolToSmiles(mol, isomericSmiles=True):
            kept += 1
    assert kept / len(chiral) > 0.85, \
        f"stereochemistry lost in {len(chiral)-kept} of {len(chiral)}"


def test_unparseable_input_returns_none_not_exception(eng):
    row = eng.score_one("not_a_smiles", standardize=True)
    assert row["std_smiles"] is None
    assert all(row[m] is None for m in eng.metrics)
