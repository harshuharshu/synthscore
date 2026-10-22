# -*- coding: utf-8 -*-
"""
Structure standardization + panel scoring.

This module exists because of a defect found in the audit (AUDIT.md, Defect 6.3):
descriptors and scores in the original pipeline were computed on whatever structure
the SMILES happened to encode -- arbitrary salt form, arbitrary protonation,
arbitrary tautomer. Two records for the same substance (e.g. a free base and its
mesylate) therefore received different scores.

Standardization order follows the RDKit MolVS pipeline:
    Cleanup -> LargestFragmentChooser -> Uncharger -> canonical tautomer
Cleanup itself sanitizes, disconnects metals, normalizes functional groups and
reionizes. Stereochemistry is PRESERVED throughout (the original code flattened it
with isomericSmiles=False while still counting stereocentres -- Defect 6.1).

Tautomer canonicalization is the slow step and is size-guarded: above
`max_atoms_tautomer` heavy atoms it is skipped and the row is flagged, rather than
stalling a large run.
"""
import os
import sys

from rdkit import Chem, RDLogger

from standardize import Standardizer

RDLogger.DisableLog("rdApp.*")

_VENDOR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor")
if _VENDOR not in sys.path:
    sys.path.insert(0, _VENDOR)


class Engine:
    """Standardize, then score with the frozen panel."""

    def __init__(self, panel=None, standardizer=None, want=None):
        import panel as _p
        self.panel = panel or _p.Panel(want=want)
        self.std = standardizer or Standardizer()

    @property
    def metrics(self):
        return self.panel.names

    def score_one(self, smiles, standardize=True):
        if standardize:
            mol, info = self.std.standardize(smiles)
        else:
            mol = Chem.MolFromSmiles(str(smiles))
            info = {"parse_ok": mol is not None, "error": None if mol else "parse"}
        row = {f"std_{k}": v for k, v in info.items()}
        if mol is None:
            row.update({m: None for m in self.metrics})
            row["std_smiles"] = None
            row["inchikey"] = None
            return row
        row["std_smiles"] = Standardizer.smiles(mol)
        row["inchikey"] = Standardizer.inchikey(mol)
        row.update(self.panel.score_all(mol))
        return row

    def score_frame(self, df, smiles_col="smiles", standardize=True, log_every=500):
        import pandas as pd
        out = []
        for i, smi in enumerate(df[smiles_col].tolist(), 1):
            out.append(self.score_one(smi, standardize=standardize))
            if log_every and i % log_every == 0:
                print(f"  scored {i}/{len(df)}", flush=True)
        return pd.concat([df.reset_index(drop=True),
                          pd.DataFrame(out)], axis=1)
