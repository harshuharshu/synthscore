# -*- coding: utf-8 -*-
"""
Structure standardization. RDKit only -- no scorer or numeric dependencies, so
both code/molgen and code/synthesizability can import it without pulling in the
scorer panel.

This module exists because of a defect found in the audit (AUDIT.md, Defect 6.3):
descriptors and scores were computed on whatever structure the SMILES happened to
encode -- arbitrary salt form, arbitrary protonation, arbitrary tautomer. Two
records for the same substance (a free base and its mesylate) therefore received
different descriptors and different scores.

Order follows the RDKit MolVS pipeline:
    Cleanup -> LargestFragmentChooser -> Uncharger -> canonical tautomer
Cleanup itself sanitizes, disconnects metals, normalizes functional groups and
reionizes. Stereochemistry is PRESERVED throughout; the original pipeline
flattened it with isomericSmiles=False while still counting stereocentres
(Defect 6.1).

Tautomer canonicalization is the slow step and is size-guarded: above
`max_atoms_tautomer` heavy atoms it is skipped and the row flagged, rather than
stalling a large run.
"""
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")


class Standardizer:
    """Deterministic structure normalization. Returns (mol, info)."""

    def __init__(self, tautomer=True, max_atoms_tautomer=70):
        self.tautomer = tautomer
        self.max_atoms_tautomer = max_atoms_tautomer
        self._largest = rdMolStandardize.LargestFragmentChooser()
        self._uncharger = rdMolStandardize.Uncharger()
        self._taut = rdMolStandardize.TautomerEnumerator()

    def standardize(self, x):
        info = {"parse_ok": False, "fragment_stripped": False,
                "charge_changed": False, "tautomer_changed": False,
                "tautomer_skipped": False, "error": None}
        mol = x if isinstance(x, Chem.Mol) else Chem.MolFromSmiles(str(x))
        if mol is None:
            info["error"] = "parse"
            return None, info
        info["parse_ok"] = True
        try:
            before = Chem.MolToSmiles(mol)

            mol = rdMolStandardize.Cleanup(mol)

            frag = self._largest.choose(mol)
            if frag is not None:
                if frag.GetNumAtoms() != mol.GetNumAtoms():
                    info["fragment_stripped"] = True
                mol = frag

            pre_charge = Chem.MolToSmiles(mol)
            mol = self._uncharger.uncharge(mol)
            if Chem.MolToSmiles(mol) != pre_charge:
                info["charge_changed"] = True

            if self.tautomer:
                if mol.GetNumHeavyAtoms() <= self.max_atoms_tautomer:
                    pre_taut = Chem.MolToSmiles(mol)
                    mol = self._taut.Canonicalize(mol)
                    if Chem.MolToSmiles(mol) != pre_taut:
                        info["tautomer_changed"] = True
                else:
                    info["tautomer_skipped"] = True

            Chem.SanitizeMol(mol)
            info["changed"] = Chem.MolToSmiles(mol) != before
        except Exception as exc:
            info["error"] = f"{type(exc).__name__}: {str(exc)[:120]}"
            return None, info
        return mol, info

    @staticmethod
    def smiles(mol):
        """Canonical SMILES with stereochemistry retained."""
        return None if mol is None else Chem.MolToSmiles(mol, isomericSmiles=True)

    @staticmethod
    def inchikey(mol):
        try:
            return None if mol is None else Chem.MolToInchiKey(mol)
        except Exception:
            return None
