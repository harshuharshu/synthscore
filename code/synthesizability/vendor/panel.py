# -*- coding: utf-8 -*-
"""
Unified loader for the synthesizability scorer panel.

Every scorer exposes `.score(mol_or_smiles) -> float | None` and carries metadata
(`name`, `direction`, `published_threshold`, `citation`). `direction` is 'lower_easier'
or 'higher_easier', so downstream code never has to remember which way each scale runs.

Scorers are loaded lazily and independently: if one is unavailable the panel degrades
to whatever else is present and records why in `.unavailable`.
"""
import gzip
import os
import pickle
import sys

import numpy as np
from rdkit import Chem
from rdkit.Chem import RDConfig, rdMolDescriptors

_HERE = os.path.dirname(os.path.abspath(__file__))


def _as_mol(x):
    return x if isinstance(x, Chem.Mol) else Chem.MolFromSmiles(str(x))


class _Base:
    name = None
    direction = "lower_easier"
    published_threshold = None
    citation = ""

    def score(self, x):
        raise NotImplementedError


class SAScore(_Base):
    name = "sascore"
    direction = "lower_easier"
    published_threshold = 6.0          # Ertl & Schuffenhauer's own recommendation
    citation = "Ertl & Schuffenhauer, J. Cheminform. 1:8 (2009)"

    def __init__(self):
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer
        self._f = sascorer.calculateScore

    def score(self, x):
        m = _as_mol(x)
        return None if m is None else round(self._f(m), 4)


class NPLikeness(_Base):
    """Not a synthesizability score. Included as the natural-product-likeness axis,
    because the prior art (Ertl 2009) attributes SAscore's class spread to NP character."""
    name = "nplikeness"
    direction = "higher_easier"        # higher = more natural-product-like; no makeability claim
    published_threshold = None
    citation = "Ertl, Roggo & Schuffenhauer, J. Chem. Inf. Model. 48:68 (2008)"

    def __init__(self):
        sys.path.append(os.path.join(RDConfig.RDContribDir, "NP_Score"))
        import npscorer
        self._model = npscorer.readNPModel()
        self._f = npscorer.scoreMol

    def score(self, x):
        m = _as_mol(x)
        return None if m is None else round(self._f(m, self._model), 4)


class SCScore(_Base):
    name = "scscore"
    direction = "lower_easier"
    published_threshold = None         # the paper defines no cutoff
    citation = "Coley et al., J. Chem. Inf. Model. 58:252 (2018)"

    def __init__(self):
        sys.path.insert(0, os.path.join(_HERE, "scscore"))
        import scscore_adapter
        self._m = scscore_adapter.SCScorer()

    def score(self, x):
        m = _as_mol(x)
        return None if m is None else self._m.score_mol(m)


class RAScore(_Base):
    name = "rascore"
    direction = "higher_easier"        # probability that a route exists
    published_threshold = 0.5
    citation = "Thakkar et al., Chem. Sci. 12:3339 (2021)"

    # NOTE: the published model.pkl was serialized by XGBoost 1.x. It loads only under
    # xgboost 1.7 (environment `rascore`), and even there the sklearn wrapper is missing
    # attributes removed after 1.5, so predict_proba() raises. We therefore call the raw
    # Booster. The model's objective is binary:logistic with base_score 0.5 and 2048
    # features, so Booster.predict() returns P(class=1) directly -- the same quantity the
    # authors' predict_proba(...)[0][1] returns.
    def __init__(self):
        import xgboost
        self._xgb = xgboost
        path = os.path.join(_HERE, "rascore", "XGB_chembl_ecfp_counts_model.pkl")
        with open(path, "rb") as fh:
            clf = pickle.load(fh)
        self._bst = clf.get_booster()
        cfg = __import__("json").loads(self._bst.save_config())
        obj = cfg["learner"]["objective"]["name"]
        if obj != "binary:logistic":
            raise RuntimeError(f"unexpected RAscore objective {obj!r}; "
                               "Booster.predict() may not return P(route exists)")

    @staticmethod
    def _ecfp(mol):
        from rdkit.Chem import AllChem
        fp = AllChem.GetMorganFingerprint(mol, 3, useCounts=True, useFeatures=False)
        arr = np.zeros((2048,), np.int32)
        for idx, v in fp.GetNonzeroElements().items():
            arr[idx % 2048] += int(v)
        return arr

    def score(self, x):
        m = _as_mol(x)
        if m is None:
            return None
        d = self._xgb.DMatrix(self._ecfp(m).reshape(1, -1).astype(np.float32))
        return round(float(self._bst.predict(d)[0]), 4)


class BRSAScore(_Base):
    name = "brsascore"
    direction = "lower_easier"
    published_threshold = None
    citation = "Chen et al., J. Cheminform. 16:83 (2024)"

    def __init__(self, frag_penalty=-6.0, complexity_buffer=1.0):
        import math
        self._math = math
        path = os.path.join(_HERE, "brsascore", "pickle",
                            "BRScores_uspto_emolecules.pkl.gz")
        with gzip.open(path, "rb") as fh:
            self._fscores = pickle.load(fh)
        self.frag_penalty = frag_penalty
        self.max_score = 0
        self.min_score = frag_penalty - complexity_buffer

    def score(self, x):
        m = _as_mol(x)
        if m is None:
            return None
        math = self._math
        bi = {}
        fp = rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)
        fps = fp.GetNonzeroElements()

        score1, nf = 0.0, 0
        for bitId, vs in bi.items():
            if vs[0][1] != 2:
                continue
            fscore = self._fscores.get(bitId, self.frag_penalty)
            if fscore < 0:
                nf += 1
                score1 += fscore
        if nf:
            score1 /= nf

        nAtoms = m.GetNumAtoms()
        nChiral = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
        nBridge = rdMolDescriptors.CalcNumBridgeheadAtoms(m)
        nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(m)
        ri = m.GetRingInfo()
        nMacro = sum(1 for r in ri.AtomRings() if len(r) > 6)
        seen = {}
        for r in ri.AtomRings():
            for at in r:
                seen[at] = seen.get(at, 0) + 1
        nMulti = sum(v - 1 for v in seen.values() if v > 1)

        score2 = -(nAtoms ** 1.005 - nAtoms) \
                 - math.log10(nChiral + 1) - math.log10(nSpiro + 1) \
                 - math.log10(nBridge + 1) \
                 - (math.log10(2) if nMacro > 0 else 0) \
                 - math.log10(nMulti + 1)

        score3 = math.log(float(nAtoms) / len(fps)) * 0.5 if nAtoms > len(fps) else 0.0

        raw = score1 + score2 + score3
        s = (raw - self.min_score) / (self.max_score - self.min_score)
        s = min(max(s, 0.0), 1.0)
        return round(10 - s * 9, 4)


_ALL = [SAScore, NPLikeness, SCScore, RAScore, BRSAScore]


class Panel:
    def __init__(self, want=None):
        self.scorers, self.unavailable = {}, {}
        for cls in _ALL:
            if want and cls.name not in want:
                continue
            try:
                self.scorers[cls.name] = cls()
            except Exception as exc:
                self.unavailable[cls.name] = f"{type(exc).__name__}: {str(exc)[:200]}"

    @property
    def names(self):
        return list(self.scorers)

    def score_all(self, x):
        m = _as_mol(x)
        return {n: s.score(m) for n, s in self.scorers.items()}

    def meta(self):
        return {n: dict(direction=s.direction,
                        published_threshold=s.published_threshold,
                        citation=s.citation)
                for n, s in self.scorers.items()}
