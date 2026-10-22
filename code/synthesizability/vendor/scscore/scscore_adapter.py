# -*- coding: utf-8 -*-
"""
Thin loader for the published SCScore model (Coley et al., J. Chem. Inf. Model. 2018).

Weights are the authors' own `model.ckpt-10654.as_numpy.json.gz` from
github.com/connorcoley/scscore, vendored unmodified under vendor/scscore/.

The forward pass below is identical to the authors' `standalone_model_numpy.py`
(alternating weight/bias pairs, ReLU on all but the last layer, then
1 + (scale-1)*sigmoid). Two compatibility fixes only, no behavioural change:
  * `np.bool` -> `bool`   (np.bool was removed in NumPy 1.24)
  * `six.iteritems` -> `.items()`
Fingerprint settings are the authors': Morgan radius 2, 1024 bits, useChirality=True.
"""
import gzip
import json
import os

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

SCORE_SCALE = 5.0
FP_LEN = 1024
FP_RAD = 2

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_WEIGHTS = os.path.join(_HERE, "model.ckpt-10654.as_numpy.json.gz")


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


class SCScorer:
    def __init__(self, weight_path=DEFAULT_WEIGHTS, fp_len=FP_LEN, fp_rad=FP_RAD,
                 score_scale=SCORE_SCALE):
        self.fp_len = fp_len
        self.fp_rad = fp_rad
        self.score_scale = score_scale
        with gzip.GzipFile(weight_path, "r") as fin:
            self.vars = [np.array(x) for x in json.loads(fin.read().decode("utf-8"))]

    def mol_to_fp(self, mol):
        if mol is None:
            return np.zeros((self.fp_len,), dtype=np.float32)
        bv = AllChem.GetMorganFingerprintAsBitVect(
            mol, self.fp_rad, nBits=self.fp_len, useChirality=True)
        return np.array(bv, dtype=bool).astype(np.float32)

    def apply(self, x):
        for i in range(0, len(self.vars), 2):
            last = (i == len(self.vars) - 2)
            x = np.matmul(x, self.vars[i]) + self.vars[i + 1]
            if not last:
                x = x * (x > 0)          # ReLU
        return float(1.0 + (self.score_scale - 1.0) * _sigmoid(x).item())

    def score(self, smiles):
        """Return SCScore in [1, 5], or None if the SMILES does not parse."""
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        fp = self.mol_to_fp(mol)
        if fp.sum() == 0:
            return None
        return round(self.apply(fp), 4)

    def score_mol(self, mol):
        if mol is None:
            return None
        fp = self.mol_to_fp(mol)
        if fp.sum() == 0:
            return None
        return round(self.apply(fp), 4)
