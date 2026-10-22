# Scorer panel (frozen)

Five metrics, all obtained and smoke-tested in this environment. Nothing is cited here that
we could not actually run. Runtime environment: `rascore` (python 3.10, rdkit, xgboost 1.7) --
all five load there, so scoring is single-pass.


## In the panel

| Metric | Direction | Published global threshold | Provenance |
|---|---|---|---|
| SAscore | lower = easier | 6.0 (authors' recommendation) | RDKit contrib, bundled |
| NP-likeness | higher = more NP-like | none; not a makeability score | RDKit contrib, bundled |
| SCScore | lower = easier | none defined by the paper | vendored from connorcoley/scscore |
| RAscore (XGB) | higher = easier, P(route exists) | 0.5 | vendored from reymond-group/RAscore |
| BR-SAScore | lower = easier | none defined by the paper | vendored from awadell1/BR-SAScore |

Citations: SAscore - Ertl & Schuffenhauer, J. Cheminform. 1:8 (2009). NP-likeness - Ertl, Roggo &
Schuffenhauer, JCIM 48:68 (2008). SCScore - Coley et al., JCIM 58:252 (2018). RAscore - Thakkar
et al., Chem. Sci. 12:3339 (2021). BR-SAScore - Chen et al., J. Cheminform. 16:83 (2024).

NP-likeness is not a synthesizability predictor and is never used as one. It is carried because the
prior art (Ertl 2009) attributes SAscore's class-to-class spread to natural-product character, so
NP-likeness is the covariate that explanation must be tested against.


## Excluded, with the reason

| Metric | Why |
|---|---|
| SYBA | Fragment score table ships via Git LFS; the blobs API returns a 133-byte pointer, not the table. Retraining would yield 'our reimplementation of SYBA', not SYBA, so it is excluded rather than approximated. |
| GASA | Requires DGL. Not installable (PyPI unreachable). |
| DFRscore | Checkpoint is obtainable (1.9 MB) but needs PyTorch, and its retro-analysis data is a 39 MB tarball. Deferred as optional, not ruled out. |
| FSscore | Not obtainable without PyPI. |
| AiZynthFinder (route-search ground truth) | PyPI-only package. Zenodo (weights) is reachable now; only the package blocks it. |

### Network and install constraints, as measured

- PyPI unreachable: `502 Bad Gateway` on the CONNECT tunnel then read timeouts, 5 retries, confirmed twice ~25 min apart.
- `git clone` and `raw.githubusercontent.com` blocked; `api.github.com` reachable and token-authenticated, so vendored files came through the git-blobs endpoint.
- conda-forge reachable (all environments and rdkit/xgboost/sklearn came from there).
- zenodo.org granted this session.


## Compatibility fixes to vendored code

Published weights are used unmodified. All three fixes are mechanical.

- **SCScore**: the authors' `standalone_model_numpy.py` uses `np.bool` (removed in NumPy 1.24) and
  `six.iteritems`. `vendor/scscore/scscore_adapter.py` reproduces their forward pass exactly
  (alternating weight/bias pairs, ReLU on all but the last layer, then `1 + (scale-1)*sigmoid`).
  Loaded weight shapes confirm the published architecture: 1024 -> 300 -> 300 -> 300 -> 300 -> 1.
- **BR-SAScore**: the authors' class resolves its pickle via `pkg_resources.resource_filename`,
  which fails outside an installed package. The adapter loads the same
  `BRScores_uspto_emolecules.pkl.gz` from the vendored path; scoring arithmetic is verbatim.
- **RAscore**: `model.pkl` was serialized by XGBoost 1.x. Modern XGBoost refuses it outright
  (`Invalid serialization file`); under xgboost 1.7 it unpickles but the sklearn wrapper raises
  `AttributeError` on attributes removed after 1.5 (`use_label_encoder`, then `callbacks`).
  The adapter therefore calls the raw `Booster`. The model's saved config reports
  `binary:logistic`, `base_score` 0.5, `num_feature` 2048, so `Booster.predict()` returns
  P(class=1) -- the same quantity the authors' `predict_proba(...)[0][1]` returns. Verified
  identical to the authors' code path on a test molecule.


## Smoke test

18 approved drugs sampled at even intervals across the SAscore range, plus the generated molecule
`cand_004995`. Structures read from `data/approved_drugs_scored.csv` and `data/jak_library.csv`;
no SMILES hand-entered. All five metrics returned a value for all 19 molecules (zero nulls).
Full output: `data/panel_smoketest.csv`.

Three observations that shape the rest of the study:

1. **The metrics are not redundant.** At nearly identical SAscore, ketoconazole (3.437) and
   `cand_004995` (3.447) separate on the others: SCScore 4.329 vs 4.824, BR-SAScore 3.958 vs 4.975,
   RAscore 0.995 vs 0.352. A panel sees a difference that SAscore alone does not.
2. **RAscore rejects marketed drugs.** At its published 0.5 cutoff it calls ceftriaxone (0.011),
   pergolide (0.041), trilostane (0.051), oritavancin (0.052) and oxymetholone (0.103)
   route-inaccessible. Ceftriaxone is a high-volume generic antibiotic. These are false negatives
   on chemistry that demonstrably exists, and they cluster in beta-lactam and steroid scaffolds --
   which is the target-class calibration question this study is built to measure.
3. **BR-SAScore spreads wider than SAscore and saturates.** It hits its ceiling of 10 for both
   oritavancin and trilostane where SAscore gives 8.98 and 5.32. Whether the extra spread is signal
   or scale artefact is a calibration question, not a claim to make from 19 molecules.
