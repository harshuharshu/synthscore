# synthscore — target-class calibrated synthesizability scoring

Synthesizability scores ship with a single global cutoff. That cutoff is wrong for
most target classes, and measurably so. Scored against 26,947 ChEMBL compounds that
were demonstrably synthesized (each carries a measured activity), SAscore at its
published cutoff of 6.0 discards **9.5% of known antibacterial chemistry but only
0.29% of known nuclear-receptor chemistry** — a 33-fold difference driven by the
target class, not by the molecule. RAscore at its published 0.5 discards **41.6%
of antibacterial chemistry** and 12.8% of oxidoreductase chemistry.

This package scores molecules with five published metrics and judges each against
a cutoff calibrated for the target class you name.

## What "calibrated" means here

The cutoff that discards a declared fraction (default 5%) of compounds in that
target class which were actually made. The positive label derives from
*experimental record*, never from a synthesizability score, so the calibration is
not circular.

There is no negative class. Nothing available is *proven* unmakeable, so the tool
reports a false-reject rate on real chemistry and does not claim a false-accept
rate. See `LIMITATIONS` below.

## Install

```
conda create -n synthscore -c conda-forge python=3.10 rdkit pandas numpy scikit-learn xgboost=1.7 matplotlib pytest
conda activate synthscore
python code/synthesizability/fetch_vendor.py
```

`fetch_vendor.py` downloads three model-weight files (SCScore, RAscore,
BR-SAScore) from their authors' own repositories. They are not redistributed
here and remain under their authors' licences. SAscore and NP-likeness ship with
RDKit, so the tool runs — with a warning and a reduced panel — even if the fetch
fails.

`xgboost=1.7` is load-bearing: RAscore's published model was serialized by
XGBoost 1.x and modern XGBoost refuses to unpickle it. Everything else runs on
any recent RDKit. Metrics that fail to load are skipped with a warning rather
than aborting the run.

## Use

```
python code/synthesizability/cli.py --list-classes

python code/synthesizability/cli.py \
    --in molecules.csv --smiles-col smiles \
    --target-class kinase --out scored.csv
```

Options: `--panel sascore,rascore` to restrict the metric set,
`--false-reject 0.10` for a looser operating point, `--no-standardize` to
reproduce pre-standardization numbers, `--id-col` to carry an identifier through.

Output carries, per molecule: the standardized SMILES and InChIKey, every metric,
the class cutoff used, a per-metric reject flag, a reject flag against the
published global cutoff where one exists, and `reject_any`.

## The metric panel

| Metric | Direction | Published cutoff | Provenance |
|---|---|---|---|
| SAscore | lower = easier | 6.0 | RDKit contrib (bundled) |
| NP-likeness | higher = more NP-like | none; **never used as a gate** | RDKit contrib (bundled) |
| SCScore | lower = easier | none | vendored, connorcoley/scscore |
| RAscore (XGB) | higher = easier | 0.5 | vendored, reymond-group/RAscore |
| BR-SAScore | lower = easier | none | vendored, awadell1/BR-SAScore |

Full provenance, the compatibility fixes applied to vendored code, and the
metrics that could not be obtained (SYBA, GASA, FSscore, DFRscore,
AiZynthFinder) are in `code/synthesizability/scorer_panel.md`.

## Standardization

Every molecule passes through `standardize.py` before scoring: RDKit `Cleanup`,
largest-fragment selection, charge neutralization, canonical tautomer, with
stereochemistry retained. This is not cosmetic. Across 730 approved-drug
salt/parent pairs, **0 of 730 scored identically before standardization** (max
difference 2.59 SAscore units, 5.35 BR-SAScore units); 724 of 730 collapse to one
structure and one score after it. The 19 covalent ester prodrugs in the same
comparison correctly stay distinct. See `data/standardization_effect.csv` and
`code/synthesizability/tests/`.

## Reproducing

```
python code/synthesizability/build_corpus.py     # ~40 min, ChEMBL API
python code/synthesizability/score_corpus.py     # ~37 min, 5 metrics x 27k compounds
python code/synthesizability/make_figures.py     # tables + 4 figures, ~40 s
python -m pytest                                  # 5 regression tests
```

`build_corpus.py` caches `data/target_class_map.csv`, so a re-run skips the
12,912-component classification scan.

`make_figures.py` regenerates `class_summary.csv`, `calibration_thresholds.csv`,
`global_vs_calibrated.csv`, `discrimination_auc.csv`, `jak_recalibrated.csv` and
the figures `class_distributions.png`, `roc_panel.png`,
`cost_of_global_threshold.png`, `jak_recalibrated.png`. It reads only the scored
CSVs, so it does **not** need the scorer panel and will run in any environment
with pandas, matplotlib, scikit-learn and RDKit — useful because the
`xgboost=1.7` pin is only required for scoring.

## LIMITATIONS

- **No proven-unmakeable negatives.** Calibration is one-class. The tool bounds
  how much real chemistry a cutoff discards; it cannot tell you how much
  unmakeable chemistry it admits.
- **No route-search ground truth.** AiZynthFinder is PyPI-only and PyPI was
  unreachable in this environment, so no metric here is validated against an
  actual retrosynthesis search.
- **Corpus membership is evidence of synthesis, not of ease.** A compound in
  ChEMBL was made by someone, possibly over 14 steps.
- **Class assignment is ChEMBL's.** Compounds active against several classes
  appear in each (513 of 26,947).
- **The panel does not reproduce expert rejection.** On the one molecule a
  synthetic chemist rejected by hand (`cand_004995`), every class-calibrated
  cutoff passes it. A cutoff strict enough to catch it costs 7.4% (RAscore) to
  29.5% (BR-SAScore) of known kinase chemistry. Reported in full rather than
  tuned away.
