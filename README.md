# synthscore

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23601853.svg)](https://doi.org/10.5281/zenodo.23601853)

A benchmark for synthesizability prediction, built from reference sets whose makeability is known
independently of any score, together with the code that produced the results in the accompanying paper.

Synthesizability scores gate most generative and virtual-screening workflows, almost always at a
published global cutoff. What those cutoffs discard had not been measured, because measuring it needs
molecules whose route length is known from something other than a score. This repository supplies six
such sets and reports what four published scores and one retrosynthesis planner do with them.

## The headline result

Scoring the six sets with one frozen panel inverts the expected relationship.

| Reference set | n | What is known | SAscore | RAscore | BR-SAScore | RAscore < 0.5 |
|---|---|---|---|---|---|---|
| ZINC purchasable | 3,581 | purchasable, so made | 3.084 | 0.980 | 4.032 | 1.3% |
| Enamine make-on-demand | 1,001 | 2 to 3 steps by construction | 3.505 | 0.940 | 4.508 | 8.0% |
| ChEMBL kinase actives | 3,012 | made, route length unknown | 2.757 | 0.952 | 3.969 | 15.9% |
| **THQ library** | **10,000** | **one step, demonstrated experimentally** | **4.847** | **0.044** | **8.392** | **99.1%** |
| BRICS recombination | 5,000 | no evidence of existence | 3.152 | 0.894 | 4.515 | 16.1% |
| GDB-17 cages | 3,581 | no evidence of existence | 5.730 | 0.096 | 10.000 | 94.8% |

A library demonstrated to be one aza-Diels-Alder step from commercial building blocks is rejected by
RAscore at its published cutoff in 9,907 of 10,000 cases. A fragment-recombination library that has
never been made scores like purchasable compounds. The scores rank molecules by resemblance to
historical catalogue chemistry, not by synthetic accessibility.

Retrosynthesis search reproduces the inversion rather than correcting it. AiZynthFinder, run on 175
molecules sampled evenly from the same sets:

| Reference set | Complete route found | Median steps |
|---|---|---|
| ZINC purchasable | 22/25 (88%) | 1 |
| ChEMBL kinase actives | 14/25 (56%) | 3 |
| Enamine make-on-demand | 11/25 (44%) | 4 |
| BRICS recombination, never made | 10/25 (40%) | 4 |
| THQ made and assayed | 1/25 (4%) | 5 |
| THQ library | 0/25 (0%) | 5 |
| GDB-17 cages | 0/25 (0%) | 5 |

Of 125 molecules known to have been made, 77 were declared unroutable.

## This is an applicability-domain failure, not a badly chosen threshold

Before criticising RAscore's cutoff we tested it on the authors' own published test set of 19,935
AiZynthFinder-labelled ChEMBL compounds. It reaches AUC 0.9530 there, and 0.5 retains 96.47% of
route-solved compounds. The cutoff is sound for the chemistry it was built for. Every tool here rests on
the USPTO reaction corpus, and chemistry absent from that corpus is scored as though it cannot be made.

Note also that RAscore predicts AiZynthFinder's verdicts at AUC 0.909 because it was trained on
AiZynthFinder labels. Their agreement shows a surrogate reproducing its teacher, blind spots included,
rather than either tracking makeability.

## Reference set provenance

| Set | Source |
|---|---|
| ZINC purchasable, GDB-17 cages | SYBA test sets, Vorsilak et al., *J. Cheminform.* 2020 |
| Enamine make-on-demand | sample shipped with synformer |
| THQ library and its 57 assayed members | Papidocha et al., chemRxiv 2026, doi:10.26434/chemrxiv.15010755/v1, CC-BY 4.0 |
| ChEMBL kinase actives | built here from ChEMBL, see `code/synthesizability/build_corpus.py` |
| BRICS recombination | generated here, see `code/molgen/` |

## Layout

```
code/molgen/              de novo generation and medicinal-chemistry filters
code/synthesizability/    scorer panel, standardization, corpus build, route search
data/                     scored reference sets and every result table in the paper
figures/                  the paper's five figures and the graphical abstract
manuscript/               manuscript source
```

## Reproducing

Third-party model weights are not redistributed here. Fetch them from their authors first:

```
python code/synthesizability/fetch_vendor.py
pip install -r requirements.txt
pytest
python code/synthesizability/score_corpus.py
python code/synthesizability/make_figures.py
```

Route search runs separately, in `code/synthesizability/route_search_colab.ipynb`. It builds its own
Python 3.11 environment because AiZynthFinder does not support 3.13.

## Limitations worth reading before reuse

No set here is proven unmakeable. The two sets with no record of synthesis, BRICS recombination and
GDB-17 cages, mean exactly that: no evidence of existence. A true negative set would need compounds that
chemists attempted and failed to make, and no such set is public.

`solved` from a route search is a soft label. AiZynthFinder fails on 44% of Enamine compounds whose
makeability is contractually guaranteed, so its non-solution is not evidence a molecule cannot be made.

Per-target-class operating points are included in `data/calibration_thresholds.csv` and are reported in
the paper as a remedy that does not work. Nine of 36 class thresholds sit at their metric's own scale
limit and admit essentially the entire class. Do not use them as recommended cutoffs.

`code/synthesizability/artifact_analysis.py` and `data/library_structural_scored.csv` implement a
structural signature we tested and rejected; the paper reports why in its discussion. They are kept for
completeness, not as a recommended filter.

## Citation

See `CITATION.cff`. Released under the terms in `LICENSE`. The THQ reference set is reused under
CC-BY 4.0 with attribution to its original authors.
