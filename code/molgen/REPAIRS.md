# molgen correctness repairs

Four defects from AUDIT.md section 6, fixed in `code/molgen/` and `code/synthesizability/`.


## What was changed

| Defect | Fix | File |
|---|---|---|
| 6.1 Stereochemistry flattened while `n_stereo` was counted | Both emit paths now use `isomericSmiles=True`. `compute_properties` additionally reports `n_stereo_unassigned`. | `generators.py`, `medchem.py` |
| 6.2 RNN not reproducible | `CharRNNGenerator` takes `seed` (default 42), seeds `torch`, CUDA and its own `random.Random`; training shuffles with that instance instead of the global RNG; `fit()` reseeds before weight init; optional `deterministic=True` constrains cuDNN. | `generators.py` |
| 6.3 No structure standardization | `compute_properties(mol, standardize=True)` routes through `standardize.Standardizer` (Cleanup, largest fragment, uncharge, canonical tautomer, stereochemistry retained). | `medchem.py`, `standardize.py` |
| 6.4 Overselling docstrings | 'This is where candidate quality is actually decided', 'IP-clean', 'the whole pipeline is reproducible and yours' removed. Replaced with the measured BRICS polycyclic-enrichment figure and a note that `sa_max` is an arbitrary cutoff. | `generators.py`, `medchem.py` |

`Standardizer` was moved out of `engine.py` into `standardize.py` so `molgen` can import it
without pulling in the scorer panel (numpy, xgboost, vendored models).


## Effect on the published JAK numbers: smaller than expected

All 5,000 molecules were re-scored with `standardize=False` (pre-audit behaviour) and
`standardize=True`. Thirteen of the fourteen property columns are **bit-identical**:

| Property | Molecules changed | Mean delta | Max abs delta |
|---|---|---|---|
| mw, clogp, tpsa, hbd, hba, rotb, fsp3, n_aromatic_rings, n_rings, n_stereo, qed, n_heavy | 0 of 5,000 | 0 | 0 |
| sa (SAscore) | 1,112 of 5,000 (22.2%) | -0.0005 | 0.5299 |

The reason is structural: BRICS recombination emits neutral, single-component molecules, so
there is no salt to strip and no charge to neutralize. The only standardization step that bites
is canonical tautomerization, and it moves SAscore for about a fifth of the library with a mean
shift of essentially zero.

**This is a negative result and it should be reported as one.** The standardization defect is
real and it demonstrably matters -- on the approved-drug reference set it moved SAscore by up to
2.59 and BR-SAScore by up to 5.35 between a salt and its own parent (`standardization_effect.csv`).
It did not materially bias this particular library. The property columns published in
`PROGRESS.md` stand.


## Stereochemistry: the defect is confirmed, and the fix is forward-looking only

Checked directly on the first 2,000 library molecules:

- SMILES containing an explicit stereo marker (`@`): **0 of 2,000**.
- `CalcNumAtomStereoCenters` and `CalcNumUnspecifiedAtomStereoCenters` return the identical
  value for **2,000 of 2,000** molecules. In this RDKit build the former counts *potential*
  centres, and every potential centre in the library is unspecified.
- 372 of 2,000 molecules (18.6%) carry at least one potential centre; the maximum is 6.

So the published `n_stereo` column never counted defined stereochemistry -- it counted
stereocentres that the generator had flattened. The `isomericSmiles=True` fix cannot
retroactively restore that information: the existing 5,000 SMILES carry none. It applies to
future generation runs. Any stereochemical claim about the current library is unsupported and
should be withdrawn rather than recomputed.
