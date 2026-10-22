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


---

## Addendum: defects found in `pipeline.py` after the first pass

`pipeline.py` was not read during the first repair pass. It was, and it partly undid three of
the four fixes, because the repairs were made in the library modules and the CLI never picked
them up. Found and fixed 2026-10-22, each verified by running the pipeline end to end.

| # | Defect | Fix |
|---|---|---|
| 7 | **The reproducibility fix did not actually work for the BRICS backend.** RDKit's `BRICS.BRICSBuild` shuffles its reagent and seed lists with `random.shuffle(..., random=random.random)` -- the *module-level* RNG -- and exposes no seed argument. Moving `BRICSGenerator` to an instance RNG left the build order nondeterministic: two runs at `--seed 7` shared only **30 of 40** molecules. | `BRICSGenerator.__init__` now also seeds the module RNG. Verified: same seed gives byte-identical output, a different seed gives 0/40 overlap. |
| 8 | **`--seed` was ignored on the RNN path.** `CharRNNGenerator.load()` took no seed, so every loaded prior sampled at the default 42 regardless of the flag. | `load(path, device=None, seed=42)`, and `pipeline.py` passes `args.seed`. |
| 9 | **Properties described a structure the output file did not contain.** `compute_properties` standardizes internally, but the CSV's `smile` column, the Murcko scaffold, the novelty fingerprints and the LigPrep `.smi` all came from the raw generated SMILES. For a tautomer-canonicalized molecule those are different structures. | The pipeline standardizes once, explicitly, and everything downstream uses that molecule. Output now carries `smile` (standardized), `smile_as_generated`, and `standardized_changed`. |
| 10 | **`n_stereo_unassigned` never reached the output.** The column was added to `compute_properties` but not to the pipeline's row dict, so the stereochemistry fix was invisible to anyone reading the CSV. | Emitted. |

Two further defects, not in the original audit, surfaced while verifying:

- **The script would not import under `PYTHONSAFEPATH=1`** (or `python -P`), because Python then
  does not prepend the script's own directory to `sys.path` and the sibling imports fail. Fixed by
  inserting `_HERE` explicitly. This is why nothing worked on the first run attempt.
- **A BRICS-only run required PyTorch.** `generators.py` imported torch at module level, so
  fragment recombination -- which needs nothing but RDKit -- could not run without a deep-learning
  stack. The torch import is now guarded, the `nn.Module` subclass is built on first use, and
  `@torch.no_grad()` is replaced by a shim that is a pass-through when torch is absent.
  `CharRNNGenerator` raises a clear, actionable ImportError instead.

The overselling docstring in `pipeline.py` ("generation is the cheap part, the value is in
discrimination") was also removed, and replaced with a note that the SAscore ceiling in `PROFILES`
is a chosen number, which is what makes any claim about SAscore's behaviour on the resulting
library circular.

### Verification run

```
python pipeline.py generate --backend brics --frags <400 library SMILES> --decompose \
       --target 40 --oversample 20 --profile balanced --seed 7 -o out.csv
```
Exit 0. 40 molecules, 0 structural alerts in the top 10, standardization changed 0 of 40
structures (BRICS output is neutral and single-component, so only tautomers can move, and none
did in this sample). `n_stereo` equals `n_stereo_unassigned` for every molecule, which is the
same flattened-stereochemistry signature reported for the published library.
