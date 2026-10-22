# Audit of the Synthesizability Finding

Independent re-check of `PROGRESS.md` sections 3, 4b and 4c against the code and data as committed.

Date: 2026-10-22. Scope: `code/synthesizability/artifact_analysis.py`, `code/synthesizability/synth_predict.py`, `code/molgen/medchem.py`, `code/molgen/generators.py`, `data/*.csv`.


## Summary

The arithmetic reproduces exactly. The *interpretation* does not survive a control.
Two defects are severe enough to stop publication in the current framing; four are fixable bugs.


## What reproduces (no defects found)

| Claim | Verified |
|---|---|
| 5,000 valid molecules | 5,000 |
| SAscore EASY / MODERATE / HARD | 3,666 / 1,334 / 0 |
| Artifact signature flagged | 702 (14.0%) |
| Failure population (artifact AND EASY) | 329 (6.6%) |
| `cand_004995` | SA 3.450, 6 rings, 2 fused bicyclic units, EASY |
| Stored `sa_score` column trustworthy | recomputed, max abs diff 0.005 (rounding) |

The union-find in `artifact_analysis.py` merges ring systems sharing **one** atom, which conflates
spiro junctions with true fusion. I re-ran with the correct rule (>=2 shared atoms). **Zero molecules**
change class. The shortcut did not bias anything here.


## Defect 1 (critical): the headline is circular

`PROGRESS.md` 4b: *"SAscore never once returns HARD on this library."* True, and uninformative.

`medchem.py` `PROFILES["balanced"]` sets `sa_max=5.0`, and `passes_windows()` rejects any molecule above it.
Measured maximum SAscore in `jak_library.csv` is **exactly 5.0000**, with zero molecules above.
The library is the *output* of an SAscore filter. HARD is defined as >6.0, so HARD count is
structurally forced to zero before any analysis runs.

SAscore is used twice in the pipeline: as a hard gate (`passes_windows`) and as a ranking term
(`desirability()`, `_d_min(sa, 2.5, 6.0)`). Evaluating SAscore on a library SAscore selected and ranked
cannot support a claim about SAscore's behaviour.

Control: on 3,198 approved small-molecule drugs from ChEMBL (`max_phase=4`), unfiltered, SAscore returns
HARD **111 times** and reaches **8.98**. SAscore does flag hard molecules. It was simply never allowed to here.


## Defect 2 (critical): the artifact signature does not mark unmakeability

The signature is `>=2 fused ring systems each containing >=2 rings`. Applied unchanged to approved drugs:

| | JAK library | Approved drugs |
|---|---|---|
| Artifact signature rate | 14.0% | **6.2%** (199 drugs) |
| ...of those, SAscore EASY | 329 | **53** |
| SAscore EASY rate | 73.3% | 61.0% |
| SAscore HARD count | 0 | 111 |
| SAscore median | 3.15 | 3.17 |
| SAscore max | 5.00 | 8.98 |

Marketed drugs carrying the signature and rated EASY include brexpiprazole (SA 2.39), belumosudil (2.42),
domperidone (2.44), fruquintinib (2.45), telmisartan (2.52), ziprasidone (2.61) and capmatinib (2.67).
These are synthesized at tonne scale. For them SAscore's EASY verdict is **correct**.

So "SAscore fails to penalize route-implausible structures" cannot rest on this signature: the signature
selects a common medicinal-chemistry motif, not a synthetic impossibility. The chemist rejected
`cand_004995` specifically; the signature generalises that single rejection to 702 molecules without
evidence that the generalisation holds.

What *is* real: the signature is **2.3x enriched** in the generated library versus approved drugs
(14.0% vs 6.2%). BRICS recombination does over-produce polycyclic chaining relative to real drug space.
That is a defensible, control-anchored observation. It is a statement about output distribution,
not a makeability classifier.


## Defect 3 (critical for the benchmark): labels derive from the tested signature

`benchmark_input.csv` assigns 50 molecules `artifact_pred_unmakeable` using the Defect-2 signature.
Since 6.2% of approved drugs satisfy the same criterion, an unknown fraction of those 50 are makeable.
A high AiZynthFinder solve rate in the artifact arm would then be a labelling error, not a result.
The benchmark cannot answer its question as constructed.


## Defect 4 (bug): Tier 3 will return empty verdicts

`synth_predict.py::run_aizynth` keys results on the SMILES string `aizynthcli` echoes back
(`results[smi] = ...`), then merges with `if smi in res` against the **raw CSV string**.
aizynthcli canonicalizes its targets, so the two strings differ for many inputs and the merge
silently drops them, leaving `aizynth_solved` None for every row.

`PROGRESS.md` 4e attributes this risk to the `orient="table"` JSON schema. The schema may also be wrong,
but string-identity merging is an independent defect that would produce the same empty output.
Fix: canonicalize both sides with `Chem.MolToSmiles(Chem.MolFromSmiles(s))` before keying, or merge on row order.

Related, same function: all 100 molecules go to one `subprocess.run(..., timeout=None)` call. One stalled
route search hangs the entire run with no partial output. Chunk the call and set a per-chunk timeout.


## Defect 5 (bug): the analysis script as committed does not run

`artifact_analysis.py` line 19 reads `pd.read_csv("Files/jak_library.csv")`. No `Files/` directory exists;
the data is in `data/`. `PROGRESS.md` section 7 documents a different path again (`../../data/`).
Anyone reproducing from the repo hits an immediate FileNotFoundError.


## Defect 6: the four correctness items, confirmed in source

All four from `PROGRESS.md` section 5 check out line by line:

1. **Stereochemistry flattened.** `generators.py` uses `isomericSmiles=False` in both the BRICS path and the
   CharRNN sampling path, while `medchem.py::compute_properties` computes and stores `n_stereo`. Stereocenters
   are counted on structures that no longer carry them.
2. **RNN not reproducible.** No `torch.manual_seed` anywhere. `CharRNNGenerator.fit` calls `random.shuffle(seqs)`
   on the global RNG; the class has no seeded RNG of its own (only `BRICSGenerator` has `self.rng`).
3. **No standardization.** `compute_properties` runs on the directly parsed molecule. No salt stripping,
   no neutralization, no tautomer canonicalization. cLogP, TPSA, HBD and HBA are therefore computed on an
   arbitrary protonation and salt state, and the filter thresholds are applied to those numbers.
4. **Overselling docstrings, verbatim.** `medchem.py`: "This is where candidate quality is actually decided".
   `generators.py`: "You own the weights -> IP-clean" and "the whole pipeline is reproducible and yours".


## What survives, and the non-circular experiment

Two findings survive the audit:

1. The JAK docking campaign result (hinge-motif filter did not enrich; ligand-efficiency correction
   equalised queue and control) is untouched by any defect above.
2. BRICS recombination is 2.3x enriched for chained polycyclic systems relative to approved drugs,
   and the SAscore gate at 5.0 passes all of it.

The question worth testing is non-circular and still open:

> **At equal SAscore and equal heavy-atom count, do generated molecules route less often than approved drugs?**

Design: matched pairs. For each generated molecule, draw an approved drug with the same SAscore (+/- 0.1)
and similar size. Run AiZynthFinder on both arms blind. The drug arm is the positive control with a known
answer — every one of them was actually synthesized — so it also validates the route-search setup itself.
Neither arm's label comes from the signature under test. Either outcome is publishable: a gap quantifies
how far generative output sits from makeable space at matched complexity; no gap is a clean negative
result about SAscore-matched comparisons that the field currently assumes without testing.
