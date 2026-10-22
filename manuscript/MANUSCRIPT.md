# Target-class calibration of synthesizability scores

*Draft skeleton. Numbers are read from the saved tables, not retyped. No conclusions section
is written beyond what the tables support.*


## Working title

One cutoff, nine answers: synthesizability scores are miscalibrated across target classes,
and antibacterial chemistry pays the most.


## Abstract (skeleton)

Synthesizability scores are deployed in generative and virtual-screening pipelines as hard
filters, each with a single global cutoff recommended by its authors. We ask what those cutoffs
cost on chemistry that is known to be makeable. Using 26,947 ChEMBL compounds that each
carry a measured activity -- so each was demonstrably synthesized -- spread over 9
target classes, we score five published metrics (SAscore, NP-likeness, SCScore, RAscore,
BR-SAScore) and measure the fraction of real chemistry each global cutoff discards, per class.
The rate is not a property of the molecule but of its target class: SAscore at 6.0 discards
9.48% of whole-cell antibacterial compounds and 0.29% of nuclear-receptor
compounds (33-fold), and RAscore at 0.5 discards 41.64% and 12.84%
(3.2-fold). We publish per-class operating points calibrated to a declared
false-reject rate and a tool that applies them. We also report two negative results: at matched
molecular size no metric separates known chemistry from generative output better than
AUC 0.67, and none of the class-calibrated cutoffs reproduces a synthetic
chemist's rejection of a specific generated molecule.


## Positioning against prior art

Three prior results bound what is new here, and the framing must concede each explicitly.

1. **Ertl & Schuffenhauer (2009)** introduced SAscore and in the same paper showed that natural
   products, bioactive molecules and catalogue compounds occupy different SAscore distributions.
   *That distributions differ by compound provenance is therefore not a finding of ours.* Our
   claim is narrower and operational: the error rate induced by the recommended cutoff differs
   by target class, and we quantify it.
2. **CoPriNet and related purchasability work** compared SAscore, SCScore, SYBA and RAscore on
   purchasable versus non-purchasable compounds. That is a different label (availability) from
   ours (demonstrated synthesis), and it is not resolved by target class.
3. **Skoraczynski et al. (2023), J. Cheminform. 15:6** benchmarked fast scores against a
   retrosynthesis planner and concluded the scores mostly discriminate feasible from infeasible
   molecules but must be chosen carefully for use inside planning algorithms. Our work supplies
   the missing operating points rather than re-running that comparison; we have no
   route-search ground truth (see Limitations).

The contribution is therefore: (a) a non-circular ground-truth corpus labelled by experimental
record rather than by any score, stratified by target class; (b) measured per-class false-reject
rates at published cutoffs; (c) calibrated per-class operating points and a tool that uses them;
(d) two honest negative results.


## Methods

### Ground truth

Every compound with a measured ChEMBL activity was synthesized and assayed. We take corpus
membership as evidence of makeability that derives from no synthesizability score, which is what
breaks the circularity that otherwise afflicts this comparison: a library pre-filtered by a score
cannot be used to evaluate that score.

We assembled 27,513 (compound, class) rows covering 26,947 unique compounds
across 9 classes: eight ChEMBL protein classes (kinase, protease, GPCR family A,
ion channel, nuclear receptor, transporter, epigenetic, oxidoreductase) plus a whole-cell
antibacterial bin, which ChEMBL models as organism-type rather than protein-class targets and
which was collected separately from 377 bacterial targets. 513 compounds are active
against more than one class and appear in each.

**No potency or assay-type filter was applied.** The label is *was made*, not *was potent*. This
matters materially: only 42.7% of corpus rows carry a pChEMBL value at all, and
just 0.6% of the antibacterial rows, because MIC readouts sit under functional
assays. A conventional pChEMBL >= 6 filter would have erased the antibacterial class entirely.

Targets were visited in a seeded shuffle with a 150-compound ceiling per target, so that no
heavily screened target dominates its class, and each class was capped at ~3,000 compounds.
Heavy-atom medians span 25-32 across classes (Figure 1), so the class
differences reported below are not trivially size-driven.


### Structure standardization

All molecules pass through RDKit Cleanup, largest-fragment selection, charge neutralization and
canonical tautomer generation, with stereochemistry retained. This is load-bearing, not cosmetic:
across 730 approved-drug salt/parent pairs, 0 scored identically before
standardization (maximum divergence 2.59 SAscore units and 5.35
BR-SAScore units for the same substance), and 724 of 730 collapse to one
structure and one score after it. The 19 covalent ester prodrugs in the same comparison
correctly remain distinct (0 merged). Regression tests enforce both behaviours,
including the case that distinguishes an acetic-acid salt from an acetate ester.


### Metric panel

Five metrics, all obtained and run locally; provenance, the compatibility fixes applied to
vendored code, and the metrics that could not be obtained are documented in `scorer_panel.md`.
NP-likeness is reported as a covariate and never used as a gate. SYBA, GASA, FSscore and
DFRscore could not be obtained in this environment and are named rather than silently omitted.


### Calibration

There is no proven-unmakeable negative set, so calibration is one-class. For each (class, metric)
the operating point is the quantile of known-makeable chemistry that yields a declared
false-reject rate (we publish 5% and 10%). This bounds how much real chemistry a cutoff discards;
it makes no claim about how much unmakeable chemistry it admits.


## Results

### R1. Published cutoffs have class-dependent error rates (Figure 3)

| Class | SAscore @ 6.0 | RAscore @ 0.5 |
|---|---|---|
| antibacterial wholecell | 9.48% | 41.64% |
| nuclear receptor | 0.29% | 30.54% |
| ion channel | 2.97% | 28.51% |
| protease | 0.6% | 28.38% |
| transporter | 3.93% | 22.1% |
| gpcr familyA | 2.54% | 18.51% |
| kinase | 0.76% | 15.9% |
| epigenetic | 2.75% | 14.03% |
| oxidoreductase | 1.36% | 12.84% |

The spread is 33-fold for SAscore and 3.2-fold for RAscore. Antibacterial
whole-cell chemistry is the worst-served class on both, and is rated hardest on three of the four
synthesizability metrics (Figure 2); SCScore's hardest-rated class is GPCR family A instead.

### R2. A single cutoff cannot be rescued by recalibrating it globally

Pooling all classes and recalibrating SAscore to a 5% target gives 5.178. That number
still discards 16.0% of antibacterial chemistry while discarding only
1.63% of kinase chemistry. The defect is not that the published number is wrong;
it is that one number cannot serve nine classes.

### R3. Calibrated per-class operating points (Table 2, `calibration_thresholds.csv`)

At a 5% false-reject rate the required cutoffs span:

- SAscore 4.528 (kinase) to 7.151 (antibacterial_wholecell) -- straddling the published 6.0
- BR-SAScore 7.069 to 10.000, the upper value being the scale ceiling: BR-SAScore saturates on antibacterial chemistry
- RAscore 0.001 to 0.098 -- two orders of magnitude below its published 0.5
- SCScore 4.878 to 5.000, i.e. pinned at the top of its scale, leaving no usable headroom

The RAscore and SCScore rows are findings in their own right. RAscore's published 0.5 is not
merely class-inappropriate; it is far stricter than real chemistry tolerates in any class.
SCScore, on this corpus, cannot be operated as a filter at a tolerable false-reject rate at all.

### R4. Negative result: no metric separates known from generated chemistry (Figure 4a)

Size-matched comparison of corpus compounds against 5,000 BRICS-generated molecules:

| Metric | AUC [95% CI] |
|---|---|
| SAscore | 0.673 [0.667-0.681] |
| RAscore | 0.636 [0.628-0.643] |
| SCScore | 0.555 [0.547-0.563] |
| BR-SAScore | 0.545 [0.537-0.552] |
| NP-likeness | 0.535 [0.526-0.542] |

The generated arm carries no evidence of existence, but is not *proven* unmakeable, so a low AUC
does not by itself convict a metric -- the generated set may genuinely be about as makeable as
the corpus. What the result does rule out is using any of these metrics as a confident
real-versus-generated discriminator at matched molecular size.

### R5. Negative result: calibration does not reproduce expert rejection (Figure 5)

A synthetic chemist rejected one generated molecule (`cand_004995`) by inspection. It passes
every kinase-calibrated cutoff: SAscore 3.485 (cutoff 4.528), SCScore 4.918 (5.000), RAscore
0.160 (0.098), BR-SAScore 4.914 (7.069). Only RAscore's published 0.5 rejects it, and that same
cutoff discards 15.9% of known kinase chemistry. Setting any cutoff exactly at the molecule's own
score costs 7.4% of real kinase chemistry on RAscore (224 compounds), 14.2% on SCScore, 16.9% on
SAscore and 29.5% on BR-SAScore. The structural signature that motivated the original hypothesis
does shift every metric in the harder direction (median SAscore +0.489, SCScore +0.510,
BR-SAScore +0.550, RAscore -0.276 for flagged versus unflagged), but the shift is small relative
to within-class spread and is not a separating criterion. We report this rather than tuning it away.


## Figures

1. `corpus_class_profile.png` -- corpus composition: compounds retained per class and
   heavy-atom distributions, establishing that classes are size-comparable.
2. `class_distributions.png` -- per-class distribution of all five metrics with published cutoffs marked.
3. `cost_of_global_threshold.png` -- known chemistry discarded per class at SAscore 6.0 and RAscore 0.5.
4. `roc_panel.png` -- (a) known versus generated discrimination; (b) per-class SAscore cutoff at a 5% false-reject rate against the published 6.0.
5. `jak_recalibrated.png` -- `cand_004995` and the generated library positioned against known kinase chemistry on four metrics.


## Limitations

- **No proven-unmakeable negatives.** Calibration is one-class; false-accept rates are not estimated.
- **No route-search ground truth.** AiZynthFinder was unobtainable in this environment, so no
  metric here is validated against an actual retrosynthesis search. This is the single most
  valuable addition for a revision.
- **Corpus membership proves synthesis, not ease.** A ChEMBL compound may have taken 14 steps.
- **Class assignment is ChEMBL's**, and is coarse; multi-class compounds appear in each class.
- **Panel incomplete.** SYBA, GASA, FSscore and DFRscore are absent for environment reasons, not
  scientific ones. SYBA in particular is a fragment-based classifier whose behaviour might differ.
- **The antibacterial bin is whole-cell, not target-based**, so permeability-driven chemistry is
  confounded with synthetic complexity. This is a plausible partial explanation for R1 and should
  be stated in the discussion rather than resolved by assertion.


## Withdrawn from the earlier draft

Two claims from the prior write-up do not survive and must not appear:

1. *"SAscore never returns HARD on this library"* -- circular. The library was generated behind an
   SAscore <= 5.0 gate and HARD is defined as > 6.0, so the count was forced to zero before any
   analysis ran. On unfiltered approved drugs SAscore returns HARD 111 times and reaches 8.98.
2. *The structural artifact signature marks unmakeability* -- it fires on 6.2% of approved drugs
   including brexpiprazole, telmisartan, ziprasidone and capmatinib, which are made at scale. What
   survives is enrichment only: 14.0% in the generated library versus 6.2% in approved drugs.
Stereochemistry claims about the existing JAK library are also withdrawn: no library SMILES
carries an explicit stereo marker, so the reported stereocentre counts were counts of flattened
potential centres. The generator is fixed going forward; the old numbers are not recomputed.
