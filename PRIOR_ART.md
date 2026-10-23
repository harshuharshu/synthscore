# Prior-art check: has this been published?

Searched 2026-10-22 across PubMed (6 natural-language + 4 field-tagged queries), Europe PMC
(12 field-tagged TITLE_ABS queries, 75 unique records retrieved) and one web search. OpenAlex was
skipped: no API key is configured, and keyless requests are not permitted.

**Short answer: no. The core claim is unpublished.** Two adjacent papers come close enough that they
must be cited prominently, and one framing from the earlier draft is genuinely pre-empted.


## Claim-by-claim verdict

| # | Claim | Prior art | Verdict |
|---|---|---|---|
| 1 | Synthesizability-score cutoffs have **target-class-dependent** error rates | Europe PMC `"synthetic accessibility" OR synthesizability` AND `"target class" OR "protein class" OR "target family"` returns **1 hit**, on kinase fragment-library filtering, unrelated | **Novel** |
| 2 | Measuring what a cutoff **discards** on chemistry known to be makeable | `... AND "false negative" OR "false rejection" OR discard` returns **0 hits** | **Novel** |
| 3 | **Per-class operating points** at a declared false-reject rate | No hit in any query; threshold work is global | **Novel** |
| 4 | Antibacterial chemistry is systematically penalized | 22 hits all *apply* SA scores to antibacterials; none measures score bias against them | **Novel** |
| 5 | ChEMBL activity record as a makeability label | Partly pre-empted (see A2) | **Partly new** |
| 6 | Global cutoffs are arbitrary / need optimization | Pre-empted (see A1, A3) | **Not novel - concede and cite** |
| 7 | Natural products vs synthetics score differently | Ertl & Schuffenhauer 2009, Fig. 4 | **Not novel - already withdrawn** |

## The two papers that come closest

**A1. Multiclass Synthetic Accessibility Prediction.** *J. Chem. Inf. Model.* 2025.
doi:10.1021/acs.jcim.4c01663

This is the closest statement of our *problem*. Its own abstract motivates the work by noting that
"predetermined thresholds and data set imbalances present challenges" for binary synthetic-accessibility
classification. It then solves that problem differently: it predicts the **minimum number of synthetic
steps** as a multiclass task, ensembled over stratified subsampled folds, trained and evaluated on two
reaction benchmark datasets, and introduces fuzzy metrics that tolerate near-miss predictions.

Why it does not scoop us: it builds a new predictor and its ground truth is reaction-dataset step
counts. It neither stratifies by biological target class nor measures how much real, previously
synthesized chemistry an existing score's published cutoff rejects. **We must cite it as the paper that
named the threshold problem, and state plainly that our contribution is the complementary one:
calibrating the scores people already use, rather than replacing them.**

**A2. Integrating synthetic accessibility with AI-based generative drug design.** *J. Cheminform.*
2023, 15:83. doi:10.1186/s13321-023-00742-8

Compares a commercial retrosynthesis-derived score against RAscore, SCScore and SAscore on a sample of
5,000 molecules drawn from a pre-processed ChEMBL set. So **ChEMBL-as-comparison-set is not new**. The
difference is what the set is used for: there, to justify treating one score as ground truth for the
others; here, as a per-target-class positive label with which to measure each cutoff's false-reject
rate. We should stop describing the ChEMBL corpus itself as a contribution and describe only its
**stratification and the error-rate measurement** as new.


## Threshold work that already exists (cite, do not re-derive)

**A3. Critical assessment of synthetic accessibility scores in computer-assisted synthesis planning.**
*J. Cheminform.* 2023. doi:10.1186/s13321-023-00678-z. Finds optimal thresholds by ROC analysis,
balancing sensitivity against specificity, with AiZynthFinder outcomes as the label. Two differences
matter and both favour us: the thresholds are **global**, and the routed set is small - AiZynthFinder
found routes for 22 compounds in their analysis, against our 26,947-compound positive set.

**A4. Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism (GASA).**
*J. Chem. Inf. Model.* 2022. doi:10.1021/acs.jcim.2c00038. Selects its own cutoff by Youden index
because SCScore has no natural separating threshold - direct precedent that these cutoffs are chosen,
not derived.

**A5.** A de novo design preprint (ReACT-Drug, arXiv:2512.20958) reports mean SAscore varying across six
protein targets (2.66 for M2 receptor to 3.62) in *generated* molecules. This is the closest published
**observation** of per-target variation, but it is descriptive, on generated rather than real chemistry,
and draws no calibration or error-rate conclusion. Worth one sentence as the nearest observation.


## Score-method papers to cite as the panel

| Metric | Citation |
|---|---|
| SAscore | Ertl & Schuffenhauer, *J. Cheminform.* 2009, 1:8. doi:10.1186/1758-2946-1-8 |
| SYBA | *J. Cheminform.* 2020. doi:10.1186/s13321-020-00439-2 |
| RAscore | *Chem. Sci.* 2021. doi:10.1039/d0sc05401a |
| GASA | *J. Chem. Inf. Model.* 2022. doi:10.1021/acs.jcim.2c00038 |
| BR-SAScore | *J. Cheminform.* 2024. doi:10.1186/s13321-024-00879-0 |
| RetroGNN | *J. Chem. Inf. Model.* 2022. doi:10.1021/acs.jcim.1c01476 |
| Reaction knowledge graph | *Molecules* 2022. doi:10.3390/molecules27031039 |
| Molecular complexity in GDB space | *J. Chem. Inf. Model.* 2025. doi:10.1021/acs.jcim.5c00334 |

SCScore (Coley 2018) and DFRscore (2023) are in the panel documentation and should be cited from the
originals; they were not re-retrieved in this search.


## What this means for the framing

Keep, because nothing contests them:

1. The false-reject rate of a published cutoff is a function of target class, measured at 33-fold
   (SAscore) and 3.2-fold (RAscore) spread.
2. Per-class operating points at a declared false-reject rate.
3. RAscore's published 0.5 is stricter than real chemistry tolerates in any class; SCScore has no
   usable headroom on this corpus. Neither observation appears anywhere in the retrieved literature.
4. Antibacterial chemistry is the worst-served class.

Drop or demote:

1. "These cutoffs are arbitrary" as a *finding* - it is already stated in A1, A3 and A4. Use it as
   motivation with citations, not as a result.
2. "A ChEMBL ground-truth corpus" as a contribution - A2 did this first. The stratification is ours;
   the corpus idea is not.


## What a reviewer will probably say

1. *"Is the antibacterial effect just whole-cell permeability chemistry, not synthetic complexity?"*
   This is the strongest attack and we cannot currently refute it - the bin is organism-level by
   construction. Pre-empt it in the discussion; do not wait to be asked.
2. *"Without retrosynthesis ground truth, how do you know the scores are wrong rather than the
   chemistry being genuinely hard?"* Partly answered - every corpus compound was made - but a routed
   subset would close it. A1 and A3 both had route-level labels and we do not.
3. *"Class imbalance / target redundancy."* Answered: ~3,000 per class by design, 150-compound ceiling
   per target, seeded shuffle, 513 multi-class compounds disclosed.
4. *"Why not just use a better score?"* Answer with R4: at matched size none of the five separates
   known from generated chemistry above AUC 0.673, so better-score-shopping does not solve it.


## Suggested venue

*Journal of Cheminformatics* is the natural home: A2, A3, SYBA and BR-SAScore all appeared there, so
the audience is exactly the one that sets these cutoffs, and it is fully open access with code and
data expectations this work already meets. *J. Chem. Inf. Model.* is the alternative and carries A1,
A4 and RetroGNN - but publishing a calibration of A1's problem statement in A1's own venue invites
the reviewer who wrote A1.
