# Review of findings — what goes in the paper, what does not

Written 2026-10-23, before the manuscript rewrite. Every number below was read back from the
saved tables, not recalled. The test applied to each finding: *does a working scientist change
what they do because of this?* If not, it is cut regardless of how much effort it took.


## The short version

One result explains everything else we found, and it arrived last. Reorganising around it cuts
roughly half the accumulated material — including the per-target-class calibration that was the
plan's centrepiece. That is the right call, and the reasoning is below.


## KEEP — the paper


### K1. Synthesizability scores track resemblance to catalogue chemistry, not whether a molecule was made

Six sets, one panel, ordered by what is actually known about each:

| Set | n | What is known | SAscore | RAscore | BR-SAScore | RAscore < 0.5 |
|---|---|---|---|---|---|---|
| ZINC purchasable | 3,581 | purchasable, so made | 3.084 | 0.980 | 4.032 | 1.3% |
| Enamine make-on-demand | 1,001 | 2-3 steps by construction | 3.505 | 0.940 | 4.508 | 8.0% |
| **THQ library** | 10,000 | **1 step, demonstrated** | 4.847 | **0.044** | 8.392 | **99.1%** |
| ChEMBL kinase actives | 3,012 | made, steps unknown | 2.757 | 0.952 | 3.969 | 15.9% |
| **BRICS recombination** | 5,000 | **never made** | 3.152 | **0.894** | 4.515 | **16.1%** |
| GDB-17 cages | 3,581 | never made | 5.730 | 0.096 | 10.000 | 94.8% |

The inversion runs both ways. A library demonstrated to be one step from catalogue building
blocks scores like enumerated strained cages; a library that has never been made scores like
purchasable compounds. Fragment recombination of ChEMBL actives produces catalogue-*looking*
molecules, which is exactly what these scores reward.

This single result explains, without further assumption, why no metric and no calibrated cutoff
flagged the one compound a synthetic chemist actually rejected: that compound is built from
ChEMBL fragments and therefore looks like catalogue chemistry.


### K2. The external validation at library scale

RAscore at its published 0.5 rejects 9,907 of 10,000 one-step THQs (99.07%); median 0.044, 99th
percentile 0.486, only 93 compounds above the cutoff. SAscore at 6.0 rejects none. The original
pipeline's own `sa_max=5.0` gate would have rejected 2,360 (23.6%).

This is someone else's data, prospective, with makeability demonstrated experimentally
(73 of 80 feasibility reactions >10% crude yield, 51 of 80 >30%). It answers the one reviewer
question the ChEMBL corpus cannot: *is the chemistry genuinely hard, or is the score wrong?*


### K3. The fairness check — this is what makes K2 publishable rather than a complaint

On RAscore's **own** published test set (19,935 AiZynthFinder-labelled ChEMBL compounds) the
vendored model gives AUC 0.9530, and at 0.5 retains 96.5% of solved compounds with 73.0%
specificity. Youden-optimal on that set is 0.8495, i.e. 0.5 is *permissive* there.

So 0.5 is defensible on the chemistry RAscore was built for. The THQ failure is an
applicability-domain failure, not a badly chosen number. Stating this pre-empts the obvious
rebuttal and is simply the honest reading.


## CUT — things that add length without changing a decision


### C1. The known-vs-generated ROC analysis. Cut entirely.

AUCs: SAscore 0.673, RAscore 0.636, SCScore 0.555, BR-SAScore 0.545, NP-likeness 0.535. We
already caveated that the generated arm is 'no evidence of existence', not 'proven unmakeable'.
K1 now shows the generated arm scores *like purchasable compounds* — so a weak AUC here is the
same phenomenon seen through a blunter instrument. Keeping both invites a reviewer to ask which
one we believe.


### C2. NP-likeness. Drop from every claim.

It is a natural-product-likeness descriptor, not a synthesizability metric. It has no published
synthesizability cutoff and no makeability interpretation. Reporting it alongside four
synthesizability scores implies a comparison that does not exist. Keep it in the released tool,
drop it from the paper's panel.


### C3. Standardization validation, molgen repairs, property deltas. Methods and repo, not results.

The salt/parent merge evidence (730 pairs, 0/730 identical before, 724/730 after) and the four
generator fixes are correctness work. They belong in a methods paragraph and the repository.
Re-scoring the JAK library both ways changed 13 of 14 property columns not at all, so there is
no result to report — only a statement that the published numbers stand.


### C4. The artifact-signature material. Already retracted; keep it retracted and brief.

One paragraph in the limitations saying the earlier framing was tested against 3,198 approved
drugs and failed (signature fires on 199 of them, 53 rated EASY, including marketed drugs made
at tonne scale). Do not reconstruct the argument to then demolish it.


## DEMOTE — the hard call


### D1. Per-target-class calibration is not the contribution. It is the remedy that does not work.

This is the uncomfortable one, because it is the largest single piece of work here: 26,947
compounds, nine classes, four figures, a calibration table, and a CLI built around it.

The finding it rests on is real and survives: at published cutoffs the false-reject rate on
chemistry known to exist is a property of the target class, spanning 32.7x on SAscore (9.48%
antibacterial vs 0.29% nuclear receptor) and 3.2x on RAscore (41.64% vs 12.84%).

But two things undercut it as a *solution*, and both are our own measurements:

1. **Ten of 36 class thresholds are pinned to the metric's own scale limit.** SCScore in 5 of 9
   classes, RAscore in 4 of 9, BR-SAScore in 1. BR-SAScore's antibacterial threshold is 10.000 —
   the top of its scale — which admits 100% of the class. RAscore's is 0.005, admitting 95.0%.
   A threshold that accepts everything is not an operating point; it is the metric declining to
   discriminate. Calling it 'calibrated' dresses up a failure.
2. **Calibration assumes the metric measures the right quantity and only the threshold is off.**
   K1 says the metric measures catalogue resemblance. Re-tuning a threshold cannot fix that, and
   our own data shows it does not: the kinase-calibrated cutoffs reject 94.7% (SAscore), 100%
   (BR-SAScore) and 80.2% (RAscore) of the one-step THQ library.

**Recommendation:** keep the per-class false-reject measurement as one figure showing that the
problem is not uniform across chemistry, and keep the demonstration that calibration fails on a
novel chemotype. Retire the calibration table as a *deliverable* — shipping thresholds we have
just shown to be inadequate would be advising the community to do something we know does not
work. The tool keeps the measurement function and loses the prescription.

A side benefit: this removes the antibacterial-permeability confound from the critical path. That
confound could not be refuted with the data we have, and it was the strongest attack available
against the old framing.


### D2. The docking selection inversion. Motivation, not result.

Correlations are modest and the two scoring functions disagree about size: better XP GScore goes
with *lower* MW (rho +0.414) while better MM-GBSA goes with *higher* MW (rho -0.275). Against
synthesizability, the significant terms are XP vs SAscore -0.264, vs BR-SAScore -0.286, vs
RAscore +0.159; MM-GBSA vs SCScore -0.232, vs RAscore +0.227, vs BR-SAScore -0.201. SAscore is
not significant against MM-GBSA, SCScore not against XP.

That is n=196, one target, one docking program, and |rho| below 0.3. It is a real observation
and it explains how `cand_004995` reached a chemist — it had the best MM-GBSA of 126 while
ranking last of 5,000 on the pipeline's own composite. But K1 explains the *rejection* better
than the docking inversion does. Use it in the introduction as the case history; do not ask a
correlation of 0.23 to carry a claim.


## What the paper becomes

Four figures, one argument:

1. The makeability gradient (K1) — the claim, as data.
2. The one-step library against published cutoffs (K2), with the RAscore fairness panel (K3).
3. Per-class false-reject rates — the problem is not uniform, and calibration does not fix it (D1).
4. The case history: how a chemist-rejected compound passed every filter (D2), as a supporting panel.

Everything else is supplementary or repository.

**What we do not have, and must say so:** route-search ground truth of our own, and any set of
compounds that chemists *attempted and failed* to make. The second does not exist publicly — it
sits inside pharma. Without it, every 'never made' set in K1 (BRICS, GDB-17) means 'no evidence
of existence', never 'proven unmakeable', and the paper must use that wording throughout.
