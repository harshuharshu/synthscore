# What Synthesizability Cutoffs Discard: Benchmarking Fast Scores and Retrosynthesis Search on Chemistry of Known Route Length

*Draft for* Journal of Chemical Information and Modeling. *All numbers in this draft were read back
from the saved result tables; none are quoted from memory.*

---

## Abstract

Fast synthesizability scores are used as gates in essentially every generative and virtual-screening
workflow, and each is applied at a published global cutoff. We ask what those cutoffs actually
measure. Scoring six reference sets with a frozen four-metric panel, ordered by what is independently
known about each set, inverts the expected relationship: a 10,000-member tetrahydroquinoline library
accessible in a single multicomponent step — 57 members of which were synthesized and assayed —
is rejected by RAscore at its published cutoff in 99.07% of cases, while a 5,000-member BRICS
recombination library that has never been made scores like purchasable compounds. The consistent reading of both directions is that these scores rank molecules by resemblance to historical catalogue chemistry rather than by synthetic accessibility. Retrosynthesis
search does not correct this. Running AiZynthFinder on 175 molecules spanning the same sets, complete
routes were found for 40% of never-made BRICS molecules against 1 of 50 one-step tetrahydroquinolines
(Fisher exact *p* = 3.4 × 10⁻⁵), and for only 48 of 125 molecules known to have been made. Because
RAscore is trained on AiZynthFinder labels, its agreement with that search (AUC 0.909) is evidence it
reproduces its teacher, blind spots included — not evidence either is right. We show the failure is
not a badly chosen threshold: on RAscore's own published test set the 0.5 cutoff retains 96.5% of
route-solved compounds at AUC 0.953. It is an applicability-domain failure of the shared USPTO
reaction corpus. We release the six labelled reference sets as a benchmark that discriminates
makeability from catalogue resemblance, which the field's current easy/hard benchmark cannot.

**Keywords:** synthetic accessibility; synthesizability prediction; retrosynthesis; applicability
domain; de novo design; benchmark datasets

---

## 1. Introduction

Every de novo design campaign must decide which generated molecules are worth making. A generative run
produces millions of candidates and a real retrosynthetic route search costs seconds to minutes per
molecule, so that decision is almost always delegated to a fast synthesizability score applied at a
published global cutoff.

Four such scores are in common use, and they reach their verdicts in fundamentally different ways.
**SAscore** (1) sums fragment-contribution terms learned from how often each substructure appears in
PubChem, then adds penalties for ring complexity, stereocentres and molecular size; it returns a value
from 1 to 10 where lower means easier, and its authors suggest 6 as a dividing line. **SCScore** (3)
is a neural network trained on 12 million USPTO reactions under a single constraint — that a product
be at least as complex as its reactants — and returns 1 to 5 with no published cutoff. **RAscore** (2)
takes a different route again: it is a binary classifier trained to predict whether AiZynthFinder
(13, 14) can find a complete route for a molecule, so it is explicitly a fast surrogate for route
search, returning a probability from 0 to 1 with 0.5 as the published decision boundary.
**BR-SAScore** (4) revises SAscore's fragment term to account for purchasable building blocks and
known reactions, and like SAscore returns a lower-is-easier scale. Two of the four therefore encode
the frequency of substructures in historical compound collections, and two encode what a
reaction-template model can plan from a historical reaction corpus. That distinction becomes the
subject of this paper.

That the global cutoffs applied to these scores are chosen rather than derived is already established.
Li et al. (5) name the problem directly, motivating their multiclass step-count predictor by noting
that "predetermined thresholds and data set imbalances present challenges" for binary classification;
Skoraczyński et al. (6) select thresholds by ROC analysis against AiZynthFinder outcomes; Yu et al.
(7) resort to the Youden index because SCScore has no natural separating value. We take the
arbitrariness of the cutoffs as given and cited, not as a finding.

The question we ask instead is what the scores measure when they are applied. Answering it requires
reference sets about which something is independently known — not sets labelled by another score. The
field's conventional benchmark, introduced with SYBA (8), contrasts purchasable ZINC compounds (9)
against enumerated GDB-17 molecules (10). Those two classes differ simultaneously in *two* respects:
whether the molecules have been made, and whether they resemble catalogue chemistry. A score that
captures only the second would appear to succeed on that benchmark. Distinguishing the two requires
sets where the two properties disagree.

Two such sets now exist. Reaction-based make-on-demand libraries are makeable by construction: their
members are defined as the products of a specified reaction applied to catalogue building blocks, so
route length is known before any scoring (16, 17). Conversely, fragment-recombination libraries built
by BRICS decomposition (11) of known drugs are catalogue-*looking* by construction while never having
been made. Placing both alongside the conventional benchmark populates the cells it omits.

---

## 2. Methods

### 2.1 Scorer panel

Four synthesizability metrics — SAscore (1), SCScore (3), RAscore (2) and BR-SAScore (4) — were each
obtained from the authors' own published implementation and weights, loaded through a common interface
and smoke-tested in the analysis environment. Nothing is cited that we could not run. Direction and
published threshold are recorded per metric so downstream code never has to infer which way a scale
runs.

Two further metrics were considered and are not part of the panel. SYBA (8) was attempted and abandoned
on evidence: its published fragment-contribution table resolves to a version-control pointer rather than
data, so the published model could not be reconstructed. NP-likeness is a natural-product-likeness
descriptor rather than a synthesizability predictor — it has no published synthesizability cutoff and no
makeability interpretation — so although the released tool computes it, it supports no claim in this
paper and is not reported in any table.

All molecules pass through a shared standardization step before scoring, implemented with RDKit (19):
largest-fragment selection, neutralization and tautomer canonicalization. The step was validated on
salt-and-parent drug pairs, which it merges correctly. Residual non-merging pairs were examined
individually and are covalent prodrugs, which must not merge.

### 2.2 Reference sets

| Set | n | What is independently known | Source |
|---|---|---|---|
| ZINC purchasable | 3,581 | purchasable, therefore made | SYBA benchmark, "easy" arm (8, 9) |
| Enamine make-on-demand | 1,001 | 2–3 steps by construction | synformer repository (16) |
| THQ library | 10,000 | one aza-Diels–Alder step, demonstrated | Papidocha et al. (17) |
| THQ made + assayed | 57 | synthesized, CDK2 IC50 measured | Papidocha et al. (17) |
| ChEMBL kinase actives | 3,012 | made; step count unknown | this work, from ChEMBL (12) |
| BRICS recombination | 5,000 | no evidence of existence | this work (11) |
| GDB-17 enumerated | 3,581 | no evidence of existence | SYBA benchmark, "hard" arm (8, 10) |

Throughout, sets with no record of synthesis are described as having **no evidence of existence**, never
as proven unmakeable. No public dataset of compounds that chemists attempted and failed to synthesize
exists; such records remain internal to industry. Their absence is this field's missing negative, and it
constrains what any benchmark — including ours — can claim.

The THQ library is the external validation set and is not ours. Papidocha et al. (17) enumerated 74,000
tetrahydroquinolines from 207 aldehydes × 361 amines plus norbornene, accessed by a single
multicomponent aza-Diels–Alder reaction, and demonstrated the chemistry experimentally: 73 of 80
feasibility reactions gave >10% crude yield and 51 of 80 >30%. They screened a 10,000-member subset and
synthesized and assayed 57 compounds. We use their deposited data unchanged.

### 2.3 Route search

AiZynthFinder (13, 14) with the published USPTO expansion policy, USPTO filter policy and ZINC stock,
run on 175 molecules — 25 sampled with a fixed seed from each of seven cells. Reported per molecule:
whether a complete route to stock compounds was found, the length of the best route, and search time.
Targets already in stock are solved in zero steps and are not route successes; one of 175 fell in that
category and is reported separately.

### 2.4 What we do not assess

Assumption diagnostics, residual analyses and alternative-estimator comparisons were not performed.
Route search was run once per molecule at default settings; search stochasticity was not characterized.

---

## 3. Results

### 3.1 Score rank inverts against what is known about the molecules

![Makeability gradient]({{artifact:art_a9d267fa-166b-4392-9a9a-103b288b2659}})

**Figure 1.** Six reference sets scored with the same panel, ordered by what is independently known
about each. Blue, demonstrably made; orange, no evidence of existence.

| Set | n | SAscore | SCScore | RAscore | BR-SAScore | RAscore < 0.5 |
|---|---|---|---|---|---|---|
| ZINC purchasable | 3,581 | 3.084 | 3.072 | 0.980 | 4.032 | 1.3% |
| Enamine make-on-demand | 1,001 | 3.505 | 4.062 | 0.940 | 4.508 | 8.0% |
| **THQ library (1 step)** | 10,000 | 4.847 | 4.346 | **0.044** | 8.392 | **99.1%** |
| ChEMBL kinase actives | 3,012 | 2.757 | 4.295 | 0.952 | 3.969 | 15.9% |
| **BRICS (never made)** | 5,000 | 3.152 | 3.899 | **0.894** | 4.515 | **16.1%** |
| GDB-17 enumerated | 3,581 | 5.730 | 3.289 | 0.096 | 10.000 | 94.8% |

The inversion runs in both directions. A library whose members are one documented step from catalogue
building blocks scores like enumerated strained polycyclic cages; a library that has never been made
scores like purchasable compounds. Both observations are consistent with a single reading: these scores
rank by resemblance to historical catalogue chemistry. BRICS recombination of fragments taken from known
drugs produces molecules that look like the catalogue, and is rewarded accordingly.

### 3.2 The external validation at library scale

![One-step library against published cutoffs]({{artifact:art_0a9d7f21-f5b8-4d39-99e0-a67877c6b7d4}})

**Figure 2.** The 10,000-member one-step library against published and class-calibrated cutoffs.

RAscore at its published 0.5 rejects 9,907 of the 10,000 one-step tetrahydroquinolines (99.07%); the
median is 0.044 and the 99th percentile 0.486, so the entire library sits below the cutoff and only 93
compounds clear it. SAscore at 6.0 rejects none (median 4.847), although a cutoff of 5.0 — a gate we had
ourselves used in earlier work — would reject 2,360 (23.6%). Restricting to the 57 compounds that were
actually synthesized and assayed gives the same picture: RAscore rejects 56 of 57.

### 3.3 Retrosynthesis search reproduces the inversion rather than correcting it

![Route search validation]({{artifact:art_d6d45f41-2178-494d-9899-12f677c9cd87}})

**Figure 3.** AiZynthFinder on 175 molecules spanning the same reference sets.

| Set | Known status | Routes found | Median steps |
|---|---|---|---|
| ZINC purchasable | made | 88% (22/25) | 1 |
| ChEMBL kinase actives | made | 56% (14/25) | 3 |
| Enamine make-on-demand | made, 2–3 steps | 44% (11/25) | 4 |
| **BRICS recombination** | **never made** | **40% (10/25)** | 4 |
| THQ made + assayed | made, 1 step | 4% (1/25) | 5 |
| THQ library | made, 1 step | 0% (0/25) | 5 |
| GDB-17 enumerated | never made | 0% (0/25) | 5 |

Complete routes were found for 10 of 25 never-made BRICS molecules against 1 of 50 one-step
tetrahydroquinolines (Fisher exact *p* = 3.4 × 10⁻⁵, odds ratio 32.7). Across the 125 molecules known to
have been made, 77 (61.6%) were declared unroutable.

Two controls bound the interpretation. Only one of 175 targets was already in stock, so the 88% solve
rate on purchasable ZINC compounds reflects genuine route finding rather than catalogue lookup. And
GDB-17 cages solve at 0%, confirming the search does reject chemistry that is genuinely hard. The search
is working; it is blind in a specific place.

The sharpest statement available from this experiment is that route search cannot distinguish
tetrahydroquinolines made in one step (0–4% solved) from GDB-17 strained virtual cages (0% solved).
These two sets sit at opposite extremes of real makeability and are indistinguishable to it.

### 3.4 The agreement between score and route search is not validation

RAscore and AiZynthFinder return the same verdict for 78.9% of the 175 molecules, and RAscore predicts
route solvability at AUC 0.909. This is expected by construction: RAscore is a surrogate trained on
AiZynthFinder labels (2). Their agreement is evidence that the student reproduces the teacher —
including the teacher's blind spots — not that either tracks makeability.

The practical consequence for a reader choosing a filter is that "use a better score" and "use
retrosynthesis instead" are the same answer. Both rest on the USPTO reaction corpus, and both fail on
precisely the chemistry that reaction-based make-on-demand libraries are built from.

### 3.5 The cutoff is defensible on the chemistry it was built for

Before criticising a published threshold we tested it on its authors' own benchmark. On RAscore's
published test set — 19,935 AiZynthFinder-labelled ChEMBL compounds, 14,994 route-solved and 4,941 not —
the vendored model achieves AUC 0.9530, and at 0.5 retains 96.47% of solved compounds with 73.0%
specificity on unsolved. The Youden-optimal threshold on that set is 0.8495, so 0.5 is in fact
permissive there.

The 0.5 cutoff is therefore well chosen for ChEMBL-like chemistry. The 99.07% rejection of a one-step
library is an applicability-domain failure, not a threshold-selection error. This distinction is the
difference between a useful finding and a complaint.

### 3.6 Per-class calibration is not the remedy

![Per-class score distributions]({{artifact:art_789edfc9-2c69-413c-8ba0-a22754d268c1}})

**Figure 4.** False-reject rates at published cutoffs across nine target classes, 27,513 ChEMBL
compounds.

At published cutoffs the false-reject rate on chemistry known to exist is a property of the target class,
not of the molecule: SAscore at 6.0 rejects 9.48% of whole-cell antibacterials against 0.29% of
nuclear-receptor ligands, and RAscore at 0.5 rejects 41.64% against 12.84% for oxidoreductase ligands.

We computed per-class operating points at a declared false-reject rate and report them as a measurement,
not as a recommendation, for two reasons that are themselves measurements.

First, 9 of 36 class thresholds are pinned at the metric's own scale limit. BR-SAScore's whole-cell
antibacterial threshold is 10.000 — the top of its scale — which admits 100% of the class; RAscore's is
0.005, admitting 95.0%; SCScore's kinase threshold is 5.000, admitting 95.2%. A threshold that accepts
everything is not an operating point. It is the metric declining to discriminate, and calling it
"calibrated" presents a failure as a solution.

Second, calibration assumes the metric measures the right quantity and only the threshold is misplaced.
Sections 3.1 and 3.3 indicate the quantity itself is wrong, and our own data confirm re-tuning does not
help: cutoffs calibrated on historical ChEMBL kinase actives reject 94.67% (SAscore), 100.0%
(BR-SAScore) and 80.16% (RAscore) of the one-step library. Historical chemistry is sp2-rich; a novel
sp3-rich chemotype falls outside it whatever threshold is applied.

We therefore retire the calibration table as a shipped deliverable. Distributing thresholds we have
shown to be inadequate would advise the community to do something we know does not work.

---

## 4. A case history

![Docking selection inversion]({{artifact:art_bc69cc83-ff2c-47b0-bdbb-58a0f6dc8c03}})

**Figure 5.** Docking rank against synthesizability in a JAK1 campaign (n = 196 XP, 126 MM-GBSA).

The work began with a conventional campaign: a BRICS (11) and RNN generator seeded with ChEMBL JAK
inhibitors, a medicinal-chemistry filter stack including an SAscore gate, Glide XP docking and MM-GBSA
rescoring against PDB 6N7A. A synthetic chemist reviewed the top-ranked compounds and declined them,
either as unmakeable or as requiring an impractically long sequence.

No panel metric had flagged them. The compound the chemist singled out ranked 4,995th of 5,000 on the
pipeline's own composite desirability while carrying the best MM-GBSA interaction energy of 126 scored
compounds. Selection pressure and synthesizability pressure pointed in opposite directions: better
MM-GBSA correlates with harder-to-make chemistry (SCScore ρ = +0.23, *p* = 0.009; RAscore ρ = +0.23,
*p* = 0.011; BR-SAScore ρ = +0.20, *p* = 0.024, oriented so positive means harder).

These correlations are modest, come from one target and one docking program, and the two scoring
functions disagree about molecular size. We present this as the case history that motivated the
investigation, not as a result that carries weight. Section 3.1 provides the better explanation for why
nothing flagged the compound: built from ChEMBL fragments, it resembles catalogue chemistry.

---

## 5. Discussion

The results support a corpus-level rather than a parameter-level claim. Fast synthesizability scores and
retrosynthesis search both rest on the USPTO reaction corpus. A reaction class absent from that corpus
is invisible to the whole stack, however few steps it takes at the bench. The aza-Diels–Alder
multicomponent reaction used by Papidocha et al. (17) is such a class, which is why a library
demonstrated to be one step from commercial building blocks is rejected by a fast score at 99.07% and by
route search at 96–100%.

This reframes what these tools can be asked to do. They are reliable interpolators within documented
reaction space and silent outside it, and they give no signal distinguishing the two conditions. A
molecule scored as hard may be hard, or may be outside the corpus; the score cannot tell a user which,
and neither can a route search that shares the corpus.

The direction this suggests is not a better score. It is to calibrate on chemistry one can actually
make — which reaction-based enumerated libraries supply for free, since their makeability is known by
construction (16). A group with a defined reaction set can measure its own score distribution and set
its own operating point without needing a general-purpose threshold at all.

Papidocha et al. (17) support this reading from the opposite direction. They note that prospective
applications of cofolding affinity models have been "largely been restricted to screening commercial
catalogs", which are built from a handful of robust transformations and so span a synthetically narrow,
sp2-rich space. That is the field's workaround for the synthesizability problem, stated in its own
words. Our data quantify what the workaround costs.

### 5.1 Relation to prior work

Li et al. (5) name the threshold problem and solve it by predicting step counts from reaction-dataset
labels; our contribution is complementary — we measure what the existing cutoffs discard rather than
building a replacement. Parrot et al. (15) established ChEMBL as a comparison set for these scores, so
the corpus idea is not ours; the stratification and the error-rate measurement are. Skoraczyński et al.
(6) derive global thresholds by ROC against AiZynthFinder on a small routed set. Yadunandan and Ghosh
(18) report per-target variation in the SAscore of generated molecules without drawing a calibration
conclusion. We claim none of these as novel and cite them as the foundation they are.

### 5.2 An approach we tested and rejected

Before arriving at the reference-set design above, we tested a structural proxy: the hypothesis that
fused-polycyclic systems chained together are a recognisable signature of unmakeable fragment-recombined
molecules, and could serve as a filter in their own right. The hypothesis fails two independent checks
and we report it so that others do not repeat it.

First, the test was circular on our own library. That library's generator applies an SAscore ceiling
upstream, so its maximum SAscore sits exactly at that ceiling by construction, and no statement about
SAscore's behaviour on that library can be evidence of anything.

Second, the signature does not discriminate. Applied to 3,198 approved drugs it fires on 199 of them, 53
of which SAscore rates as easy — including marketed drugs manufactured at scale, for which the lenient
verdict is correct. A signature that flags successful commercial products is not measuring
manufacturability.

We also note that stereochemistry cannot be analysed in that library: none of its structures carried
defined stereochemistry, so stereocentre counts computed from it describe flattened potential centres
rather than real ones.

### 5.3 Limitations

The external validation rests on a single tetrahydroquinoline core. Ten thousand substituent variations
of 207 aldehydes, 361 amines and norbornene are not ten thousand independent molecules, and we claim
only that an entire reaction-accessible chemotype is discarded — not that ten thousand independent
compounds are misclassified.

A failure to find a route is not evidence that a molecule cannot be made. AiZynthFinder found routes for
only 44% of Enamine make-on-demand compounds and 56% of ChEMBL kinase actives, both known-makeable, so
"solved" is a soft label and we do not use route-search failure as proof of anything.

No set used here is a true negative. Every set we describe as never made means only that no record of
synthesis exists. Compounds that chemists attempted and failed to make would be the real negative
control, and those records are not public.

The whole-cell antibacterial result carries a confound we cannot refute. That bin is defined at the
level of the organism rather than a protein, so its distinctive score distribution may reflect
permeability chemistry rather than synthetic complexity. With per-class calibration demoted (Section
3.6) the result no longer lies on this paper's critical path, and it should not be read as a
synthesizability finding.

The route search was run in one configuration: a single run per molecule, default settings, one
expansion policy. Neither search stochasticity nor the choice of policy was explored.

---

## 6. Availability

Code, the frozen scorer panel, the six labelled reference sets and all result tables are released under
MIT (derived ChEMBL data under CC BY-SA). Third-party model weights are re-fetched by a provided script
rather than redistributed. The primary deliverable is the benchmark: the six sets with their
provenance, plus the scoring entry point, so that a new metric can be tested for whether it tracks
makeability or catalogue resemblance.

## Supporting Information

The following files are available free of charge:

- Per-molecule scores for all six reference sets, with standardized SMILES (CSV)
- Per-molecule route-search output for the 175 sampled molecules, including route length, search time and precursor counts (CSV)
- Rejection rates at published and class-calibrated cutoffs for the one-step library (CSV)
- RAscore threshold check on its own published test set (CSV)
- Per-target-class score distributions and false-reject rates across 27,513 ChEMBL compounds (CSV)
- Docking and synthesizability correlations for the JAK1 case history (CSV)

---

## References

(1) Ertl, P.; Schuffenhauer, A. Estimation of synthetic accessibility score of drug-like molecules based on molecular complexity and fragment contributions. *J. Cheminform.* **2009**, *1*, 8. DOI: 10.1186/1758-2946-1-8.

(2) Thakkar, A.; Chadimová, V.; Bjerrum, E. J.; Engkvist, O.; Reymond, J. L. Retrosynthetic accessibility score (RAscore) - rapid machine learned synthesizability classification from AI driven retrosynthetic planning. *Chem. Sci.* **2021**, *12*, 3339–3349. DOI: 10.1039/d0sc05401a.

(3) Coley, C. W.; Rogers, L.; Green, W. H.; Jensen, K. F. SCScore: Synthetic Complexity Learned from a Reaction Corpus. *J. Chem. Inf. Model.* **2018**, *58*, 252–261. DOI: 10.1021/acs.jcim.7b00622.

(4) Chen, S.; Jung, Y. Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore. *J. Cheminform.* **2024**, *16*, 83. DOI: 10.1186/s13321-024-00879-0.

(5) Li, X.; Walsh, R.; Abbas, W.; Pascual-Diaz, S.; Hand, C.; Garland, R.; Khan, F. M.; Mohan Das, N.; Desai, V.; AbouZleikha, M.; et al. Multiclass Synthetic Accessibility Prediction. *J. Chem. Inf. Model.* **2025**, *65*, 1155–1165. DOI: 10.1021/acs.jcim.4c01663.

(6) Skoraczyński, G.; Kitlas, M.; Miasojedow, B.; Gambin, A. Critical assessment of synthetic accessibility scores in computer-assisted synthesis planning. *J. Cheminform.* **2023**, *15*, 6. DOI: 10.1186/s13321-023-00678-z.

(7) Yu, J.; Wang, J.; Zhao, H.; Gao, J.; Kang, Y.; Cao, D.; Wang, Z.; Hou, T. Organic Compound Synthetic Accessibility Prediction Based on the Graph Attention Mechanism. *J. Chem. Inf. Model.* **2022**, *62*, 2973–2986. DOI: 10.1021/acs.jcim.2c00038.

(8) Voršilák, M.; Kolář, M.; Čmelo, I.; Svozil, D. SYBA: Bayesian estimation of synthetic accessibility of organic compounds. *J. Cheminform.* **2020**, *12*, 35. DOI: 10.1186/s13321-020-00439-2.

(9) Irwin, J. J.; Tang, K. G.; Young, J.; Dandarchuluun, C.; Wong, B. R.; Khurelbaatar, M.; Moroz, Y. S.; Mayfield, J.; Sayle, R. A. ZINC20-A Free Ultralarge-Scale Chemical Database for Ligand Discovery. *J. Chem. Inf. Model.* **2020**, *60*, 6065–6073. DOI: 10.1021/acs.jcim.0c00675.

(10) Ruddigkeit, L.; van Deursen, R.; Blum, L. C.; Reymond, J. L. Enumeration of 166 billion organic small molecules in the chemical universe database GDB-17. *J. Chem. Inf. Model.* **2012**, *52*, 2864–2875. DOI: 10.1021/ci300415d.

(11) Degen, J.; Wegscheid-Gerlach, C.; Zaliani, A.; Rarey, M. On the art of compiling and using 'drug-like' chemical fragment spaces. *ChemMedChem* **2008**, *3*, 1503–1507. DOI: 10.1002/cmdc.200800178.

(12) Zdrazil, B.; Felix, E.; Hunter, F.; Manners, E. J.; Blackshaw, J.; Corbett, S.; de Veij, M.; Ioannidis, H.; Lopez, D. M.; Mosquera, J. F.; et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. *Nucleic Acids Res.* **2024**, *52*, D1180–D1192. DOI: 10.1093/nar/gkad1004.

(13) Genheden, S.; Thakkar, A.; Chadimová, V.; Reymond, J. L.; Engkvist, O.; Bjerrum, E. AiZynthFinder: a fast, robust and flexible open-source software for retrosynthetic planning. *J. Cheminform.* **2020**, *12*, 70. DOI: 10.1186/s13321-020-00472-1.

(14) Saigiridharan, L.; Hassen, A. K.; Lai, H.; Torren-Peraire, P.; Engkvist, O.; Genheden, S. AiZynthFinder 4.0: developments based on learnings from 3 years of industrial application. *J. Cheminform.* **2024**, *16*, 57. DOI: 10.1186/s13321-024-00860-x.

(15) Parrot, M.; Tajmouati, H.; da Silva, V. B. R.; Atwood, B. R.; Fourcade, R.; Gaston-Mathé, Y.; Do Huu, N.; Perron, Q. Integrating synthetic accessibility with AI-based generative drug design. *J. Cheminform.* **2023**, *15*, 83. DOI: 10.1186/s13321-023-00742-8.

(16) Gao, W.; Luo, S.; Coley, C. W. Generative Artificial Intelligence for Navigating Synthesizable Chemical Space. *arXiv* **2024**, arXiv:2410.03494.

(17) Papidocha, S. M.; Ogura, Y.; Ren, J.; Saykin, M.; Green, S.; Bernales, V.; Aspuru-Guzik, A. Discovery of Nanomolar CDK2 Inhibitors Using Boltz-2. *ChemRxiv* **2026**. DOI: 10.26434/chemrxiv.15010755/v1.

(18) Yadunandan, R.; Ghosh, N. ReACT-Drug: Reaction-Template Guided Reinforcement Learning for de novo Drug Design. *arXiv* **2025**, arXiv:2512.20958.

(19) RDKit: Open-Source Cheminformatics. https://www.rdkit.org (accessed 2026-10-23).
