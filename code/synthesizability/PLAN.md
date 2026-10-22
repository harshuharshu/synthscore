# Synthesizability Prediction Tool - Runnable Plan

A general-purpose, target-agnostic tool to predict whether generated molecules
can actually be made, *before* they reach a synthetic chemist. Built because
your JAK campaign proved the point: molecules with great docking scores and good
SA scores were declared unmakeable at the bench.

## Why this design (the science)

The ASAP benchmark (Skoraczynski et al., J. Cheminform. 2023) tested whether the
common fast scores predict real retrosynthesis outcomes. Key findings that shape
this tool:

- **SAscore and RAscore** discriminate makeable from unmakeable reasonably well
  (AUC ~0.81 against AiZynthFinder). **SCScore and SYBA are ~20 points worse.**
- **RAscore is essentially a fast approximation of AiZynthFinder's verdict** -
  useful, but it inherits AiZynthFinder's blind spots.
- The **actual ground truth** is a retrosynthesis route search (AiZynthFinder):
  does a route to purchasable building blocks exist, and how long is it.

Proven on your own data: the tris-fused molecule your chemist rejected
(cand_004995) scores SAscore 3.47 = "easy" - *better than imatinib*. SAscore
alone will NOT catch it. Only the retrosynthesis tier will. That single result
is the justification for a tiered tool rather than a one-score filter.

## Architecture: three tiers, graceful degradation

| Tier | Tool | Speed | What it tells you | Install cost |
|------|------|-------|-------------------|--------------|
| 1 | SAscore | instant | complexity proxy (weak on fused/chained systems) | none (bundled with RDKit) |
| 2 | RAscore | ~ms | prob. a route exists (best fast score) | moderate (separate env) |
| 3 | AiZynthFinder | ~30-120 s/mol | **ground truth**: solved + n_steps | higher (conda env + model download) |

Run Tier 1 on everything, Tier 2 on everything, Tier 3 only on the top-N
survivors (it is slow). The tool auto-detects which tiers are installed and uses
what it can.

## Install

### Tier 1 only (works today, no new environment)
```
pip install -r requirements-tier1.txt
python synth_predict.py --in your_molecules.csv --smiles-col smile --out scored.csv
```
Use this immediately to get a first pass. But remember: Tier 1 alone missed
cand_004995. It is a coarse pre-filter, not the answer.

### Tier 3 - AiZynthFinder (the ground truth; do this for the real analysis)

AiZynthFinder needs **Python 3.10-3.12 in its own conda env** (do NOT try to put
it in your 3.14 - it will not work). Your laptop already has 3.12 via `py -3.12`,
but AiZynthFinder is cleanest in a dedicated conda environment:

```
conda create -n aizynth python=3.11 -y
conda activate aizynth
conda install -c conda-forge rdkit -y
pip install aizynthfinder

# download the public USPTO model + ZINC stock (one time, ~1 GB):
download_public_data my_aizynth_data
# this writes config.yml pointing at the model + stock
```

Then run the retrosynthesis directly to sanity-check one molecule:
```
aizynthcli --config my_aizynth_data/config.yml --smiles "CC(=O)Oc1ccccc1C(=O)O"
```

### Tier 2 - RAscore (optional fast approximation)

RAscore is finicky (older TF/sklearn pins). Install in its own env if you want
the fast pre-filter for Tier 3:
```
conda create -n rascore python=3.10 -y
conda activate rascore
pip install rdkit
git clone https://github.com/reymond-group/RAscore
cd RAscore && pip install -e .
```

## Running the full pipeline

Because the tiers may live in different conda envs, the practical workflow is:

**Step 1 - Tier 1 + 2 (in the rascore env, or Tier 1 alone anywhere):**
```
python synth_predict.py --in ranked.csv --smiles-col smile \
       --rascore --out scored_fast.csv
```

**Step 2 - Tier 3 on survivors (in the aizynth env):**
```
python synth_predict.py --in scored_fast.csv --smiles-col smile \
       --aizynth --aizynth-config my_aizynth_data/config.yml \
       --aizynth-topn 50 --out scored_final.csv
```

The tool runs AiZynthFinder on only the top 50 by RAscore, so you spend the
expensive compute where it counts. Output adds `aizynth_solved` (True/False) and
`aizynth_steps`, and a fused `synth_tier`:
EASY / MODERATE / HARD / LIKELY_UNMAKEABLE.

## How to read the output

- **LIKELY_UNMAKEABLE** = AiZynthFinder found no route (or RAscore < 0.25). These
  are your cand_004995-type molecules. Drop them before they reach a chemist.
- **EASY / MODERATE** = a route exists; n_steps tells you how painful. Advance these.
- Sort your library by `synth_tier` then by your binding/desirability score to get
  a shortlist that is both potent-predicted AND makeable.

## For the paper (if you go that route)

You now have the ingredients for a genuine methods finding:
1. Generate library (done - JAK, BRICS).
2. Score with SAscore, RAscore, AiZynthFinder (this tool).
3. Show SAscore misclassifies chemist-rejected molecules as easy (done for
   cand_004995 - reproduce across the set).
4. Show the retrosynthesis tier catches them.
5. Optional strong version: have 1-2 chemists label a defined set for makeability
   and route length; benchmark each score against that human ground truth. That
   turns it from "we applied tools" into "we measured which tool matches reality."

The honest, publishable message: complexity != infeasibility, fast scores miss
chained/fused artifacts, and retrosynthesis-in-the-loop is needed for generative
pipelines. Your own molecules are the case study.

## Extending to other targets

The tool is fully target-agnostic - it only reads SMILES. Point it at any
generated library (RIPK1, or anything from your generator). Nothing here is
JAK-specific. The same scored.csv workflow applies to every target you test.
```
