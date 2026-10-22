# -*- coding: utf-8 -*-
"""
Score the target-class corpus with the frozen five-metric panel.

Each unique compound is standardized and scored ONCE, then joined back to its
class rows -- a compound active against two classes contributes to both, but is
not scored twice.

Must run in the `rascore` environment: RAscore's published model only loads
under xgboost 1.7 (see scorer_panel.md).
"""
import os
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))
for p in (HERE, os.path.join(HERE, "vendor")):
    if p not in sys.path:
        sys.path.insert(0, p)

import engine  # noqa: E402


def main():
    cor = pd.read_csv(os.path.join(DATA, "chembl_targetclass_corpus.csv"))
    uniq = cor.drop_duplicates("molecule_chembl_id")[["molecule_chembl_id", "smiles"]]
    print(f"corpus rows {len(cor)}, unique compounds {len(uniq)}", flush=True)

    eng = engine.Engine()
    print("metrics:", eng.metrics, flush=True)
    assert len(eng.metrics) == 5, f"panel incomplete: {eng.metrics}"

    t0 = time.time()
    scored = eng.score_frame(uniq, smiles_col="smiles", standardize=True, log_every=1000)
    print(f"scoring took {time.time()-t0:.0f}s", flush=True)

    keep = ["molecule_chembl_id", "std_smiles", "inchikey"] + eng.metrics + \
           [c for c in scored.columns if c.startswith("std_") and c != "std_smiles"]
    out = cor.merge(scored[keep], on="molecule_chembl_id", how="left")
    path = os.path.join(DATA, "corpus_scored.csv")
    out.to_csv(path, index=False)
    print(f"written -> {path}  rows={len(out)}", flush=True)

    for m in eng.metrics:
        print(f"  {m:12s} non-null {out[m].notna().sum():6d}/{len(out)}", flush=True)
    print("\nstandardization flags:", flush=True)
    for c in ("std_fragment_stripped", "std_charge_changed", "std_tautomer_changed",
              "std_tautomer_skipped"):
        if c in out:
            print(f"  {c:26s} {int(out[c].fillna(False).sum()):6d}", flush=True)


if __name__ == "__main__":
    main()
