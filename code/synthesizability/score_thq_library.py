# -*- coding: utf-8 -*-
"""
Score the external THQ validation library with the frozen five-metric panel.

Source: Papidocha et al., "Discovery of Nanomolar CDK2 Inhibitors Using Boltz-2",
chemRxiv 2026, doi:10.26434/chemrxiv.15010755/v1 (CC-BY 4.0). Library file
THQ_virtual_library_10000_seed42.csv from github.com/papidocha/CDK2-Inhibitor-Discovery.

Why this set: every member is one aza-Diels-Alder step from commercial building
blocks, and the authors demonstrated that chemistry experimentally (73 of 80
feasibility reactions gave >10% crude yield, 51 of 80 >30%). Route length is
therefore known to be ONE, which no ChEMBL-derived set can claim.

Must run in the `rascore` environment (RAscore's published model loads only
under xgboost 1.7).
"""
import os
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, "vendor")):
    if p not in sys.path:
        sys.path.insert(0, p)

import engine  # noqa: E402

IN_CSV = os.environ.get("THQ_IN", "thq_library_10k.csv")
OUT_CSV = os.environ.get("THQ_OUT", "thq_library_10k_scored.csv")


def main():
    df = pd.read_csv(IN_CSV)
    assert "Product_SMILES" in df.columns, df.columns.tolist()
    eng = engine.Engine()
    print(f"metrics: {eng.metrics}", flush=True)
    if len(eng.metrics) != 5:
        print(f"WARNING: only {len(eng.metrics)} metrics loaded: {eng.unavailable}", flush=True)

    t0 = time.time()
    out = eng.score_frame(df, smiles_col="Product_SMILES", standardize=True, log_every=500)
    out.to_csv(OUT_CSV, index=False)
    print(f"\nwritten -> {OUT_CSV}  rows={len(out)}  elapsed={time.time()-t0:.0f}s", flush=True)
    for m in eng.metrics:
        print(f"  {m:12s} non-null {out[m].notna().sum():6d}/{len(out)}", flush=True)
    for c in ("std_fragment_stripped", "std_charge_changed", "std_tautomer_changed",
              "std_tautomer_skipped"):
        if c in out:
            print(f"  {c:26s} {int(out[c].fillna(False).sum()):6d}", flush=True)


if __name__ == "__main__":
    main()
