# -*- coding: utf-8 -*-
"""
synthscore - target-class calibrated synthesizability scoring.

Scores molecules with up to five published metrics and judges each against a
cutoff calibrated for the target class you name, rather than the metric's
published global cutoff.

Why that matters: at its published cutoff, SAscore discards 9.5% of known
antibacterial chemistry but only 0.29% of known nuclear-receptor chemistry, and
RAscore discards 41.6% vs 12.8%. The cutoff that is right for one class is wrong
for another. Calibration here means: the cutoff that discards a declared fraction
(default 5%) of compounds in that class which were demonstrably synthesized,
because they carry a measured ChEMBL activity.

Examples
--------
  python cli.py --list-classes
  python cli.py --in molecules.csv --target-class kinase --out scored.csv
  python cli.py --in molecules.csv --smiles-col SMILES --target-class antibacterial_wholecell \\
                --false-reject 0.10 --panel sascore,rascore --out scored.csv
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))
for p in (HERE, os.path.join(HERE, "vendor")):
    if p not in sys.path:
        sys.path.insert(0, p)

CALIBRATION = os.path.join(DATA, "calibration_thresholds.csv")

DIRECTION = {"sascore": "lower_easier", "nplikeness": "higher_easier",
             "scscore": "lower_easier", "rascore": "higher_easier",
             "brsascore": "lower_easier"}
PUBLISHED = {"sascore": 6.0, "rascore": 0.5}
# NP-likeness is not a synthesizability predictor; it is never used as a gate.
GATING = [m for m in DIRECTION if m != "nplikeness"]


def load_calibration(path=CALIBRATION):
    if not os.path.exists(path):
        raise SystemExit(f"calibration table not found: {path}\n"
                         "Run build_corpus.py then score_corpus.py to regenerate it.")
    return pd.read_csv(path)


def cutoff_for(cal, target_class, metric, false_reject):
    col = f"thr_fr{int(round(false_reject * 100))}"
    if col not in cal.columns:
        raise SystemExit(f"no calibration at a {false_reject:.0%} false-reject rate; "
                         f"available: {[c for c in cal.columns if c.startswith('thr_fr')]}")
    row = cal[(cal.target_class == target_class) & (cal.metric == metric)]
    return None if row.empty else float(row[col].iloc[0])


def rejects(value, cutoff, metric):
    if value is None or cutoff is None or pd.isna(value):
        return None
    return bool(value > cutoff) if DIRECTION[metric] == "lower_easier" else bool(value < cutoff)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="synthscore", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="infile", help="input CSV")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--id-col", default=None, help="optional identifier column to carry through")
    ap.add_argument("--target-class", help="calibration class; see --list-classes")
    ap.add_argument("--panel", default="all",
                    help="comma-separated metric subset, or 'all' (default)")
    ap.add_argument("--false-reject", type=float, default=0.05,
                    help="tolerated fraction of real chemistry discarded (default 0.05)")
    ap.add_argument("--no-standardize", action="store_true",
                    help="skip salt stripping / neutralization / tautomer canonicalization")
    ap.add_argument("--out", help="output CSV")
    ap.add_argument("--list-classes", action="store_true")
    a = ap.parse_args(argv)

    cal = load_calibration()
    if a.list_classes:
        for cls in sorted(cal.target_class.unique()):
            sub = cal[cal.target_class == cls]
            n = int(sub.n.iloc[0])
            print(f"{cls:26s} n={n:6d}  SAscore 5% cutoff "
                  f"{float(sub[sub.metric=='sascore'].thr_fr5.iloc[0]):.3f}")
        return 0

    missing = [f for f in ("infile", "target_class", "out") if not getattr(a, f)]
    if missing:
        ap.error("required unless --list-classes: " + ", ".join("--" + m.replace("_", "-")
                                                                for m in missing))
    if a.target_class not in set(cal.target_class):
        raise SystemExit(f"unknown target class {a.target_class!r}; see --list-classes")

    want = None if a.panel == "all" else [m.strip() for m in a.panel.split(",")]
    if want:
        unknown = [m for m in want if m not in DIRECTION]
        if unknown:
            raise SystemExit(f"unknown metric(s): {unknown}; known: {list(DIRECTION)}")

    import engine
    eng = engine.Engine(want=want)
    if not eng.metrics:
        raise SystemExit("no metrics could be loaded; see scorer_panel.md")
    absent = [m for m in (want or DIRECTION) if m not in eng.metrics]
    if absent:
        print(f"[warn] unavailable, continuing without: {absent}", file=sys.stderr)
        for m in absent:
            print(f"       {m}: {eng.panel.unavailable.get(m, 'not requested')}", file=sys.stderr)

    df = pd.read_csv(a.infile)
    if a.smiles_col not in df.columns:
        raise SystemExit(f"column {a.smiles_col!r} not in {a.infile}; "
                         f"columns are {list(df.columns)}")

    scored = eng.score_frame(df, smiles_col=a.smiles_col,
                             standardize=not a.no_standardize, log_every=500)

    gating = [m for m in eng.metrics if m in GATING]
    cutoffs = {m: cutoff_for(cal, a.target_class, m, a.false_reject) for m in gating}
    for m in gating:
        scored[f"cutoff_{m}"] = cutoffs[m]
        scored[f"reject_{m}"] = [rejects(v, cutoffs[m], m) for v in scored[m]]
        if m in PUBLISHED:
            scored[f"reject_{m}_published"] = [rejects(v, PUBLISHED[m], m) for v in scored[m]]
    scored["reject_any"] = scored[[f"reject_{m}" for m in gating]].fillna(False).any(axis=1)

    scored.to_csv(a.out, index=False)

    n = len(scored)
    print(f"\nscored {n} molecules against {a.target_class} cutoffs "
          f"({a.false_reject:.0%} false-reject calibration)")
    for m in gating:
        r = scored[f"reject_{m}"].fillna(False)
        line = (f"  {m:11s} cutoff {cutoffs[m]:8.3f}  rejected {int(r.sum()):5d} "
                f"({100*r.mean():5.1f}%)")
        if m in PUBLISHED:
            rp = scored[f"reject_{m}_published"].fillna(False)
            line += (f"   | published {PUBLISHED[m]:g}: {int(rp.sum()):5d} "
                     f"({100*rp.mean():5.1f}%)")
        print(line)
    print(f"  rejected by any metric: {int(scored.reject_any.sum())} "
          f"({100*scored.reject_any.mean():.1f}%)")
    print(f"written -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
