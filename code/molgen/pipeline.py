#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pipeline.py - generate -> filter -> ranked shortlist.

Candidates flow cheap->expensive through structural alerts, property windows,
synthesizability and novelty, then rank by a Derringer desirability composite.

Note on the synthesizability filter: PROFILES in medchem.py gates on an SAscore
ceiling, and that ceiling is a chosen number, not a derived one. Any statement
about SAscore's behaviour on a library produced by this pipeline is therefore
circular (see AUDIT.md, Defect 1). Target-class calibrated cutoffs, measured
against chemistry known to have been synthesized, are in code/synthesizability/.

GENERATION BACKENDS
  --backend brics   : fragment recombination from a fragment library (no training)
  --backend rnn     : sample a trained char-RNN SMILES prior (--prior model.pt)
  --backend both    : pool BRICS + RNN candidates

TRAINING THE RNN PRIOR (do this once, on your GPU)
  python pipeline.py train --train-smiles chembl_druglike.smi --prior my_prior.pt \
         --epochs 25 --hidden 512 --layers 3

GENERATING
  # BRICS:
  python pipeline.py generate --backend brics --frags fragments.csv --smiles-col SMILES \
         --target 5000 --profile balanced -o shortlist.csv

  # RNN prior:
  python pipeline.py generate --backend rnn --prior my_prior.pt --temperature 1.0 \
         --target 5000 --profile balanced -o shortlist.csv

  # both, with novelty scoring against your known actives:
  python pipeline.py generate --backend both --frags fragments.csv --prior my_prior.pt \
         --known actives.smi --target 8000 --profile balanced -o shortlist.csv

OUTPUT
  Ranked CSV (best first) with full property panel, desirability + components,
  structural-alert names, Murcko scaffold, and (if --known given) max Tanimoto to
  known actives + a novelty flag. A LigPrep-ready .smi is written alongside.
"""

import argparse
import os
import sys

# Python does not prepend the script's own directory to sys.path when
# PYTHONSAFEPATH is set (or under `python -P`), so sibling imports below would
# fail. Add it explicitly.
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import pandas as pd
from rdkit import Chem, RDLogger

import medchem as mc  # noqa: F401  (also puts code/synthesizability on sys.path)
from generators import BRICSGenerator, CharRNNGenerator
from standardize import Standardizer

RDLogger.DisableLog("rdApp.*")

# Standardize once per candidate, here, so that the structure written to the
# output file is the same structure the properties were computed on. Calling
# mc.compute_properties(standardize=True) would standardize internally and
# silently return properties for a structure the caller never sees.
_STD = Standardizer()


# ---------------------------------------------------------------------------
def read_smiles_file(path, smiles_col="SMILES"):
    if path.lower().endswith((".smi", ".txt")):
        with open(path) as fh:
            return [ln.strip().split()[0] for ln in fh if ln.strip()]
    df = pd.read_csv(path)
    if smiles_col not in df.columns:
        cols = list(df.columns)
        sys.exit("ERROR: column " + repr(smiles_col) + " not in " + path + ". Columns: " + str(cols))
    return df[smiles_col].dropna().astype(str).tolist()


# ---------------------------------------------------------------------------
def cmd_train(args):
    smiles = read_smiles_file(args.train_smiles, args.smiles_col)
    print(f"Training char-RNN prior on {len(smiles)} SMILES "
          f"(epochs={args.epochs}, hidden={args.hidden}, layers={args.layers})...")
    g = CharRNNGenerator()
    print(f"  device: {g.device}")
    g.fit(smiles, epochs=args.epochs, batch_size=args.batch_size, lr=args.lr,
          hidden=args.hidden, layers=args.layers, max_len=args.max_len)
    g.save(args.prior)
    print(f"Saved prior -> {args.prior}")


# ---------------------------------------------------------------------------
def get_candidate_stream(args):
    """Return an iterator of SMILES from the chosen backend(s)."""
    streams = []
    if args.backend in ("brics", "both"):
        frags = read_smiles_file(args.frags, args.smiles_col)
        if args.decompose:
            from rdkit.Chem import BRICS
            fset = set()
            for s in frags:
                m = Chem.MolFromSmiles(s)
                if m:
                    fset.update(BRICS.BRICSDecompose(m))
            frags = sorted(fset)
            print(f"  decomposed to {len(frags)} BRICS fragments")
        bg = BRICSGenerator(frags, n_sample=args.n_sample, seed=args.seed)
        # oversample: many candidates die in the gauntlet
        streams.append(bg.sample(args.target * args.oversample))
    if args.backend in ("rnn", "both"):
        g = CharRNNGenerator.load(args.prior, seed=args.seed)
        print(f"  RNN prior loaded on {g.device} (seed={args.seed})")
        streams.append(g.sample(args.target * args.oversample,
                                temperature=args.temperature))

    def chained():
        for st in streams:
            for smi in st:
                yield smi
    return chained()


# ---------------------------------------------------------------------------
def cmd_generate(args):
    # reference fingerprints for novelty, if known actives provided
    ref_fps, known_scaffolds = [], set()
    if args.known:
        known = read_smiles_file(args.known, args.smiles_col)
        for s in known:
            m = Chem.MolFromSmiles(s)
            if m:
                ref_fps.append(mc.morgan_fp(m))
                known_scaffolds.add(mc.murcko_scaffold(m))
        print(f"Loaded {len(ref_fps)} known actives for novelty scoring.")

    print(f"Generating (backend={args.backend}, target={args.target}, "
          f"profile={args.profile})...")
    stream = get_candidate_stream(args)

    seen = set()
    rows = []
    n_seen = n_invalid = n_dup = n_failed = 0

    for smi in stream:
        n_seen += 1
        if smi in seen:
            n_dup += 1
            continue
        seen.add(smi)
        raw = Chem.MolFromSmiles(smi)
        if raw is None:
            n_invalid += 1
            continue
        mol, _std_info = _STD.standardize(raw)
        if mol is None:
            n_invalid += 1
            continue
        std_smi = Chem.MolToSmiles(mol, isomericSmiles=True)

        props = mc.compute_properties(mol, standardize=False)
        alerts = mc.structural_alerts(mol)
        ok, reasons = mc.passes_windows(props, alerts, args.profile)

        # in permissive mode we keep everything and just annotate
        if not ok and args.profile != "permissive":
            n_failed += 1
            continue

        D, comp = mc.desirability(props)
        scaf = mc.murcko_scaffold(mol)
        row = {
            "smile": std_smi,          # the structure the properties below describe
            "smile_as_generated": smi,  # pre-standardization, for provenance
            "standardized_changed": std_smi != Chem.MolToSmiles(raw, isomericSmiles=True),
            "desirability": D,
            "qed": round(props["qed"], 3),
            "mw": round(props["mw"], 1),
            "clogp": round(props["clogp"], 2),
            "tpsa": round(props["tpsa"], 1),
            "hbd": props["hbd"], "hba": props["hba"],
            "rotb": props["rotb"],
            "fsp3": round(props["fsp3"], 3),
            "n_aromatic_rings": props["n_aromatic_rings"],
            "n_stereo": props["n_stereo"],
            "n_stereo_unassigned": props["n_stereo_unassigned"],
            "sa_score": round(props["sa"], 2),
            "n_alerts": len(alerts),
            "alerts": "; ".join(alerts[:5]),
            "murcko_scaffold": scaf,
            "passes_windows": ok,
        }
        if ref_fps:
            sim = mc.max_tanimoto(mol, ref_fps)
            row["max_tanimoto_known"] = round(sim, 3)
            row["novel_scaffold"] = scaf not in known_scaffolds
            row["novelty_ok"] = sim <= args.max_similarity
            # optionally drop near-duplicates of known actives
            if args.drop_similar and sim > args.max_similarity:
                continue
        rows.append(row)

        if len(rows) >= args.target:
            break
        if len(rows) % 500 == 0 and len(rows) > 0:
            print(f"  kept {len(rows)} / target {args.target} "
                  f"(seen {n_seen}, gauntlet-rejected {n_failed})", flush=True)

    if not rows:
        sys.exit("No candidates survived. Loosen --profile or check inputs.")

    df = pd.DataFrame(rows).sort_values("desirability", ascending=False).reset_index(drop=True)
    df.insert(0, "compound_id", [f"cand_{i:06d}" for i in range(len(df))])
    df.to_csv(args.out, index=False)

    smi_path = args.out.rsplit(".", 1)[0] + ".smi"
    df[["smile", "compound_id"]].to_csv(smi_path, sep=" ", index=False, header=False)

    print(f"\n=== DONE ===")
    print(f"Seen {n_seen} | invalid {n_invalid} | dup {n_dup} | "
          f"gauntlet-rejected {n_failed} | kept {len(df)}")
    print(f"Ranked shortlist -> {args.out}")
    print(f"LigPrep-ready    -> {smi_path}")
    print(f"\nTop 10 by desirability:")
    cols = ["compound_id", "smile", "desirability", "qed", "mw", "clogp", "fsp3", "sa_score", "n_alerts"]
    print(df[cols].head(10).to_string(index=False))


# ---------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description="Med-chem-aware molecule generator.")
    sub = p.add_subparsers(dest="cmd", required=True)

    # train
    t = sub.add_parser("train", help="train the char-RNN SMILES prior")
    t.add_argument("--train-smiles", required=True, help="SMILES file/CSV to train on")
    t.add_argument("--smiles-col", default="SMILES")
    t.add_argument("--prior", required=True, help="output checkpoint path (.pt)")
    t.add_argument("--epochs", type=int, default=25)
    t.add_argument("--batch-size", type=int, default=128)
    t.add_argument("--lr", type=float, default=1e-3)
    t.add_argument("--hidden", type=int, default=512)
    t.add_argument("--layers", type=int, default=3)
    t.add_argument("--max-len", type=int, default=120)
    t.set_defaults(func=cmd_train)

    # generate
    g = sub.add_parser("generate", help="generate + filter + rank")
    g.add_argument("--backend", choices=["brics", "rnn", "both"], default="brics")
    g.add_argument("--frags", help="fragment library (for brics/both)")
    g.add_argument("--decompose", action="store_true", help="BRICS-decompose --frags first")
    g.add_argument("--prior", help="trained RNN checkpoint (for rnn/both)")
    g.add_argument("--temperature", type=float, default=1.0, help="RNN sampling temperature")
    g.add_argument("--smiles-col", default="SMILES")
    g.add_argument("--known", help="known actives SMILES for novelty scoring")
    g.add_argument("--max-similarity", type=float, default=0.7,
                   help="Tanimoto above this = too similar to a known active")
    g.add_argument("--drop-similar", action="store_true",
                   help="drop candidates above --max-similarity to known actives")
    g.add_argument("--profile", choices=["strict", "balanced", "permissive"], default="balanced")
    g.add_argument("--target", type=int, default=2000, help="size of final shortlist")
    g.add_argument("--oversample", type=int, default=10,
                   help="generate target*oversample raw candidates (most die in gauntlet)")
    g.add_argument("--n-sample", type=int, default=4, help="BRICS fragments per batch")
    g.add_argument("--seed", type=int, default=42)
    g.add_argument("-o", "--out", default="shortlist.csv")
    g.set_defaults(func=cmd_generate)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    if args.cmd == "generate" and args.backend in ("brics", "both") and not args.frags:
        sys.exit("ERROR: --frags required for backend 'brics'/'both'.")
    if args.cmd == "generate" and args.backend in ("rnn", "both") and not args.prior:
        sys.exit("ERROR: --prior required for backend 'rnn'/'both'.")
    args.func(args)
