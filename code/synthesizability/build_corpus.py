# -*- coding: utf-8 -*-
"""
Build the ChEMBL target-class ground-truth corpus.

Ground truth premise: every compound with a measured ChEMBL activity was actually
synthesized and assayed. Membership in this corpus is therefore evidence of
makeability that does not derive from any synthesizability score -- which is what
makes it a non-circular benchmark.

Design notes
------------
* No potency threshold. The label is "this molecule was made", not "this molecule
  was potent". Filtering on pchembl would restrict the corpus to optimized leads,
  which plausibly correlates with synthetic tractability and would bias the very
  thing we are measuring. pchembl_value is recorded per row so a potency-stratified
  sensitivity check remains possible.
* No assay-type filter, for the same reason. Binding vs functional readout is
  irrelevant to whether the compound exists. (It also matters practically: whole-cell
  antibacterial MIC data is functional and often carries no pchembl_value at all.)
* Targets within a class are visited in a seeded shuffle and each contributes at most
  MAX_PER_TARGET compounds, so no single heavily-screened target dominates its class.

Output: data/corpus_<class>.csv per class plus data/chembl_targetclass_corpus.csv.
"""
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

import pandas as pd

BASE = "https://www.ebi.ac.uk/chembl/api/data"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")

# Level-1 / level-2 ChEMBL protein classes. Ids verified against the live
# protein_classification endpoint, not assumed.
CLASSES = {
    "kinase":           6,
    "protease":         1028,
    "gpcr_familyA":     1020,
    "ion_channel":      1005,
    "nuclear_receptor": 1033,
    "transporter":      1006,
    "epigenetic":       836,
    "oxidoreductase":   10,
}

CAP_PER_CLASS = 3000    # compounds retained per class
MAX_PER_TARGET = 150    # ceiling per target, so one screened target cannot dominate
SHUFFLE_SEED = 0

# Whole-cell antibacterial bin. Not a protein class: ChEMBL models these as
# target_type=ORGANISM, so they are collected separately and labelled as such.
BACTERIAL_GENERA = [
    "Staphylococcus", "Streptococcus", "Escherichia", "Pseudomonas", "Klebsiella",
    "Acinetobacter", "Enterococcus", "Mycobacterium", "Salmonella", "Bacillus",
    "Haemophilus", "Neisseria", "Clostridium", "Clostridioides", "Enterobacter",
    "Proteus", "Shigella", "Vibrio", "Helicobacter", "Campylobacter", "Listeria",
    "Legionella", "Moraxella", "Serratia", "Burkholderia", "Stenotrophomonas",
    "Bacteroides", "Fusobacterium", "Corynebacterium", "Micrococcus", "Nocardia",
    "Mycobacteroides", "Mycolicibacterium",
]


def get(ep, **params):
    params.setdefault("format", "json")
    url = f"{BASE}/{ep}?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if 400 <= e.code < 500:
                raise
            last = e
            time.sleep(3 * (attempt + 1))
        except Exception as e:
            last = e
            time.sleep(3 * (attempt + 1))
    raise last


def paged(ep, key, limit=1000, **params):
    off = 0
    while True:
        j = get(ep, limit=limit, offset=off, **params)
        rows = j[key]
        if not rows:
            return
        for r in rows:
            yield r
        off += limit
        if not j["page_meta"]["next"]:
            return


def log(*a):
    print(*a, flush=True)


def build_class_ancestry():
    """Map every protein_classification_id to the set of its ancestor class ids."""
    pc = pd.DataFrame(list(paged("protein_classification", "protein_classifications")))
    parent = dict(zip(pc.protein_class_id, pc.parent_id))
    anc = {}
    for cid in pc.protein_class_id:
        chain, cur, guard = set(), cid, 0
        while cur is not None and guard < 20:
            chain.add(cur)
            cur = parent.get(cur)
            guard += 1
        anc[cid] = chain
    log(f"protein_classification rows: {len(pc)}")
    return anc


def build_target_class_map(anc):
    """target_chembl_id -> set of our class labels."""
    want = {cid: lab for lab, cid in CLASSES.items()}
    tmap = {}
    n = 0
    for tc in paged("target_component", "target_components", component_type="PROTEIN",
                    only="component_id,accession,organism,protein_classifications,targets"):
        n += 1
        labs = set()
        for pcl in tc.get("protein_classifications") or []:
            for a in anc.get(pcl["protein_classification_id"], ()):
                if a in want:
                    labs.add(want[a])
        if not labs:
            continue
        for t in tc.get("targets") or []:
            tmap.setdefault(t["target_chembl_id"], set()).update(labs)
        if n % 2000 == 0:
            log(f"  components scanned {n}, classified targets so far {len(tmap)}")
    log(f"protein components scanned: {n}; targets with a class: {len(tmap)}")
    return tmap


def fetch_target(tid, seen, label, extra=None):
    """Distinct new compounds tested against one target, capped at MAX_PER_TARGET."""
    out, got = [], 0
    for a in paged("activity", "activities", limit=1000, target_chembl_id=tid,
                   only=("molecule_chembl_id,canonical_smiles,pchembl_value,"
                         "target_chembl_id,standard_type,assay_type,assay_chembl_id")):
        mid, smi = a.get("molecule_chembl_id"), a.get("canonical_smiles")
        if not mid or not smi or mid in seen:
            continue
        seen.add(mid)
        row = dict(target_class=label, molecule_chembl_id=mid, smiles=smi,
                   pchembl_value=a.get("pchembl_value"),
                   target_chembl_id=a.get("target_chembl_id"),
                   standard_type=a.get("standard_type"),
                   assay_type=a.get("assay_type"),
                   assay_chembl_id=a.get("assay_chembl_id"))
        if extra:
            row.update(extra.get(tid, {}))
        out.append(row)
        got += 1
        if got >= MAX_PER_TARGET:
            break
    return out


def fetch_class(label, target_ids, extra=None):
    """Walk a class's targets in seeded-shuffle order until the class cap is reached."""
    import random
    ids = sorted(target_ids)
    random.Random(SHUFFLE_SEED).shuffle(ids)
    seen, rows, failed = set(), [], 0
    for n, tid in enumerate(ids, 1):
        if len(seen) >= CAP_PER_CLASS:
            break
        for attempt in range(4):
            try:
                rows += fetch_target(tid, seen, label, extra)
                break
            except Exception as e:
                if attempt == 3:
                    failed += 1
                    log(f"  [{label}] target {tid} gave up: {type(e).__name__}")
                else:
                    time.sleep(5 * (attempt + 1))
        if n % 25 == 0:
            log(f"  [{label}] targets {n}/{len(ids)}  compounds {len(seen)}  failed {failed}")
    log(f"  [{label}] DONE compounds={len(seen)} targets_used={n} failed={failed}")
    return rows


def write(df, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    p = os.path.abspath(os.path.join(OUT_DIR, name))
    df.to_csv(p, index=False)
    log(f"  written -> {p}  rows={len(df)}")
    return p


def main():
    map_path = os.path.abspath(os.path.join(OUT_DIR, "target_class_map.csv"))
    if os.path.exists(map_path):
        tm = pd.read_csv(map_path)
        log(f"target_class_map.csv reused ({len(tm)} rows)")
    else:
        anc = build_class_ancestry()
        tmap = build_target_class_map(anc)
        tm = pd.DataFrame([{"target_chembl_id": t, "target_class": lab}
                           for t, labs in tmap.items() for lab in labs])
        write(tm, "target_class_map.csv")

    by_class = {lab: set(g.target_chembl_id)
                for lab, g in tm.groupby("target_class")}
    log("\ntargets per class:")
    for lab in CLASSES:
        log(f"  {lab:18s} {len(by_class.get(lab, ())):6d}")

    frames = []
    for lab in CLASSES:
        tids = by_class.get(lab)
        if not tids:
            log(f"\n[{lab}] no targets, skipping")
            continue
        log(f"\n=== {lab}: {len(tids)} targets ===")
        d = pd.DataFrame(fetch_class(lab, tids))
        if len(d):
            write(d, f"corpus_{lab}.csv")
            frames.append(d)

    # whole-cell antibacterial bin (target_type=ORGANISM)
    log("\n=== antibacterial_wholecell ===")
    orgs = pd.DataFrame(list(paged("target", "targets", target_type="ORGANISM",
                                   only="target_chembl_id,pref_name,organism,tax_id")))
    bact = orgs[orgs.organism.fillna("").str.split().str[0].isin(BACTERIAL_GENERA)]
    log(f"  bacterial ORGANISM targets: {len(bact)}")
    extra = {r.target_chembl_id: {"organism": r.organism} for r in bact.itertuples()}
    d = pd.DataFrame(fetch_class("antibacterial_wholecell",
                                 set(bact.target_chembl_id), extra))
    if len(d):
        write(d, "corpus_antibacterial_wholecell.csv")
        frames.append(d)

    df = pd.concat(frames, ignore_index=True)
    write(df, "chembl_targetclass_corpus.csv")
    log("\n" + df.target_class.value_counts().to_string())
    log(f"\ntotal rows={len(df)} unique compounds={df.molecule_chembl_id.nunique()}")


if __name__ == "__main__":
    main()
