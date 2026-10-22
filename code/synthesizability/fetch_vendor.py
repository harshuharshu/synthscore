# -*- coding: utf-8 -*-
"""
Fetch the third-party model weights the scorer panel needs.

The weights are published by their respective authors and are NOT redistributed
in this repository. This script downloads them from the authors' own GitHub
repositories into `vendor/`, where `panel.py` expects them.

Only three files are needed. The panel's wrapper code for SCScore, RAscore and
BR-SAScore was written here (see vendor/panel.py and vendor/scscore/
scscore_adapter.py), so none of the authors' source files are required.

Usage:
    python fetch_vendor.py
    GITHUB_TOKEN=... python fetch_vendor.py      # higher API rate limit

Set GITHUB_TOKEN if you hit a 403 rate-limit; unauthenticated requests are
capped at 60 per hour, and this script makes six.
"""
import base64
import json
import os
import sys
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
VENDOR = os.path.join(HERE, "vendor")

# (owner, repo, path_in_repo, local_destination_relative_to_vendor)
TARGETS = [
    ("connorcoley", "scscore",
     "models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz",
     os.path.join("scscore", "model.ckpt-10654.as_numpy.json.gz")),
    ("reymond-group", "RAscore",
     "RAscore/models/XGB_chembl_ecfp_counts/model.pkl",
     os.path.join("rascore", "XGB_chembl_ecfp_counts_model.pkl")),
    ("awadell1", "BR-SAScore",
     "BRSAScore/pickle/BRScores_uspto_emolecules.pkl.gz",
     os.path.join("brsascore", "pickle", "BRScores_uspto_emolecules.pkl.gz")),
]


def _headers():
    h = {"Accept": "application/vnd.github+json", "User-Agent": "synthscore-fetch"}
    tok = os.environ.get("GITHUB_TOKEN")
    if tok:
        h["Authorization"] = f"Bearer {tok}"
    return h


def _api(url):
    req = urllib.request.Request(url, headers=_headers())
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode())


def fetch(owner, repo, path, dest):
    out = os.path.join(VENDOR, dest)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print(f"  already present, skipping: {dest}")
        return True
    tree = _api(f"https://api.github.com/repos/{owner}/{repo}/git/trees/HEAD?recursive=1")
    entry = next((e for e in tree.get("tree", []) if e["path"] == path), None)
    if entry is None:
        print(f"  NOT FOUND upstream: {owner}/{repo}:{path}", file=sys.stderr)
        return False
    blob = _api(f"https://api.github.com/repos/{owner}/{repo}/git/blobs/{entry['sha']}")
    data = (base64.b64decode(blob["content"]) if blob.get("encoding") == "base64"
            else blob["content"].encode())
    if data[:7] == b"version" and b"git-lfs" in data[:200]:
        print(f"  upstream file is a Git LFS pointer, not content: {path}", file=sys.stderr)
        return False
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "wb") as fh:
        fh.write(data)
    print(f"  {len(data)/1e6:6.2f} MB  {dest}")
    return True


def main():
    print(f"fetching model weights into {VENDOR}")
    ok = 0
    for owner, repo, path, dest in TARGETS:
        try:
            ok += bool(fetch(owner, repo, path, dest))
        except urllib.error.HTTPError as e:
            print(f"  HTTP {e.code} fetching {dest}"
                  f"{' (rate limit; set GITHUB_TOKEN)' if e.code == 403 else ''}",
                  file=sys.stderr)
        except Exception as e:
            print(f"  {type(e).__name__} fetching {dest}: {str(e)[:120]}", file=sys.stderr)
    print(f"\n{ok}/{len(TARGETS)} weight files present.")
    if ok < len(TARGETS):
        print("Metrics whose weights are missing are skipped at runtime with a warning; "
              "SAscore and NP-likeness ship with RDKit and always work.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
