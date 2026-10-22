# -*- coding: utf-8 -*-
"""
artifact_analysis.py
Reproduce the SAscore-fails-on-artifacts finding across the full JAK library,
without needing the AiZynthFinder model (which is blocked in-sandbox).

Idea: BRICS fragment recombination produces molecules with chained/fused
polycyclic systems (often repeated fragments). SAscore, being an atom/fragment
complexity proxy, does NOT penalize these route-implausible artifacts. We
isolate the structural artifact population and show how many SAscore rates EASY.
This is the population that requires Tier 3 (retrosynthesis) adjudication.
"""
import os, sys
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

df = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  "..", "..", "data", "jak_library.csv"))

def ring_systems(mol):
    """Return list of fused-ring-system sizes (rings per fused cluster)."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    if not rings:
        return []
    # union-find fusion by shared atoms
    parent = list(range(len(rings)))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        parent[find(a)]=find(b)
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            if rings[i]&rings[j]:
                union(i,j)
    from collections import Counter
    comp = Counter(find(i) for i in range(len(rings)))
    return sorted(comp.values(), reverse=True)

def largest_fused_atoms(mol):
    ri = mol.GetRingInfo()
    rings=[set(r) for r in ri.AtomRings()]
    if not rings: return 0
    parent=list(range(len(rings)))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b): parent[find(a)]=find(b)
    for i in range(len(rings)):
        for j in range(i+1,len(rings)):
            if rings[i]&rings[j]: union(i,j)
    from collections import defaultdict
    clusters=defaultdict(set)
    for i,r in enumerate(rings): clusters[find(i)] |= r
    return max(len(v) for v in clusters.values())

rows=[]
for _,r in df.iterrows():
    m=Chem.MolFromSmiles(str(r.smile))
    if m is None:
        rows.append(dict(compound_id=r.compound_id, valid=False)); continue
    rs=ring_systems(m)
    n_rings=rdMolDescriptors.CalcNumRings(m)
    # "chained fused-bicyclic" signature: >=2 distinct fused ring systems that
    # each contain >=2 rings (i.e. multiple bicyclic+ units linked together)
    multi_fused_units = sum(1 for s in rs if s>=2)
    rows.append(dict(
        compound_id=r.compound_id,
        smile=r.smile,
        sa_score=r.sa_score,
        n_rings=n_rings,
        n_ring_systems=len(rs),
        multi_fused_units=multi_fused_units,
        largest_fused_ringsys_atoms=largest_fused_atoms(m),
        valid=True,
    ))
a=pd.DataFrame(rows)
a=a[a.valid].copy()

# SAscore verdict at the tool's own threshold
a["sa_verdict"]=pd.cut(a.sa_score,[-1,3.5,6.0,99],labels=["EASY","MODERATE","HARD"])

# Structural artifact signature (route-implausible BRICS recombination):
#  >=2 separate multi-ring (fused bicyclic+) systems chained in one molecule.
a["artifact_sig"] = a.multi_fused_units>=2

print("=== Library ===")
print("valid molecules:", len(a))
print("\nSAscore verdict counts:")
print(a.sa_verdict.value_counts().to_string())

print("\n=== Structural artifact signature (>=2 chained fused-bicyclic+ units) ===")
print("flagged molecules:", int(a.artifact_sig.sum()),
      f"({100*a.artifact_sig.mean():.1f}%)")

print("\n=== The failure population: artifact-signature AND SAscore=EASY ===")
fail = a[(a.artifact_sig) & (a.sa_verdict=="EASY")]
print("count:", len(fail), f"({100*len(fail)/len(a):.1f}% of library)")
print("their sa_score: min %.2f  median %.2f  max %.2f"%(fail.sa_score.min(),fail.sa_score.median(),fail.sa_score.max()))
print("  (SAscore of imatinib ~ 3.0-3.7 for reference; these look 'easier' than a real drug)")

print("\ncand_004995 in failure population?", "cand_004995" in set(fail.compound_id))
c=a[a.compound_id=="cand_004995"].iloc[0]
print("  cand_004995: sa=%.2f  n_rings=%d  n_ring_systems=%d  multi_fused_units=%d  sa_verdict=%s"%(
    c.sa_score,c.n_rings,c.n_ring_systems,c.multi_fused_units,c.sa_verdict))

print("\n=== cross-tab: artifact_sig x sa_verdict ===")
print(pd.crosstab(a.artifact_sig, a.sa_verdict).to_string())

# save the flagged failure population, most-artifact-like first
fail_sorted=fail.sort_values(["multi_fused_units","n_rings","sa_score"],ascending=[False,False,True])
fail_sorted.to_csv("sascore_failure_population.csv",index=False)
a.to_csv("library_structural_scored.csv",index=False)
print("\nwritten -> sascore_failure_population.csv (%d rows)"%len(fail_sorted))
print("written -> library_structural_scored.csv (%d rows)"%len(a))
print("\nTop 10 worst artifacts SAscore calls EASY:")
print(fail_sorted[["compound_id","sa_score","n_rings","n_ring_systems","multi_fused_units"]].head(10).to_string(index=False))
