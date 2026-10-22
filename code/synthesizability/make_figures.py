# -*- coding: utf-8 -*-
"""
Regenerate every calibration table and figure from data/corpus_scored.csv.

Inputs : data/corpus_scored.csv, data/jak_scored.csv, data/library_structural_scored.csv
Outputs: data/class_summary.csv, data/calibration_thresholds.csv,
         data/global_vs_calibrated.csv, data/discrimination_auc.csv,
         data/jak_recalibrated.csv, and the four figures.

Run after build_corpus.py and score_corpus.py.
"""
import os
import sys
import tempfile

# Must precede the matplotlib import: on some sandboxed filesystems matplotlib's
# default config directory cannot be created, which aborts the import.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "mplconfig"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))

METRICS = ["sascore", "nplikeness", "scscore", "rascore", "brsascore"]
GATING = [m for m in METRICS if m != "nplikeness"]
THRESH = {"sascore": 6.0, "rascore": 0.5}
DIRN = {"sascore": "lower_easier", "nplikeness": "higher_easier",
        "scscore": "lower_easier", "rascore": "higher_easier",
        "brsascore": "lower_easier"}
LAB = {"sascore": "SAscore", "nplikeness": "NP-likeness", "scscore": "SCScore",
       "rascore": "RAscore", "brsascore": "BR-SAScore"}
ARROW = {"sascore": "lower = easier", "nplikeness": "higher = more NP-like",
         "scscore": "lower = easier", "rascore": "higher = easier",
         "brsascore": "lower = easier"}
FOCAL = "antibacterial_wholecell"
GREY = "#8c8c8c"
BLUE = "#1f77b4"
RED = "#d62728"
GREEN = "#2ca02c"
PURPLE = "#9467bd"
RATES = [0.05, 0.10]


def style():
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
        "axes.titlelocation": "left", "axes.titleweight": "regular",
        "figure.dpi": 110, "savefig.dpi": 300, "axes.grid": False,
    })


def frame(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(.8)
        ax.spines[s].set_color("#444444")
    ax.tick_params(length=2.5, width=.8, colors="#444444")


def letter(ax, ch):
    ax.text(-0.02, 1.06, ch, transform=ax.transAxes, fontweight="bold",
            fontsize=9, ha="right", va="bottom")


def pretty(c):
    return (c.replace("_", " ").replace("gpcr familyA", "GPCR (family A)")
             .replace("antibacterial wholecell", "antibacterial (whole-cell)"))


def frej(v, thr, m):
    return float((v > thr).mean() if DIRN[m] == "lower_easier" else (v < thr).mean())


def build_tables(sc):
    summ = []
    for cls, g in sc.groupby("target_class"):
        r = {"target_class": cls, "n": len(g)}
        for m in METRICS:
            r[f"{m}_median"] = round(g[m].median(), 3)
            r[f"{m}_iqr"] = round(g[m].quantile(.75) - g[m].quantile(.25), 3)
        r["pct_sascore_gt6"] = round(100 * (g.sascore > 6.0).mean(), 2)
        r["pct_rascore_lt05"] = round(100 * (g.rascore < 0.5).mean(), 2)
        summ.append(r)
    S = pd.DataFrame(summ).sort_values("pct_rascore_lt05", ascending=False)
    S.to_csv(os.path.join(DATA, "class_summary.csv"), index=False)

    rows = []
    for cls, g in sc.groupby("target_class"):
        for m in METRICS:
            v = g[m].dropna().values
            rec = dict(target_class=cls, metric=m, direction=DIRN[m], n=len(v),
                       median=round(float(np.median(v)), 4))
            for fr in RATES:
                q = (1 - fr) if DIRN[m] == "lower_easier" else fr
                rec[f"thr_fr{int(fr*100)}"] = round(float(np.quantile(v, q)), 4)
            gt = THRESH.get(m)
            rec["published_threshold"] = gt
            rec["published_false_reject_pct"] = (None if gt is None
                                                 else round(100 * frej(v, gt, m), 2))
            rows.append(rec)
    CAL = pd.DataFrame(rows)
    CAL.to_csv(os.path.join(DATA, "calibration_thresholds.csv"), index=False)

    glob = {}
    for m in METRICS:
        v = sc[m].dropna().values
        q = .95 if DIRN[m] == "lower_easier" else .05
        glob[m] = float(np.quantile(v, q))
    rows = []
    for cls, g in sc.groupby("target_class"):
        for m in METRICS:
            v = g[m].dropna().values
            own = float(CAL[(CAL.target_class == cls) & (CAL.metric == m)].thr_fr5.iloc[0])
            r = dict(target_class=cls, metric=m, n=len(v),
                     global_cutoff=round(glob[m], 4), class_cutoff=round(own, 4),
                     global_false_reject_pct=round(100 * frej(v, glob[m], m), 2),
                     class_false_reject_pct=round(100 * frej(v, own, m), 2))
            r["n_lost_global"] = int(round(r["global_false_reject_pct"] / 100 * len(v)))
            r["n_lost_class"] = int(round(r["class_false_reject_pct"] / 100 * len(v)))
            r["excess_lost"] = r["n_lost_global"] - r["n_lost_class"]
            pub = THRESH.get(m)
            r["published_cutoff"] = pub
            r["published_false_reject_pct"] = (None if pub is None
                                               else round(100 * frej(v, pub, m), 2))
            r["n_lost_published"] = (None if pub is None
                                     else int(round(frej(v, pub, m) * len(v))))
            rows.append(r)
    G = pd.DataFrame(rows)
    G.to_csv(os.path.join(DATA, "global_vs_calibrated.csv"), index=False)
    return S, CAL, G


def fig_class_distributions(sc, S, order):
    fig, axes = plt.subplots(1, 5, figsize=(9.0, 3.5), sharey=True)
    cols = {c: (BLUE if c == FOCAL else GREY) for c in order}
    for k, m in enumerate(METRICS):
        ax = axes[k]
        data = [sc.loc[sc.target_class == c, m].dropna().values for c in order]
        bp = ax.boxplot(data, orientation="horizontal", widths=.62, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black", lw=1.0))
        for patch, c in zip(bp["boxes"], order):
            patch.set_facecolor(cols[c]); patch.set_alpha(.55 if c == FOCAL else .28)
            patch.set_edgecolor(cols[c]); patch.set_linewidth(.8)
        for wk in bp["whiskers"] + bp["caps"]:
            wk.set_color(GREY); wk.set_linewidth(.7)
        if m in THRESH:
            ax.axvline(THRESH[m], color=RED, lw=1.0, ls="--", zorder=0)
            ax.text(THRESH[m], len(order) + 0.72, f"published\ncutoff {THRESH[m]:g}",
                    color=RED, fontsize=6, ha="center", va="bottom", linespacing=.95)
        ax.set_xlabel(LAB[m]); ax.margins(x=.04)
        ax.text(.5, -.30, ARROW[m], transform=ax.transAxes, ha="center", va="top",
                fontsize=6, color=GREY)
        frame(ax)
    axes[0].set_yticks(np.arange(1, len(order) + 1))
    axes[0].set_yticklabels([pretty(c) for c in order])
    fig.suptitle("Antibacterial chemistry is rated hardest on three of the four "
                 "synthesizability metrics", x=.02, ha="left", y=1.02)
    fig.tight_layout(rect=[0, .03, 1, .97]); fig.subplots_adjust(wspace=.10)
    fig.savefig(os.path.join(DATA, "class_distributions.png"), bbox_inches="tight")
    plt.close(fig)


def fig_roc_and_thresholds(sc, jk, CAL, order):
    lo, hi = int(jk.n_heavy.quantile(.05)), int(jk.n_heavy.quantile(.95))
    cm = sc[(sc.n_heavy >= lo) & (sc.n_heavy <= hi)]
    gm = jk[(jk.n_heavy >= lo) & (jk.n_heavy <= hi)]
    rng = np.random.default_rng(0)
    rows = []
    for m in METRICS:
        a = cm[m].dropna().values; b = gm[m].dropna().values
        sign = 1.0 if DIRN[m] == "higher_easier" else -1.0
        y = np.r_[np.ones(len(a)), np.zeros(len(b))]
        auc = roc_auc_score(y, sign * np.r_[a, b])
        boots = [roc_auc_score(y, sign * np.r_[a[rng.integers(0, len(a), len(a))],
                                               b[rng.integers(0, len(b), len(b))]])
                 for _ in range(300)]
        rows.append(dict(metric=m, auc=round(auc, 4),
                         ci_lo=round(float(np.percentile(boots, 2.5)), 4),
                         ci_hi=round(float(np.percentile(boots, 97.5)), 4),
                         n_known=len(a), n_generated=len(b),
                         median_known=round(float(np.median(a)), 3),
                         median_generated=round(float(np.median(b)), 3)))
    AUC = pd.DataFrame(rows).sort_values("auc", ascending=False)
    AUC.to_csv(os.path.join(DATA, "discrimination_auc.csv"), index=False)

    PAL = {"sascore": BLUE, "rascore": "#ff7f0e", "scscore": GREEN,
           "brsascore": PURPLE, "nplikeness": GREY}
    piv = CAL.pivot(index="target_class", columns="metric", values="thr_fr5")
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), width_ratios=[1, 1.12])
    ax = axes[0]
    for m in AUC.metric:
        a = cm[m].dropna().values; b = gm[m].dropna().values
        sign = 1.0 if DIRN[m] == "higher_easier" else -1.0
        y = np.r_[np.ones(len(a)), np.zeros(len(b))]
        fpr, tpr, _ = roc_curve(y, sign * np.r_[a, b])
        r = AUC[AUC.metric == m].iloc[0]
        ax.plot(fpr, tpr, color=PAL[m], lw=1.3,
                label=f"{LAB[m]}  {r.auc:.2f} [{r.ci_lo:.2f}-{r.ci_hi:.2f}]")
    ax.plot([0, 1], [0, 1], color=GREY, lw=.8, ls=":")
    ax.set_xlabel("fraction of generated molecules accepted")
    ax.set_ylabel("fraction of known chemistry accepted")
    ax.set_title("No metric cleanly separates known from\ngenerated chemistry at matched size")
    ax.legend(frameon=False, loc="lower right", fontsize=6, title="AUC [95% CI]",
              title_fontsize=6, handlelength=1.4)
    ax.set_xlim(-.02, 1.02); ax.set_ylim(-.02, 1.02); frame(ax)

    ax = axes[1]
    for i, c in enumerate(order[::-1]):
        v = piv.loc[c, "sascore"]
        col = BLUE if c == FOCAL else GREY
        ax.plot([6.0, v], [i, i], color=col, lw=1.0, alpha=.5, zorder=1)
        ax.plot(v, i, "o", color=col, ms=5, zorder=2)
    ax.axvline(6.0, color=RED, lw=1.0, ls="--", zorder=0)
    ax.text(6.0, len(order) - 0.4, "  published\n  global cutoff 6.0", color=RED,
            fontsize=6, ha="left", va="top", linespacing=1.0)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels([pretty(c) for c in order[::-1]])
    ax.set_xlabel("SAscore cutoff that discards 5% of known chemistry")
    ax.set_title("A single global cutoff is too strict for antibacterials\nand too loose for kinases")
    ax.set_xlim(4.0, 7.8); frame(ax)
    for a, l in zip(axes, "ab"):
        letter(a, l)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA, "roc_panel.png"), bbox_inches="tight")
    plt.close(fig)
    return cm, gm


def fig_cost(G, order):
    ordP = G[G.metric == "rascore"].sort_values("published_false_reject_pct").target_class.tolist()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
    for ax, m in zip(axes, ["sascore", "rascore"]):
        sub = G[G.metric == m].set_index("target_class").loc[ordP]
        vals = sub.published_false_reject_pct.values
        ax.barh(np.arange(len(ordP)), vals,
                color=[BLUE if c == FOCAL else GREY for c in ordP],
                edgecolor="none", height=.66)
        for i, (c, v) in enumerate(zip(ordP, vals)):
            ax.text(v + max(vals) * .015, i, f"{v:.1f}%", va="center", ha="left",
                    fontsize=6, color=(BLUE if c == FOCAL else GREY))
        ax.axvline(5.0, color=RED, lw=1.0, ls="--", zorder=3)
        ax.set_xlabel("known chemistry discarded (%)")
        ax.set_title(f"{LAB[m]} at its published cutoff of {THRESH[m]:g}")
        ax.set_xlim(0, max(vals) * 1.22); frame(ax)
    axes[0].set_yticks(np.arange(len(ordP)))
    axes[0].set_yticklabels([pretty(c) for c in ordP])
    axes[1].text(5.0, len(ordP) - 0.3, "  5% target", color=RED, fontsize=6,
                 ha="left", va="top")
    for a, l in zip(axes, "ab"):
        letter(a, l)
    fig.suptitle("At a published cutoff, how much real chemistry is thrown away "
                 "depends on the target class", x=.02, ha="left", y=1.03)
    fig.tight_layout(rect=[0, 0, 1, .97])
    fig.savefig(os.path.join(DATA, "cost_of_global_threshold.png"), bbox_inches="tight")
    plt.close(fig)


def fig_jak(sc, jk, CAL):
    struct_path = os.path.join(DATA, "library_structural_scored.csv")
    struct = pd.read_csv(struct_path)[["compound_id", "multi_fused_units"]]
    jk2 = jk.merge(struct, on="compound_id", how="left")
    jk2["artifact_sig"] = jk2.multi_fused_units >= 2
    KIN = {m: float(CAL[(CAL.target_class == "kinase") & (CAL.metric == m)].thr_fr5.iloc[0])
           for m in METRICS}
    kin = sc[sc.target_class == "kinase"]
    cand = jk2[jk2.compound_id == "cand_004995"].iloc[0]

    MS = ["sascore", "scscore", "rascore", "brsascore"]
    fig, axes = plt.subplots(1, 4, figsize=(9.0, 2.9), sharey=True)
    for k, m in enumerate(MS):
        ax = axes[k]
        bp = ax.boxplot([jk2[m].dropna().values, kin[m].dropna().values],
                        orientation="horizontal", widths=.5, showfliers=False,
                        patch_artist=True, medianprops=dict(color="black", lw=1.0))
        for patch, col in zip(bp["boxes"], [PURPLE, GREY]):
            patch.set_facecolor(col); patch.set_alpha(.45)
            patch.set_edgecolor(col); patch.set_linewidth(.8)
        for wk in bp["whiskers"] + bp["caps"]:
            wk.set_color(GREY); wk.set_linewidth(.7)
        ax.axvline(KIN[m], color=GREEN, lw=1.1, zorder=0)
        if m in THRESH:
            ax.axvline(THRESH[m], color=RED, lw=1.1, ls="--", zorder=0)
        ax.plot(float(cand[m]), 1, marker="D", color=RED, ms=5.5, zorder=5,
                markeredgecolor="white", markeredgewidth=.6)
        ax.set_xlabel(LAB[m]); ax.margins(x=.05); frame(ax)
    axes[0].set_yticks([1, 2])
    axes[0].set_yticklabels(["generated\nlibrary", "known kinase\nchemistry"])
    h = [plt.Line2D([], [], color=GREEN, lw=1.1, label="kinase-calibrated cutoff (5% false reject)"),
         plt.Line2D([], [], color=RED, lw=1.1, ls="--", label="published global cutoff"),
         plt.Line2D([], [], color=RED, marker="D", ls="none", ms=5,
                    label="cand_004995 (chemist-rejected)")]
    fig.legend(handles=h, frameon=False, fontsize=6.5, loc="lower center",
               bbox_to_anchor=(.5, -.13), ncol=3)
    fig.suptitle("No class-appropriate cutoff separates the chemist-rejected molecule "
                 "from real kinase chemistry", x=.02, ha="left", y=1.04)
    fig.tight_layout()
    fig.savefig(os.path.join(DATA, "jak_recalibrated.png"), bbox_inches="tight")
    plt.close(fig)

    out = jk2[["compound_id", "std_smiles", "multi_fused_units", "artifact_sig"] + METRICS].copy()
    for m in MS:
        out[f"reject_kinase_{m}"] = (out[m] > KIN[m]) if DIRN[m] == "lower_easier" else (out[m] < KIN[m])
        if m in THRESH:
            out[f"reject_published_{m}"] = ((out[m] > THRESH[m]) if DIRN[m] == "lower_easier"
                                            else (out[m] < THRESH[m]))
    out.to_csv(os.path.join(DATA, "jak_recalibrated.csv"), index=False)


def main():
    style()
    sc = pd.read_csv(os.path.join(DATA, "corpus_scored.csv"))
    jk = pd.read_csv(os.path.join(DATA, "jak_scored.csv"))
    if "n_heavy" not in jk.columns:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        jk["n_heavy"] = [(lambda x: x.GetNumHeavyAtoms() if x else None)(
            Chem.MolFromSmiles(str(s))) for s in jk.std_smiles.fillna("")]
    S, CAL, G = build_tables(sc)
    order = S.sort_values("pct_rascore_lt05").target_class.tolist()
    fig_class_distributions(sc, S, order)
    fig_roc_and_thresholds(sc, jk, CAL, order)
    fig_cost(G, order)
    fig_jak(sc, jk, CAL)
    print("regenerated tables and 4 figures in", DATA, flush=True)


if __name__ == "__main__":
    main()
