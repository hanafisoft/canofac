"""Reproduce the tables and figures of article A (pseudo-singular values and the power-q family).

Usage:
  python3 articleA.py --potato potato.csv --friday friday.csv --tcga tcga.csv --out OUTDIR [--lang en|fr]

Data (not distributed with canofac): exchange-format CSV files exported from R, e.g.
  cf_write_blocks(cf_potato(blocks = "A", preprocess = FALSE), "potato.csv")
  cf_write_blocks(cf_friday(preprocess = FALSE), "friday.csv")
  tc <- cf_tcga(preprocess = FALSE); cf_write_blocks(tc$blocks, "tcga.csv", extra = data.frame(subtype = tc$subtype))
A missing dataset is skipped.
Outputs: tables (CSV and printed) and figures (PDF) in OUTDIR.
Note: the random starts come from the common generator of canofac (seed 2026), as in the article: the R and Python
versions give the same tables.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import canofac as cf

QGRID = [1, 1.25, 1.5, 1.75, 1.9, 2, 2.1, 2.5, 3, 4, 6, 8, 10]


def load(args):
    data = {}
    if args.potato and os.path.exists(args.potato):
        data["Potato"] = (cf.cf_potato(args.potato, blocks="A"), None)
    if args.tcga and os.path.exists(args.tcga):
        data["breast.TCGA"] = cf.cf_tcga(args.tcga)
    if args.friday and os.path.exists(args.friday):
        data["Friday"] = (cf.cf_friday(args.friday), None)
    return data


def multimodality(data, out):
    rows = []
    for name, (bl, _) in data.items():
        for ext in ("simultaneous", "sequential"):
            row = {"dataset": name, "extraction": ext}
            for q in (1, 2, 4, 6, 10):
                ms = cf.cf_multistart(bl, R=2, q=q, nstart=100, extraction=ext, seed=2026)
                star = "" if ms.deterministic_best else "*"
                row[f"q={q}"] = f"{ms.share_best:.0f}" + (f" ({ms.n_optima})" if ms.n_optima > 1 else "") + star
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "tableA_multimodality.csv"), index=False)
    print("\nMultimodality (% of 100 random starts reaching the best criterion; * CPCA start not best)\n", df.to_string(index=False))


def trajectories(data, out):
    for name, (bl, _) in data.items():
        tr = cf.cf_trajectories(bl, R=2, qgrid=QGRID, nstart=20)
        ax = tr.plot()
        ax.set_title(name)
        plt.savefig(os.path.join(out, f"figA_traj_{name.split()[0]}.pdf"))
        plt.close("all")


def by_method(data, out):
    for name, (bl, _) in data.items():
        fits = [cf.canofac(bl, 2, 1), cf.canofac(bl, 2, 1, axes="orthonormal"), cf.canofac(bl, 2, 2),
                cf.canofac(bl, 2, 4, nstart=10)]
        names = ["q=1, unconstrained", "simultaneous MAXBET", "CPCA", "INDORT"]
        print(f"\n{name}")
        df = cf.cf_compare(fits, names=names, match=False)
        df.to_csv(os.path.join(out, f"tableA_methods_{name.split()[0]}.csv"))


def seq_sim(data, out):
    if "Friday" not in data:
        return
    bl = data["Friday"][0]
    rows = []
    for lab, q, ax in [("q=1, unconstrained", 1, "free"), ("q=2, unconstrained", 2, "free"),
                       ("q=2, orthonormal", 2, "orthonormal"), ("q=4, unconstrained", 4, "free")]:
        fs = cf.canofac(bl, 2, q, ax, "sequential")
        fm = cf.canofac(bl, 2, q, ax, "simultaneous")
        a = cf.cf_angles(fs, fm)
        rows.append({"criterion": lab, "angle": round(float(a["components"].max()), 1),
                     "crit_seq": round(fs.criterion, 3), "crit_sim": round(fm.criterion, 3)})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "tableA_seq_sim_Friday.csv"), index=False)
    print("\nFriday: sequential versus simultaneous\n", df.to_string(index=False))


def tcga_individuals(data, out):
    if "breast.TCGA" not in data:
        return
    bl, sub = data["breast.TCGA"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.5))
    for a, (q, lab) in zip(axs, [(2, "CPCA"), (4, "INDORT")]):
        cf.cf_plot_scores(cf.canofac(bl, 2, q, nstart=20), group=sub, ax=a, title=lab)
    plt.savefig(os.path.join(out, "figA_indiv_tcga.pdf"))
    plt.close("all")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--potato"); p.add_argument("--friday"); p.add_argument("--tcga")
    p.add_argument("--out", default="articleA_results"); p.add_argument("--lang", default="en")
    args = p.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cf.cf_lang(args.lang)
    data = load(args)
    print("datasets:", ", ".join(data))
    multimodality(data, args.out)
    trajectories(data, args.out)
    by_method(data, args.out)
    seq_sim(data, args.out)
    tcga_individuals(data, args.out)
