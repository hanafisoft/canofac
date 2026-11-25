"""Reproduce the application and a reduced simulation study of article BB (LOBO permutation test).

Usage:
  python3 articleBB.py --potato potato.csv --out OUTDIR [--nrep 100] [--B 999] [--lang en|fr]

potato.csv: exchange-format file, e.g. exported from R with
  cf_write_blocks(cf_potato(preprocess = FALSE), "potato.csv")
or the multiblock potato table written with write.csv (columns <block>.x.<variable>).
Outputs (OUTDIR): tables (CSV, printed) and figures (PDF).
The simulated data come from cf_simulate() (common generator); the published study used another generator, so the
rejection rates agree up to sampling fluctuations. --nrep 0 skips the simulations.
"""
import argparse, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import canofac as cf

METHODS = {"MAXBET": dict(q=1, extraction="sequential"), "CPCA": dict(q=2), "INDORT": dict(q=4, nstart=10),
           "MFA": dict(q=2, weighting="mfa"), "STATIS": dict(q=2, weighting="statis")}


def potato(file, out, B):
    bl = cf.cf_potato(file)
    rows = []
    for m, kw in METHODS.items():
        t = cf.cf_test(bl, R=2, B=B, seed=1, compare=True, **kw)
        r = {"method": m}
        for k, b in enumerate(t.block_names):
            r[f"p_sub {b}"] = t.p_subspace[k]
        r["linked (subspace)"] = t.link_subspace
        for h in range(2):
            r[f"min_k p_k{h + 1}"] = t.p[:, h].min()
            r[f"comp {h + 1} linked"] = bool(t.retained[h])
        for k, b in enumerate(t.block_names):
            r[f"RV {b}"] = t.p_compare.iloc[k, 0]
            r[f"PROTEST {b}"] = t.p_compare.iloc[k, 1]
        rows.append(r)
        if m == "CPCA":
            print(t)
            t.plot("null")
            plt.savefig(os.path.join(out, "figBB_potato_null.pdf"))
            t.plot()
            plt.savefig(os.path.join(out, "figBB_potato_pvalues.pdf"))
            plt.close("all")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "tableBB_potato.csv"), index=False)
    print("\nPotato, LOBO test (restricted RV), B =", B, "\n", df.round(3).to_string(index=False))
    # re-estimation statistic (comparison)
    tr = cf.cf_test(bl, R=2, q=2, B=min(B, 199), seed=1, statistic="refit")
    print("\nRe-estimation statistic (CPCA), subspace p-values:", np.round(tr.p_subspace, 3), "->", tr.link_subspace)
    # robustness to R
    prof = []
    for m in ("MAXBET", "CPCA", "INDORT"):
        for R in range(1, 6):
            t = cf.cf_test(bl, R=R, B=B, seed=1, **METHODS[m])
            prof.append({"method": m, "R": R, **{f"p_sub {b}": t.p_subspace[k] for k, b in enumerate(t.block_names)},
                         "linked components": "".join(str(h + 1) for h in range(R) if t.retained[h]) or "-"})
    dp = pd.DataFrame(prof)
    dp.to_csv(os.path.join(out, "tableBB_potato_R.csv"), index=False)
    print("\nPotato, robustness to R\n", dp.round(4).to_string(index=False))
    # p-values without permutation
    g = cf.cf_test(bl, R=2, q=2, seed=1, pvalue="gamma")
    h = cf.cf_test(bl, R=2, q=2, B=19999, seed=1, pvalue="hybrid")
    print("\nPotato, CPCA: min_k p_k2 by permutation / gamma / hybrid:",
          round(df.loc[df.method == "CPCA", "min_k p_k2"].item(), 4), round(g.p[:, 1].min(), 4), round(h.p[:, 1].min(), 4))


def simulations(out, nrep, B=99):
    configs = [("S0", 0.0), ("SA", 1.1), ("SA", 1.6), ("SB", 1.8), ("SC", 1.8), ("SD", 1.8), ("SE", 1.8)]
    carriers = {"S0": [], "SA": [0, 1, 2, 3], "SB": [0, 1, 2, 3], "SC": [0, 1], "SD": [0, 1, 2, 3], "SE": [0, 1, 2, 3]}
    rows = []
    for scen, th in configs:
        for n in (20, 100):
            rej = {s: [] for s in ("rv", "eq", "RV", "PROTEST")}
            exact = {s: [] for s in ("rv", "eq")}
            for rep in range(nrep):
                bl = cf.cf_simulate(scen, n=n, p=10, K=4, theta=th, seed=100000 * (n == 100) + 1000 * configs.index((scen, th)) + rep + 1)
                t = cf.cf_test(bl, R=2, q=2, B=B, seed=rep + 1, compare=True)
                rej["rv"].append(t.p_subspace_inertia <= 0.05)
                rej["eq"].append(t.p_subspace_equal <= 0.05)
                rej["RV"].append(t.p_compare["RV"].to_numpy() <= 0.05)
                rej["PROTEST"].append(t.p_compare["PROTEST"].to_numpy() <= 0.05)
                tgt = np.array([k in carriers[scen] for k in range(4)])
                for s, pv in (("rv", t.p_subspace_inertia), ("eq", t.p_subspace_equal)):
                    L = np.ones(4, bool) if pv.max() <= 0.05 else cf.holm(pv) <= 0.05
                    exact[s].append(bool((L == tgt).all()))
            r = {"scenario": scen, "theta": th, "n": n}
            for s, v in rej.items():
                r[f"rej {s} (%)"] = float(round(100 * np.mean(v), 1))
            for s, v in exact.items():
                r[f"exact {s} (%)"] = float(round(100 * np.mean(v), 1))
            rows.append(r)
            print(r, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "tableBB_simulations.csv"), index=False)
    print("\nSimulations (K = 4, p = 10, CPCA, R = 2, B =", B, ", nrep =", nrep, ")\n", df.to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--potato"); p.add_argument("--out", default="articleBB_results")
    p.add_argument("--nrep", type=int, default=100); p.add_argument("--B", type=int, default=999)
    p.add_argument("--lang", default="en")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cf.cf_lang(a.lang)
    if a.potato and os.path.exists(a.potato):
        potato(a.potato, a.out, a.B)
    if a.nrep > 0:
        simulations(a.out, a.nrep)
