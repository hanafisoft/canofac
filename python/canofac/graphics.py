"""Plots of canofac (matplotlib), with the same content and labels as the R package."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from scipy import stats
from .lang import _tr, _link_label, cf_lang
from ._core import match_dims
from .explore import method_name

PALETTE = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
           "#17becf", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]


def _pal(n):
    return [PALETTE[i % len(PALETTE)] for i in range(n)]


def _ax(ax, **kw):
    if ax is None:
        _, ax = plt.subplots(layout="constrained", **kw)
    return ax


def _ellipse(P, level=0.95, npts=100):
    m = P.mean(axis=0)
    w, v = np.linalg.eigh(np.cov(P.T))
    r = np.sqrt(stats.chi2.ppf(level, 2))
    th = np.linspace(0, 2 * np.pi, npts)
    A = v @ np.diag(np.sqrt(np.maximum(w, 0)))
    return (m[:, None] + r * A @ np.vstack([np.cos(th), np.sin(th)])).T


def cf_plot_scores(x, dims=(1, 2), group=None, ellipse=None, level=0.95, labels=False, title=None, ax=None):
    """Individuals on the reference components (multiplied by sqrt(n)), with normal ellipses by group."""
    ax = _ax(ax, figsize=(6, 6))
    U = x.U[:, [dims[0] - 1, dims[1] - 1]] * np.sqrt(x.U.shape[0])
    if ellipse is None:
        ellipse = group is not None
    if group is None:
        ax.scatter(U[:, 0], U[:, 1], color="black", s=18)
    else:
        group = np.asarray(group)
        levels = list(dict.fromkeys(group.tolist()))
        for lev, col in zip(levels, _pal(len(levels))):
            P = U[group == lev]
            ax.scatter(P[:, 0], P[:, 1], color=col, s=18, label=str(lev))
            if ellipse and P.shape[0] >= 3:
                E = _ellipse(P, level)
                ax.plot(E[:, 0], E[:, 1], color=col, lw=1.5)
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    if labels is not False:
        lab = labels if not isinstance(labels, bool) else (x.row_names or list(range(1, U.shape[0] + 1)))
        for (a, b), t in zip(U, lab):
            ax.annotate(str(t), (a, b), fontsize=7, xytext=(0, 4), textcoords="offset points", ha="center")
    ax.axhline(0, ls=":", color="grey", lw=0.8)
    ax.axvline(0, ls=":", color="grey", lw=0.8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel(f"{_tr('ref_comp')} {dims[0]}")
    ax.set_ylabel(f"{_tr('ref_comp')} {dims[1]}")
    ax.set_title(title or _tr("scores_title"))
    return ax


def _heat(M, txt, title, rows, cols, bold=None, vmin=None, vmax=None, top=None, cmap="YlOrRd", ax=None):
    ax = _ax(ax, figsize=(1.4 + 1.1 * M.shape[1], 1.2 + 0.45 * M.shape[0]))
    ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            w = "bold" if (bold is not None and bold[i, j]) else "normal"
            ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=9, fontweight=w)
    ax.set_xticks(range(M.shape[1]))
    ax.set_xticklabels(cols)
    ax.set_yticks(range(M.shape[0]))
    ax.set_yticklabels(rows)
    if top is not None:
        sec = ax.secondary_xaxis("top")
        sec.set_xticks(range(M.shape[1]))
        sec.set_xticklabels(top, fontsize=8)
        sec.tick_params(length=0)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=10, pad=24 if top is not None else 6)
    return ax


def cf_plot_saliences(x, stat="saliences", digits=2, title=None, ax=None):
    """Heat map of the saliences, pseudo-singular values or delta_kh (K x R)."""
    M = {"saliences": x.saliences, "d": x.d, "delta": x.delta}[stat]
    txt = [[f"{v:.{digits}f}" for v in row] for row in M]
    t = title or _tr({"saliences": "sal_title", "d": "d_title", "delta": "delta_title"}[stat])
    return _heat(M, txt, t, x.block_names, [f"{_tr('dim')} {h + 1}" for h in range(M.shape[1])],
                 vmin=0, vmax=max(M.max(), 1e-12), cmap="Blues", ax=ax)


def plot_fit(x, type="scores", dims=(1, 2), group=None, labels=False, block=None, ax=None, **kw):
    """Plot of a :class:`Canofac` fit."""
    if type == "scores":
        return cf_plot_scores(x, dims=dims, group=group, labels=labels, ax=ax, **kw)
    if type == "loadings":
        blocks = block if block is not None else x.block_names
        blocks = [blocks] if isinstance(blocks, str) else list(blocks)
        fig, axs = plt.subplots(1, len(blocks), figsize=(4.5 * len(blocks), 3.5), squeeze=False)
        for a, b in zip(axs[0], blocks):
            V = x.V[b][:, [d - 1 for d in dims if d <= x.U.shape[1]]]
            vn = x.variable_names.get(b)
            try:
                xx = np.array([float(v) for v in vn])
                spectral = True
            except (TypeError, ValueError):
                xx, spectral = np.arange(1, V.shape[0] + 1), False
            for j, col in enumerate(_pal(V.shape[1])):
                if spectral or V.shape[0] > 30:
                    a.plot(xx, V[:, j], color=col, lw=1.2, label=f"{_tr('dim')} {dims[j]}")
                else:
                    a.vlines(xx + (j - 0.5) * 0.2, 0, V[:, j], color=col, lw=3, label=f"{_tr('dim')} {dims[j]}")
            if spectral and xx[0] > xx[-1]:
                a.invert_xaxis()
            if not spectral and vn is not None and V.shape[0] <= 30:
                a.set_xticks(xx)
                a.set_xticklabels(vn, rotation=90, fontsize=7)
            a.axhline(0, ls=":", color="grey", lw=0.8)
            a.set_xlabel(_tr("wavelength") if spectral else _tr("variable"))
            a.set_ylabel(_tr("loading"))
            a.set_title(b)
            a.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        return axs
    M = {"d": x.d, "saliences": x.saliences, "delta": x.delta}[type]
    ax = _ax(ax, figsize=(1.5 + 0.8 * M.shape[0], 3.5))
    K, R = M.shape
    w = 0.8 / R
    for h, col in enumerate(_pal(R)):
        ax.bar(np.arange(K) + (h - (R - 1) / 2) * w, M[:, h], width=w, color=col, label=f"Dim{h + 1}")
    ax.set_xticks(range(K))
    ax.set_xticklabels(x.block_names, rotation=90)
    ax.set_ylabel({"d": _tr("pvs"), "saliences": _tr("salience"), "delta": r"$\delta_{kh}$"}[type])
    ax.legend(frameon=False, fontsize=8)
    return ax


def plot_trajectories(x, dims=(1, 2), labels=True, mark=(2, 4), title=None, legend=True, ax=None):
    """Trajectory of each block in the plane (delta_k1, delta_k2) as q varies."""
    ax = _ax(ax, figsize=(6, 6))
    D, q = x.delta, x.q
    K, nq = D.shape[0], len(q)
    cols = _pal(K)
    s = x.settings
    mq = [q[0]] + [m for m in mark if m in q[1:-1]] + [q[-1]]
    mk = ["o"] + (["s", "^", "D", "*", "+", "x"] * 3)[:len(mq) - 2] + ["o"]
    filled = [False] + [True] * (len(mq) - 2) + [True]
    ax.plot([0, 1], [0, 1], ls=":", color="grey", lw=0.8)
    i1, i2 = dims[0] - 1, dims[1] - 1
    for k in range(K):
        ax.plot(D[k, i1, :], D[k, i2, :], color=cols[k], lw=1.5)
        for v, m, f in zip(mq, mk, filled):
            j = int(np.where(q == v)[0][0])
            ax.scatter(D[k, i1, j], D[k, i2, j], marker=m, s=30,
                       facecolors=cols[k] if f else "none", edgecolors=cols[k])
        if labels:
            nm = labels[k] if isinstance(labels, (list, tuple)) else x.block_names[k]
            ax.annotate(nm, (D[k, i1, -1], D[k, i2, -1]), fontsize=7, color=cols[k],
                        xytext=(4, 0), textcoords="offset points")
    if legend:
        hs = []
        for v, m, f in zip(mq, mk, filled):
            nm = method_name(v, s["axes"], s["extraction"])
            lab = f"q = {v:g} ({nm})" if nm else f"q = {v:g}"
            hs.append(ax.scatter([], [], marker=m, facecolors="black" if f else "none", edgecolors="black", label=lab))
        ax.legend(handles=hs, frameon=False, fontsize=8, loc="lower right")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel(rf"$\delta_{{{dims[0]}}}$")
    ax.set_ylabel(rf"$\delta_{{{dims[1]}}}$")
    ax.set_title(title or _tr("traj_title"))
    return ax


cf_plot_trajectories = plot_trajectories


def plot_test(x, digits=1, subspace=True, title=None, ax=None):
    """Heat map of the p-values (%) by block and component; linked blocks in bold, outcome on top."""
    P = x.p.copy()
    fixed = x.settings["statistic"] == "lobo"
    cols = [f"{_tr('comp' if fixed else 'dim')} {h + 1}" for h in range(P.shape[1])]
    bold = x.linked.copy()
    top = _link_label(x.link)
    if subspace:
        P = np.column_stack([P, x.p_subspace])
        cols.append(_tr("subspace"))
        bold = np.column_stack([bold, x.linked_subspace])
        top = top + [_link_label([x.link_subspace])[0].upper()]
    txt = [[f"{100 * v:.{digits}f}" for v in row] for row in P]
    zmax = -np.log10(1 / (x.settings["B"] + 1)) if x.settings["pvalue"] == "permutation" else 4
    return _heat(-np.log10(np.maximum(P, 1e-300)), txt, title or _tr("pvalues_title"), x.block_names, cols,
                 bold=bold, vmin=0, vmax=zmax, top=top, ax=ax)


cf_plot_test = plot_test


def plot_null(x, title=None, bins=30):
    """Permutation distribution of the subspace statistic for each block; vertical line: observed value."""
    if not x.null:
        raise ValueError("no permutation distribution stored (pvalue = 'gamma')")
    names = list(x.null)
    fig, axs = plt.subplots(1, len(names), figsize=(3.2 * len(names), 2.8), squeeze=False)
    for a, nm in zip(axs[0], names):
        d = x.null[nm]
        a.hist(d["permuted"], bins=bins, color="#9DB9E6", edgecolor="white")
        a.axvline(d["observed"], color="#B4491F", lw=2)
        k = x.block_names.index(nm)
        a.set_title(f"{nm}  p = {x.p_subspace[k]:.3g}", fontsize=9)
        a.set_xlabel(_tr("statistic"))
        a.set_ylabel(_tr("frequency"))
    fig.suptitle(title or _tr("null_title"), fontsize=10)
    fig.tight_layout()
    return axs


def cf_plot_compare(fits, dim=1, names=None, title=None, ax=None):
    """Bar plot of the pseudo-singular values of several fits on one dimension (matched to the first fit)."""
    if names is None:
        names = []
        for f in fits:
            m = method_name(f.settings["q"], f.settings["axes"], f.settings["extraction"])
            names.append(m if m else f"q = {f.settings['q']:g}")
    ref = fits[0].U
    M = np.column_stack([f.d[:, match_dims(f.U, ref)[dim - 1]] for f in fits])
    ax = _ax(ax, figsize=(1.5 + 0.9 * M.shape[0], 3.5))
    K, m = M.shape
    w = 0.8 / m
    for j, col in enumerate(_pal(m)):
        ax.bar(np.arange(K) + (j - (m - 1) / 2) * w, M[:, j], width=w, color=col, label=names[j])
    ax.set_xticks(range(K))
    ax.set_xticklabels(fits[0].block_names, rotation=90)
    ax.set_ylabel(_tr("pvs"))
    ax.set_ylim(0, max(M.max(), 1e-12) * 1.15)
    ax.legend(frameon=False, fontsize=8, ncol=m, loc="upper center")
    ax.set_title(title or _tr("compare_title") % dim)
    return ax


def _block_colors(K):
    import colorsys
    return [colorsys.hsv_to_rgb(k / K, 0.65, 0.85) for k in range(K)]


def _rose(ax, values, vmax, r0, r1, colors, label=True, names=None):
    K = len(values)
    width = 360 / K * 0.88
    for k in range(K):
        th = k * 360 / K
        rk = r0 + (r1 - r0) * max(values[k], 0) / vmax
        if rk > r0:
            ax.add_patch(Wedge((0, 0), rk, th - width / 2, th + width / 2, width=rk - r0, color=colors[k], lw=0))
        if label:
            ax.text(rk * 1.08 * np.cos(np.radians(th)), rk * 1.08 * np.sin(np.radians(th)), f"{values[k]:.2f}",
                    fontsize=6, ha="center", va="center")
        if names is not None:
            ax.text(r1 * 1.14 * np.cos(np.radians(th)), r1 * 1.14 * np.sin(np.radians(th)), names[k],
                    fontsize=6, fontweight="bold", color=colors[k], ha="center", va="center")


def cf_rose(x, dim=1, stat="d", block_names=None, labels=True, legend=True, title=None):
    """Rose plot of the pseudo-singular values: inner, squared deviations from the mean; middle, values;
    outer ring, share of mean^2 / (mean^2 + V). ``x`` is a fit (one panel per dimension) or a list of fits."""
    val = (lambda f: f.d) if stat == "d" else (lambda f: f.saliences)
    if hasattr(x, "U"):
        M = val(x)
        panels = [M[:, h] for h in range(M.shape[1])]
        titles = [f"{_tr('dimension')} {h + 1}" for h in range(M.shape[1])]
        names = block_names or x.block_names
    else:
        ref = x[0].U
        panels = [val(f)[:, match_dims(f.U, ref)[dim - 1]] for f in x]
        titles = [f"{_tr('fit')} {i + 1} - {_tr('dimension').lower()} {dim}" for i in range(len(x))]
        names = block_names or x[0].block_names
    K = len(panels[0])
    colors = _block_colors(K)
    dev2 = [(v - v.mean()) ** 2 for v in panels]
    vmax = max(max(v.max() for v in panels), 1e-12)
    dmax = max(max(v.max() for v in dev2), 1e-12)
    fig, axs = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4.6), squeeze=False)
    for a, v, dv, t in zip(axs[0], panels, dev2, titles):
        _rose(a, dv, dmax, 0, 0.32, colors, label=labels)
        _rose(a, v, vmax, 0.40, 0.72, colors, label=labels, names=names)
        m, V = v.mean(), dv.mean()
        frac = m ** 2 / (m ** 2 + V) if m ** 2 + V > 0 else 1
        a.add_patch(Wedge((0, 0), 0.92, 0, 360 * frac, width=0.08, color="#d97706", lw=0))
        if frac < 1:
            a.add_patch(Wedge((0, 0), 0.92, 360 * frac, 360, width=0.08, color="#0f766e", lw=0))
        a.text(1.08 * np.cos(np.pi * frac), 1.08 * np.sin(np.pi * frac), f"{100 * frac:.0f}%", color="#d97706",
               fontweight="bold", fontsize=8, ha="center")
        a.set_xlim(-1.45, 1.45)
        a.set_ylim(-1.45, 1.45)
        a.set_aspect("equal")
        a.axis("off")
        a.set_title(f"{t}\n{_tr('rose_mean')} = {m:.3f}  V = {V:.3f}", fontsize=9)
    sym = "d" if stat == "d" else "d^2"
    if legend:
        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(color=c, label=n) for c, n in zip(colors, names)], loc="lower center",
                   ncol=min(K, 6), frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.04))
    fig.text(0.5, 0.005, _tr("rose_bottom") % (sym, sym), ha="center", fontsize=7)
    if title:
        fig.suptitle(title)
    return axs


def plot_path(x, ax=None):
    """EXPERIMENTAL choice of R: min_k p_kh (%) by model and dimension; retained dimensions in bold."""
    Rmax = x.Rmax
    M = np.full((Rmax, Rmax), np.nan)
    Ret = np.zeros((Rmax, Rmax), bool)
    for r in range(1, Rmax + 1):
        t = x.tests[r]
        M[r - 1, :r] = t.p.min(axis=0)
        Ret[r - 1, :r] = t.retained
    txt = [["" if np.isnan(v) else f"{100 * v:.1f}" for v in row] for row in M]
    Z = -np.log10(M)
    ax = _heat(Z, txt, _tr("experimental_title") + "\n" + _tr("path_title"),
               [f"r = {r}" for r in range(1, Rmax + 1)], [f"h = {h}" for h in range(1, Rmax + 1)],
               bold=Ret, vmin=0, ax=ax)
    ax.axhline(x.R - 0.5, ls="--", color="black")
    ax.set_xlabel(_tr("dim_h"))
    ax.set_ylabel(_tr("n_dims"))
    return ax


def plot_multistart(x, ax=None):
    """Criterion reached by each random start (sorted); dashed line: deterministic CPCA start."""
    ax = _ax(ax, figsize=(6, 3.5))
    c = np.sort(x.criteria)[::-1]
    ax.plot(np.arange(1, len(c) + 1), c, "o", ms=3, color=PALETTE[0])
    ax.axhline(x.deterministic, ls="--", color=PALETTE[1])
    ax.set_xlabel(_tr("start"))
    ax.set_ylabel(_tr("criterion_lab"))
    ax.set_title(_tr("multistart_title"), fontsize=10)
    return ax
