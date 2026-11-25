"""LOBO permutation test (article BB): which blocks are linked to the reference components?"""
import numpy as np
import pandas as pd
from ._core import cf_fit, weight_blocks, nstart_default, compress, match_dims
from ._rng import Generator
from .fit import _as_blocks, canofac
from .lang import _tr, _link_label
from .moments import gamma_pvalue, is_centred
from .preprocess import cf_preprocess


def holm(p):
    """Holm-adjusted p-values (same as R's p.adjust(method = "holm"))."""
    p = np.asarray(p, dtype=float)
    K = p.size
    o = np.argsort(p, kind="stable")
    adj = np.minimum(1.0, np.maximum.accumulate((K - np.arange(K)) * p[o]))
    out = np.empty(K)
    out[o] = adj
    return out


def _link_code(all_linked, n_linked):
    all_linked, n_linked = np.atleast_1d(all_linked), np.atleast_1d(n_linked)
    return np.where(all_linked, "all", np.where(n_linked >= 2, "several",
                                                np.where(n_linked == 1, "one", "none")))


def _group_index(groups, n):
    groups = np.asarray(groups)
    if groups.shape[0] != n:
        raise ValueError("groups must have one value per row")
    levels = list(dict.fromkeys(groups.tolist()))     # order of first appearance
    return [np.where(groups == g)[0] for g in levels]


class _Permuter:
    """Draws the permuted versions of a block from the common generator (same order as in R)."""

    def __init__(self, scheme, gen, n, groups=None):
        self.scheme, self.gen, self.n = scheme, gen, n
        self.idx = _group_index(groups, n) if scheme in ("groups", "within") else None
        if scheme == "groups":
            sizes = {len(i) for i in self.idx}
            if len(sizes) != 1:
                raise ValueError("scheme = 'groups' needs groups of equal sizes")

    def row_perms(self, B):
        """(B, n) row permutations for the row-type schemes."""
        n = self.n
        if self.scheme == "rows":
            return self.gen.permutations(B, n)
        perm = np.tile(np.arange(n), (B, 1))
        if self.scheme == "groups":
            gp = self.gen.permutations(B, len(self.idx))
            for g, ig in enumerate(self.idx):
                perm[:, ig] = np.array([self.idx[j] for j in range(len(self.idx))])[gp[:, g]]
            return perm
        for ig in self.idx:                              # within
            pg = self.gen.permutations(B, len(ig))
            perm[:, ig] = ig[pg]
        return perm


class CfTest:
    """Result of :func:`cf_test`.

    Attributes
    ----------
    p : ndarray (K, R)
        p-values by block and component of the other blocks (LOBO; descriptive) or
        observed dimension (re-estimation).
    p_holm : ndarray (K, R)
        Holm-adjusted over the blocks.
    p_subspace : ndarray (K,)
        Subspace p-values (restricted RV T_k^w by default, or T_k).
    p_subspace_equal, p_subspace_inertia : ndarray (K,)
        Both subspace statistics (LOBO).
    linked, link, linked_subspace, link_subspace, all_linked, retained
        Decisions: all blocks linked if max_k p_k <= alpha, otherwise blocks linked
        after Holm; a component is linked if min_k p_kh <= alpha / K.
    p_compare : pandas.DataFrame or None
        Full RV and PROTEST p-values (``compare=True``), same permutations.
    pmethod : ndarray (K, R + 1) of str
        "perm" or "gamma" for each p-value (components, then subspace).
    null : dict or None
        Permutation distributions of the subspace statistic, by block.
    fit : Canofac
        Fit on all the blocks.
    settings : dict
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def table(self, what="p"):
        M = getattr(self, what)
        if M.ndim == 1:
            return pd.Series(M, index=self.block_names)
        lab = "Comp" if self.settings["statistic"] == "lobo" else "Dim"
        return pd.DataFrame(M, index=self.block_names, columns=[f"{lab}{h + 1}" for h in range(M.shape[1])])

    def __repr__(self):
        s = self.settings
        fixed = s["statistic"] == "lobo"
        out = [_tr("title_test")]
        out.append((_tr("settings_test") % (s["q"], s["axes"], s["extraction"], s["scheme"], s["B"],
                                            s["alpha"])).rstrip("\n"))
        out.append(_tr("stat_fixed" if fixed else "stat_refit"))
        if fixed:
            out.append(_tr("sub_inertia" if s["subspace"] == "inertia" else "sub_equal"))
        if s["weighting"] != "none":
            out.append(_tr("weighting_line") % s["weighting"].upper())
        if s["scheme"] in ("groups", "within"):
            out.append(_tr("groups_line") % (s["n_groups"], s["scheme"]))
        if s["pvalue"] == "permutation":
            out.append(_tr("pvalue_perm") % s["B"])
        elif s["pvalue"] == "gamma":
            out.append(_tr("pvalue_gamma"))
        else:
            npm = int(np.sum(self.pmethod == "perm"))
            out.append(_tr("pvalue_hybrid") % (s["B"], s["hybrid_c"], s["hybrid_c"], npm, self.pmethod.size))
        out.append("")
        out.append(_tr("psub_head"))
        out.append((100 * self.table("p_subspace")).round(1).to_string())
        who = [b for b, l in zip(self.block_names, self.linked_subspace) if l]
        out.append(_tr("subspace_line") % (_link_label([self.link_subspace])[0],
                                           f" ({', '.join(who)})" if who else ""))
        out.append("")
        out.append(_tr("pval_head_fixed" if fixed else "pval_head"))
        out.append((100 * self.table("p")).round(1).to_string())
        out.append("")
        out.append(_tr("struct_head_fixed" if fixed else "struct_head"))
        lab = _link_label(self.link)
        w = max(len(x) for x in lab)
        cols = self.table("p").columns
        for h, c in enumerate(cols):
            who = [b for b, l in zip(self.block_names, self.linked[:, h]) if l]
            out.append(f"  {c}: {lab[h]:<{w}} " + (f"({', '.join(who)})" if who else ""))
        yn = [_tr("yes") if r else _tr("no") for r in self.retained]
        out.append("")
        out.append(_tr("linked_line" if fixed else "retained_line") +
                   ", ".join(f"{c}={y}" for c, y in zip(cols, yn)))
        if self.p_compare is not None:
            out.append("")
            out.append(_tr("compare_head"))
            out.append((100 * self.p_compare).round(1).to_string())
            out.append(_tr("compare_note"))
        return "\n".join(out)

    def plot(self, type="pvalues", **kw):
        """``type`` = "pvalues" (heat map) or "null" (permutation distributions of the subspace statistic)."""
        from . import graphics
        if type == "null":
            return graphics.plot_null(self, **kw)
        return graphics.plot_test(self, **kw)


def cf_test(blocks, R=2, q=2, axes="free", extraction="simultaneous", scheme="rows", groups=None,
            B=999, alpha=0.05, statistic="lobo", matching=True, seed=None, tol=1e-10, maxit=5000,
            nstart=None, subspace="inertia", weighting="none", pvalue="permutation", hybrid_c=3.0,
            compare=False):
    """LOBO permutation test of the links between each block and the reference components.

    For each block k, the reference components U^(-k) are estimated once on the
    other K - 1 blocks by the same procedure (leave-one-block-out) and kept
    fixed; the link between block k and U^(-k) is compared with the links
    obtained after permuting the rows of block k. Statistics: S_kh = ||X_k' u_h^(-k)||
    for each component (descriptive), and, for the subspace, the restricted RV
    T_k^w = tr(W_k P W_(-k) P) (default, ``subspace="inertia"``) or T_k^2 = tr(W_k P)
    (``subspace="equal"``). Decisions: all blocks linked if max_k p_k <= alpha
    (intersection-union test); otherwise the blocks linked after Holm adjustment;
    a component is linked if min_k p_kh <= alpha / K.

    Parameters
    ----------
    blocks : dict, list or array
        Preprocessed blocks (centred columns).
    R, q, axes, extraction, tol, maxit, nstart, weighting
        Procedure M, as in :func:`canofac`.
    scheme : {"rows", "groups", "within", "variables", "components"}
        Permutation of the rows of the tested block; of whole groups of rows
        (``groups`` of equal sizes, e.g. replicates); of the rows within each
        group; of each variable independently (H0: the block is noise); or of
        each principal component of the block independently.
    groups : array-like, optional
        Group of each row, for ``scheme="groups"`` or ``"within"``.
    B : int
        Number of permutations (B >= K / alpha - 1 at least; 10,000 to 20,000
        near a corrected threshold).
    alpha : float
        Level.
    statistic : {"lobo", "refit"}
        LOBO statistic (default; ``"fixed"`` is an alias) or re-estimation of the
        whole model on each permuted data set (comparison; less powerful).
    matching : bool
        For ``"refit"``: match the dimensions to the solution on the other blocks.
    seed : int, optional
        Seed of the common generator: the same seed gives the same permutations,
        hence the same p-values, in R and in Python.
    subspace : {"inertia", "equal"}
        Subspace statistic of the LOBO test.
    pvalue : {"permutation", "gamma", "hybrid"}
        Permutation p-values (default); gamma approximation with the exact
        permutation mean and variance, without permutation (slightly liberal);
        or hybrid: gamma, replaced by the permutation p-value when it lies
        between alpha / (hybrid_c K) and hybrid_c alpha. "gamma" and "hybrid"
        need ``statistic="lobo"`` and ``scheme="rows"``.
    hybrid_c : float
        Width of the zone of the hybrid mode (3 by default).
    compare : bool
        Also compute the full RV test (block against the concatenated other
        blocks) and PROTEST (first R principal components), with the same
        permutations. They are comparisons, not part of the LOBO test.

    Returns
    -------
    CfTest
        See :class:`CfTest`. French help: ``cf_help("cf_test", "fr")``.

    References
    ----------
    Hanafi M., Bougeard S., Galharret J.-M., Van Deun K. Which blocks are linked to
    the reference components? A leave-one-block-out permutation test for
    multiblock dimension reduction (article BB).
    """
    if statistic == "fixed":
        statistic = "lobo"
    if statistic not in ("lobo", "refit"):
        raise ValueError("statistic must be 'lobo' or 'refit'")
    if scheme not in ("rows", "groups", "within", "variables", "components"):
        raise ValueError("unknown scheme")
    if subspace not in ("inertia", "equal"):
        raise ValueError("subspace must be 'inertia' or 'equal'")
    if pvalue not in ("permutation", "gamma", "hybrid"):
        raise ValueError("pvalue must be 'permutation', 'gamma' or 'hybrid'")
    if pvalue != "permutation" and (statistic != "lobo" or scheme != "rows"):
        raise ValueError("pvalue = 'gamma' or 'hybrid' needs statistic = 'lobo' and scheme = 'rows'")
    if scheme in ("groups", "within") and groups is None:
        raise ValueError("scheme = 'groups' or 'within' needs groups")
    if compare and scheme in ("variables", "components"):
        raise ValueError("compare = True needs a row-type scheme ('rows', 'groups', 'within')")
    if statistic == "refit":
        subspace = "equal"
    names, arrs, _, _ = _as_blocks(blocks)
    K, n = len(arrs), arrs[0].shape[0]
    nstart = nstart_default(nstart, q, axes)
    if seed is None:
        seed = int(np.random.default_rng().integers(1, 2 ** 31 - 2))
    gen = Generator(seed)
    if pvalue != "permutation" and not all(is_centred(X) for X in arrs):
        raise ValueError("the gamma approximation needs column-centred blocks (see cf_preprocess)")
    work = arrs
    if scheme in ("rows", "groups", "within", "components"):
        work = [compress(X) if X.shape[1] > n - 1 else X for X in arrs]
    pcs = None
    if scheme == "components":
        pcs = []
        for X in work:
            u, s, vt = np.linalg.svd(X, full_matrices=False)
            pcs.append((u * s, vt))
    permuter = _Permuter(scheme, gen, n, groups)

    def permuted_blocks(k, Bk):
        """List (or generator) of the Bk permuted versions of block k, drawn in the R order."""
        if scheme in ("rows", "groups", "within"):
            P = permuter.row_perms(Bk)
            return P, None
        if scheme == "components":
            sc, vt = pcs[k]
            r = sc.shape[1]
            P = gen.permutations(Bk * r, n).reshape(Bk, r, n)
            return None, [np.column_stack([sc[P[b, j], j] for j in range(r)]) @ vt for b in range(Bk)]
        X = arrs[k]
        p = X.shape[1]
        P = gen.permutations(Bk * p, n).reshape(Bk, p, n)
        return None, [list(cf_preprocess([np.column_stack([X[P[b, j], j] for j in range(p)])]).values())[0]
                      for b in range(Bk)]

    Pm = np.ones((K, R))
    Psub = np.ones(K)
    Psub_eq = Psub_in = None
    pmethod = np.full((K, R + 1), "perm", dtype=object)
    null = {} if pvalue != "gamma" else None
    refs = {}
    comp_p = np.ones((K, 2)) if compare else None

    if statistic == "lobo":
        Psub_eq, Psub_in = np.ones(K), np.ones(K)
        for k in range(K):
            Xk = work[k]
            if K > 1:
                others, _ = weight_blocks([work[j] for j in range(K) if j != k], weighting)
                Uref = cf_fit(others, R, q, axes, extraction, tol, maxit, nstart)["U"]
                Cw = sum((Uref.T @ X) @ (Uref.T @ X).T for X in others)
            else:
                Uref, Cw = np.zeros((n, R)), np.zeros((R, R))
            refs[names[k]] = Uref
            A = Xk.T @ Uref
            AtA = A.T @ A
            s_obs = np.sqrt(np.diag(AtA))
            t_obs = np.sqrt(np.sum(s_obs ** 2))
            tw_obs = np.sum(AtA * Cw)
            Xr = np.hstack([work[j] for j in range(K) if j != k]) if (compare and K > 1) else None
            # gamma p-values
            need = np.ones(R + 1, dtype=bool)
            pg = None
            if pvalue != "permutation":
                Wk = Xk @ Xk.T
                pg = np.array([gamma_pvalue(s_obs[h] ** 2, Wk, np.outer(Uref[:, h], Uref[:, h])) for h in range(R)]
                              + [gamma_pvalue(t_obs ** 2, Wk, Uref @ Uref.T),
                                 gamma_pvalue(tw_obs, Wk, Uref @ Cw @ Uref.T)])
                lo, hi = alpha / (hybrid_c * K), hybrid_c * alpha
                if pvalue == "gamma":
                    need = np.zeros(R + 1, dtype=bool)
                else:
                    zone = (pg >= lo) & (pg <= hi)
                    sub_idx = R if subspace == "equal" else R + 1
                    need = np.append(zone[:R], zone[sub_idx])
                    need_eq_in = zone[R:]
            run_perm = pvalue == "permutation" or need.any() or (compare and pvalue != "gamma")
            if run_perm:
                P, _ = permuted_blocks(k, B) if scheme in ("rows", "groups", "within") else (None, None)
                if P is not None:
                    Ap = np.empty((B, Xk.shape[1], R))
                    for i in range(0, B, 2000):
                        Ap[i:i + 2000] = np.einsum("bnp,nr->bpr", Xk[P[i:i + 2000]], Uref)
                else:
                    _, Xs = permuted_blocks(k, B)
                    Ap = np.array([Xp.T @ Uref for Xp in Xs])
                Mp = np.einsum("bpr,bps->brs", Ap, Ap)
                s = np.sqrt(np.einsum("brr->br", Mp))
                tt = np.sqrt(np.sum(s ** 2, axis=1))
                tw = np.einsum("brs,rs->b", Mp, Cw)
                p_comp = (1 + np.sum(s >= s_obs - 1e-12, axis=0)) / (B + 1)
                p_eq = (1 + np.sum(tt >= t_obs - 1e-12)) / (B + 1)
                p_in = (1 + np.sum(tw >= tw_obs - 1e-12 * abs(tw_obs))) / (B + 1)
                null[names[k]] = {"observed": tw_obs if subspace == "inertia" else t_obs,
                                  "permuted": tw if subspace == "inertia" else tt}
                if compare and P is not None and K > 1:
                    Wk = Xk @ Xk.T
                    Wr = Xr @ Xr.T
                    rv_obs = np.sum(Wk * Wr)
                    rv_null = np.concatenate([np.einsum("bij,ij->b", Wk[pp[:, :, None], pp[:, None, :]], Wr)
                                              for pp in (P[i:i + 500] for i in range(0, B, 500))])
                    Zk, Zr = _pcs(Xk, R), _pcs(Xr, R)
                    pr_obs = np.linalg.svd(Zk.T @ Zr, compute_uv=False).sum()
                    pr_null = np.array([np.linalg.svd(Zk[pp].T @ Zr, compute_uv=False).sum() for pp in P])
                    comp_p[k] = [(1 + np.sum(rv_null >= rv_obs * (1 - 1e-12))) / (B + 1),
                                 (1 + np.sum(pr_null >= pr_obs - 1e-12)) / (B + 1)]
            if compare and pvalue == "gamma" and K > 1:
                Wk = Xk @ Xk.T
                Wr = Xr @ Xr.T
                comp_p[k] = [gamma_pvalue(np.sum(Wk * Wr), Wk, Wr), np.nan]
            if pvalue == "permutation":
                Pm[k], Psub_eq[k], Psub_in[k] = p_comp, p_eq, p_in
            else:
                Pm[k] = pg[:R]
                Psub_eq[k], Psub_in[k] = pg[R], pg[R + 1]
                pmethod[k] = "gamma"
                if pvalue == "hybrid" and run_perm:
                    for h in range(R):
                        if need[h]:
                            Pm[k, h] = p_comp[h]
                            pmethod[k, h] = "perm"
                    if need_eq_in[0]:
                        Psub_eq[k] = p_eq
                    if need_eq_in[1]:
                        Psub_in[k] = p_in
                    if need[R]:
                        pmethod[k, R] = "perm"
        Psub = Psub_in if subspace == "inertia" else Psub_eq
    else:
        def Mfit(bl):
            return cf_fit(weight_blocks(bl, weighting)[0], R, q, axes, extraction, tol, maxit, nstart)

        def substat(D, k):
            return np.sum(D[k]) if (axes == "orthonormal" and q == 1) else np.sum(D[k] ** 2)
        match_on = matching and K > 1
        obs = Mfit(work)
        for k in range(K):
            if match_on:
                Uref = Mfit([work[j] for j in range(K) if j != k])["U"]
                pi_obs = match_dims(obs["U"], Uref)
                s_obs = obs["d"][k, pi_obs]
                refs[names[k]] = Uref
            else:
                s_obs = obs["d"][k]
            sub_obs = substat(obs["d"], k)
            cnt, cnt_sub = np.zeros(R), 0
            P, Xs = permuted_blocks(k, B)
            permd = []
            for b in range(B):
                wp = list(work)
                wp[k] = work[k][P[b]] if P is not None else Xs[b]
                sol = Mfit(wp)
                sv = sol["d"][k, match_dims(sol["U"], Uref)] if match_on else sol["d"][k]
                cnt += sv >= s_obs - 1e-12
                ss = substat(sol["d"], k)
                cnt_sub += ss >= sub_obs - 1e-12
                permd.append(ss)
            p_ref = (1 + cnt) / (B + 1)
            if match_on:
                Pm[k, pi_obs] = p_ref
            else:
                Pm[k] = p_ref
            Psub[k] = (1 + cnt_sub) / (B + 1)
            null[names[k]] = {"observed": sub_obs, "permuted": np.array(permd)}
    fit = canofac(dict(zip(names, arrs)), R, q, axes, extraction, tol, maxit, nstart, weighting=weighting)
    P_holm = np.column_stack([holm(Pm[:, h]) for h in range(R)])
    linked = P_holm <= alpha
    all_linked = Pm.max(axis=0) <= alpha
    link = _link_code(all_linked, linked.sum(axis=0))
    linked_sub = holm(Psub) <= alpha
    link_sub = str(_link_code(Psub.max() <= alpha, linked_sub.sum())[0])
    retained = Pm.min(axis=0) <= alpha / K if statistic == "lobo" else link != "none"
    pc = None
    if compare:
        pc = pd.DataFrame(comp_p, index=names, columns=["RV", "PROTEST"])
    return CfTest(p=Pm, p_holm=P_holm, p_subspace=Psub, p_subspace_equal=Psub_eq,
                  p_subspace_inertia=Psub_in, linked=linked, link=link, linked_subspace=linked_sub,
                  link_subspace=link_sub, all_linked=all_linked, retained=retained, refs=refs, fit=fit,
                  p_compare=pc, pmethod=pmethod, null=null, block_names=names,
                  settings=dict(R=R, q=q, axes=axes, extraction=extraction, scheme=scheme, B=B, alpha=alpha,
                                statistic=statistic, matching=matching, nstart=nstart, subspace=subspace,
                                weighting=weighting, pvalue=pvalue, hybrid_c=hybrid_c, seed=seed,
                                n_groups=len(_group_index(groups, n)) if groups is not None else 0))


def _pcs(X, r):
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    Z = u[:, :r] * s[:r]
    return Z / np.linalg.norm(Z)
