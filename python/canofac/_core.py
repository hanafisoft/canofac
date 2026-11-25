"""Estimation core of canofac (same algorithms and same starts as the R package)."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from ._rng import Generator

START_SEED = 12345


def polar(G):
    """Orthonormal factor of the polar decomposition of G."""
    u, _, vt = np.linalg.svd(G, full_matrices=False)
    return u @ vt


def compress(X):
    """Principal coordinates of a block: same W = XX', cost independent of p."""
    u, s, _ = np.linalg.svd(X, full_matrices=False)
    r = int(np.sum(s > s[0] * 1e-12)) if s.size else 0
    return u[:, :r] * s[:r]


def d_free(blocks, U):
    """Pseudo-singular values for unconstrained axes, d_kh = ||X_k' u_h||."""
    return np.array([np.sqrt(np.sum((X.T @ U) ** 2, axis=0)) for X in blocks])


def d_axes(blocks, U, V):
    """Pseudo-singular values for given axes, d_kh = u_h' X_k v_kh."""
    return np.array([np.sum(U * (X @ Vk), axis=0) for X, Vk in zip(blocks, V)])


def match_dims(U, Uref):
    """pi such that dimension pi[j] of U is matched to dimension j of Uref (optimal assignment)."""
    C = np.abs(Uref.T @ U)
    rows, cols = linear_sum_assignment(-C)
    pi = np.empty(C.shape[0], dtype=int)
    pi[rows] = cols
    return pi


def effective_rank(X):
    """Effective rank D_2 = (tr W)^2 / tr(W^2)."""
    l = np.linalg.svd(X, compute_uv=False) ** 2
    return l.sum() ** 2 / np.sum(l ** 2)


def nstart_default(nstart, q, axes):
    if nstart is None:
        nstart = 10 if (axes == "free" and q > 2) else 1
    return max(1, int(nstart))


def _wnorm(X):
    return np.sqrt(np.sum((X @ X.T) ** 2)) if X.shape[1] > X.shape[0] else np.sqrt(np.sum((X.T @ X) ** 2))


def weight_blocks(blocks, weighting="none"):
    """STATIS or MFA weighting. Returns (weighted blocks, relative weights)."""
    K = len(blocks)
    if weighting == "none" or K == 0:
        return list(blocks), np.full(K, 1.0 / max(K, 1))
    if weighting == "mfa":
        s1 = np.array([np.linalg.svd(X, compute_uv=False)[0] for X in blocks])
        out = [X / s for X, s in zip(blocks, s1)]
        w = 1.0 / s1 ** 2
    elif weighting == "statis":
        nw = np.array([_wnorm(X) for X in blocks])
        Y = [X / np.sqrt(a) for X, a in zip(blocks, nw)]
        W = [y @ y.T for y in Y]
        C = np.array([[np.sum(a * b) for b in W] for a in W])
        a = np.abs(np.linalg.eigh(C)[1][:, -1])
        out = [np.sqrt(ak) * y for ak, y in zip(a, Y)]
        w = a / nw
    else:
        raise ValueError("weighting must be 'none', 'mfa' or 'statis'")
    return out, w / w.sum()


def random_starts(n, R, nstart, seed=START_SEED):
    """nstart - 1 random centred orthonormal starts, from the common generator."""
    if nstart <= 1:
        return []
    g = Generator(seed)
    out = []
    for _ in range(nstart - 1):
        U = g.normal(n * R).reshape((n, R), order="F")
        out.append(polar(U - U.mean(axis=0)))
    return out


def free_simultaneous(blocks, R, q, tol=1e-10, maxit=5000, U0=None, eps=1e-12):
    if U0 is None:
        U = np.linalg.svd(np.hstack(blocks), full_matrices=False)[0][:, :R]
    else:
        U = U0
    W = [X @ X.T for X in blocks]
    crit_old = -np.inf
    it = 0
    for it in range(1, maxit + 1):
        D = d_free(blocks, U)
        crit = np.sum(D ** q)
        if abs(crit - crit_old) <= tol * max(1.0, abs(crit)):
            break
        crit_old = crit
        G = np.zeros_like(U)
        for k, Wk in enumerate(W):
            G += (Wk @ U) * np.maximum(D[k], eps) ** (q - 2)
        U = polar(G)
    return U, d_free(blocks, U), it


def free_best(blocks, R, q, tol=1e-10, maxit=5000, nstart=1, return_all=False):
    """CPCA start plus random starts; best criterion kept."""
    U, D, _ = free_simultaneous(blocks, R, q, tol, maxit)
    best, cbest = (U, D), np.sum(D ** q)
    allc = [cbest]
    for U0 in random_starts(blocks[0].shape[0], R, nstart):
        U, D, _ = free_simultaneous(blocks, R, q, tol, maxit, U0=U0)
        c = np.sum(D ** q)
        allc.append(c)
        if c > cbest * (1 + 1e-9):
            best, cbest = (U, D), c
    if return_all:
        return best, np.array(allc)
    return best


def orth_simultaneous(blocks, R, q, tol, maxit, U0, V0):
    U, V = U0, [v.copy() for v in V0]
    crit_old = -np.inf
    for _ in range(maxit):
        D = d_axes(blocks, U, V)
        crit = np.sum(np.abs(D) ** q)
        if abs(crit - crit_old) <= tol * max(1.0, abs(crit)):
            break
        crit_old = crit
        Wt = np.abs(D) ** (q - 1) * np.sign(D)
        U = polar(sum((X @ Vk) * Wt[k] for k, (X, Vk) in enumerate(zip(blocks, V))))
        D = d_axes(blocks, U, V)
        Wt = np.abs(D) ** (q - 1) * np.sign(D)
        V = [polar((X.T @ U) * Wt[k]) for k, X in enumerate(blocks)]
    return U, V, d_axes(blocks, U, V)


def cf_fit(blocks, R, q, axes="free", extraction="simultaneous", tol=1e-10, maxit=5000, nstart=1):
    """Core estimation (no checks). Returns dict U, V (list), d (K x R)."""
    K = len(blocks)
    if extraction == "sequential":
        Xc = [X.copy() for X in blocks]
        n = blocks[0].shape[0]
        U = np.zeros((n, R))
        V = [np.zeros((X.shape[1], R)) for X in blocks]
        D = np.zeros((K, R))
        for h in range(R):
            u = free_best(Xc, 1, q, tol, maxit, nstart)[0][:, 0]
            for k in range(K):
                a = Xc[k].T @ u
                D[k, h] = np.sqrt(np.sum(a ** 2))
                V[k][:, h] = a / D[k, h] if D[k, h] > 0 else 0.0
            U[:, h] = u
            for k in range(K):
                if axes == "free":
                    Xc[k] = Xc[k] - np.outer(u, u @ Xc[k])
                else:
                    v = V[k][:, h]
                    Xc[k] = Xc[k] - np.outer(Xc[k] @ v, v)
        return {"U": U, "V": V, "d": D}
    if axes == "free":
        U, D = free_best(blocks, R, q, tol, maxit, nstart)
        V = [(X.T @ U) / np.where(D[k] > 0, D[k], 1.0) for k, X in enumerate(blocks)]
    else:
        start = cf_fit(blocks, R, 2, "orthonormal", "sequential", tol, maxit)
        U, V, D = orth_simultaneous(blocks, R, q, tol, maxit, start["U"], start["V"])
        if q == 1:                                   # identification constraint
            Mt = sum(U.T @ (X @ Vk) for X, Vk in zip(blocks, V))
            Mt = (Mt + Mt.T) / 2
            Q = np.linalg.eigh(Mt)[1][:, ::-1]
            U = U @ Q
            V = [Vk @ Q for Vk in V]
            D = d_axes(blocks, U, V)
        for k in range(K):                           # positive pseudo-singular values
            sg = np.where(D[k] < 0, -1.0, 1.0)
            V[k] = V[k] * sg
            D[k] = D[k] * sg
    order = np.argsort(-np.sum(np.abs(D) ** q, axis=0), kind="stable")
    return {"U": U[:, order], "V": [Vk[:, order] for Vk in V], "d": D[:, order]}
