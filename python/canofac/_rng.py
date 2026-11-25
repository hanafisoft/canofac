"""Common random number generator of canofac (identical in R and Python).

Park-Miller "minimal standard" generator, multiplier 48271, modulus 2^31 - 1
(MINSTD). All products stay below 2^53, so the R version, in double precision,
and this version, in 64-bit integers, produce the same sequence for the same
seed. Permutations are drawn by the Fisher-Yates algorithm from this sequence;
B permutations are drawn at once, permutation b using the draws
b(n - 1) + 1, ..., (b + 1)(n - 1) of the stream (skip-ahead), so that the
result does not depend on how the computation is vectorised.

The generator is used for every random draw of canofac: permutations of the
test and random starts of the estimation. It never touches the global random
number generator of the language.
"""
import numpy as np

M = 2147483647          # 2^31 - 1
A = 48271


def _mulmod(a, b):
    """(a * b) mod M for integers (or int64 arrays) below M."""
    return (np.asarray(a, dtype=np.int64) * np.asarray(b, dtype=np.int64)) % M


def _powmod(a, e):
    """a^e mod M (scalar)."""
    r, a = 1, int(a) % M
    e = int(e)
    while e > 0:
        if e & 1:
            r = (r * a) % M
        a = (a * a) % M
        e >>= 1
    return r


def seed_state(seed):
    """Initial state for an integer seed (any integer)."""
    return int(seed) % (M - 1) + 1


class Generator:
    """Stream of the common generator.

    Parameters
    ----------
    seed : int
        Seed; the same seed gives the same draws in R and in Python.
    """

    def __init__(self, seed=1):
        self.state = seed_state(seed)

    def uniform(self, size):
        """Uniform draws in (0, 1), as a float array of the given size."""
        size = int(size)
        out = np.empty(size)
        x = self.state
        for i in range(size):
            x = (x * A) % M
            out[i] = x / M
        self.state = x
        return out

    def normal(self, size):
        """Standard normal draws (Box-Muller), as a float array of the given size."""
        size = int(size)
        m = (size + 1) // 2
        u = self.uniform(2 * m).reshape(m, 2)
        r = np.sqrt(-2.0 * np.log(u[:, 0]))
        z = np.column_stack([r * np.cos(2 * np.pi * u[:, 1]), r * np.sin(2 * np.pi * u[:, 1])]).ravel()
        return z[:size]

    def permutations(self, B, n):
        """B permutations of 0, ..., n - 1, as a (B, n) integer array."""
        B, n = int(B), int(n)
        perm = np.tile(np.arange(n), (B, 1))
        if n < 2 or B < 1:
            return perm
        step = _powmod(A, n - 1)                      # multiplier for n - 1 draws
        x = np.empty(B, dtype=np.int64)               # state before the draws of permutation b
        s = self.state
        for b in range(B):
            x[b] = s
            s = (s * step) % M
        rows = np.arange(B)
        for i in range(n - 1, 0, -1):
            x = _mulmod(x, A)
            j = np.floor(x / M * (i + 1)).astype(np.int64)
            j = np.minimum(j, i)
            tmp = perm[rows, i].copy()
            perm[rows, i] = perm[rows, j]
            perm[rows, j] = tmp
        self.state = int(s)
        return perm

    def permutation(self, n):
        """One permutation of 0, ..., n - 1."""
        return self.permutations(1, n)[0]


def cf_permutations(B, n, seed=1):
    """B permutations of 0, ..., n - 1 drawn by the common generator (same as R, which returns 1, ..., n).

    Parameters
    ----------
    B, n : int
    seed : int

    Returns
    -------
    ndarray (B, n)
    """
    return Generator(seed).permutations(B, n)
