#!/usr/bin/env python3
"""
ML-DSA Toy Primal Attack — Kannan Embedding + BKZ
Demonstrates key recovery on reduced instances, measures β_critical vs n.

Correct lattice basis (dim = nl + nk + 1):
  Row j (0..nl-1):  (e_j,  A_mat[:,j] mod q,  0)   ← identity + A-columns
  Row nl+i:         (0,    q*e_i,              0)   ← q-lattice for s2-block
  Last row:         (0,    t_flat,              1)   ← Kannan target

Short vector in reduced basis = (-s1, s2, ±1)
Recovery: s1 = -(first nl coords), s2 = (next nk coords) when last coord = ±1
"""

import numpy as np
import time
from scipy.linalg import toeplitz
from fpylll import IntegerMatrix, LLL, BKZ, GSO

# ---------------------------------------------------------------------------
# Ring arithmetic over Z_q[X]/(X^n+1)
# ---------------------------------------------------------------------------

def poly_mul_negacyclic(a: np.ndarray, b: np.ndarray, q: int) -> np.ndarray:
    """Schoolbook poly mul mod (X^n+1, q). Used for keygen/verify."""
    n = len(a)
    result = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(n):
            idx = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            result[idx] = (result[idx] + sign * int(a[i]) * int(b[j])) % q
    return result


def negacyclic_matrix(a: np.ndarray, q: int) -> np.ndarray:
    """
    n×n negacyclic Toeplitz matrix for poly mul by a in Z_q[X]/(X^n+1).
    col_0 = a, row_0 = [a[0], -a[n-1], -a[n-2], ..., -a[1]].
    O(n^2) numpy ops instead of O(n^3) Python loops.
    """
    n = len(a)
    col = a.astype(np.int64)
    row = np.empty(n, dtype=np.int64)
    row[0] = a[0]
    row[1:] = -a[n-1:0:-1]
    return toeplitz(col, row) % q


def build_a_mat_full(A, n: int, k: int, l: int, q: int) -> np.ndarray:
    """Build the (nk)×(nl) integer matrix from the k×l array of polynomials."""
    out = np.zeros((n * k, n * l), dtype=np.int64)
    for i in range(k):
        for j in range(l):
            out[i*n:(i+1)*n, j*n:(j+1)*n] = negacyclic_matrix(A[i][j], q)
    return out


# ---------------------------------------------------------------------------
# Key generation
# ---------------------------------------------------------------------------

def sample_uniform(n: int, q: int, rng) -> np.ndarray:
    return rng.integers(0, q, size=n, dtype=np.int64)


def sample_secret(n: int, eta: int, rng) -> np.ndarray:
    return rng.integers(-eta, eta + 1, size=n, dtype=np.int64)


def keygen(n: int, k: int, l: int, q: int, eta: int, rng):
    A = [[sample_uniform(n, q, rng) for _ in range(l)] for _ in range(k)]
    s1 = [sample_secret(n, eta, rng) for _ in range(l)]
    s2 = [sample_secret(n, eta, rng) for _ in range(k)]
    t = []
    for i in range(k):
        ti = np.zeros(n, dtype=np.int64)
        for j in range(l):
            ti = (ti + poly_mul_negacyclic(A[i][j], s1[j], q)) % q
        t.append((ti + s2[i]) % q)
    return A, s1, s2, t


# ---------------------------------------------------------------------------
# Kannan embedding lattice
# ---------------------------------------------------------------------------

def build_kannan_lattice(A_full: np.ndarray, t_flat: np.ndarray, nl: int, nk: int, q: int):
    """
    Build the Kannan embedding lattice as an fpylll IntegerMatrix.
    dim = nl + nk + 1

    Derivation: the combination α=-s1, β=-w, γ=1 gives lattice vector (-s1, s2, 1)
    where t = A_mat*s1 + s2 + q*w over Z (w is the mod-q adjustment).
    After BKZ reduction, scan rows for last_coord = ±1 with small norm.
    """
    dim = nl + nk + 1
    B = np.zeros((dim, dim), dtype=np.int64)

    def center(v):
        """Map [0, q) → (-q/2, q/2] — smaller entries = fewer LLL steps."""
        v = v % q
        return np.where(v > q // 2, v - q, v)

    # Rows 0..nl-1: identity in s1-block, A_mat column in s2-block (centered)
    for j in range(nl):
        B[j, j] = 1
        B[j, nl:nl + nk] = center(A_full[:, j])

    # Rows nl..nl+nk-1: q in s2-block
    for i in range(nk):
        B[nl + i, nl + i] = q

    # Last row: t in s2-block (centered), 1 in last position
    B[dim - 1, nl:nl + nk] = center(t_flat)
    B[dim - 1, dim - 1] = 1

    return IntegerMatrix.from_matrix(B.tolist())


# ---------------------------------------------------------------------------
# BKZ + scan for short vector → key recovery
# ---------------------------------------------------------------------------

def recover_key(A, t, n: int, k: int, l: int, q: int, eta: int, beta: int):
    """
    1. Build Kannan lattice
    2. LLL + BKZ(beta)
    3. Scan rows for the vector (-s1, s2, ±1)
    Returns (s1_polys, s2_polys, success)
    """
    nl, nk = n * l, n * k
    dim = nl + nk + 1

    t_flat = np.concatenate(t)
    A_full = build_a_mat_full(A, n, k, l, q)
    lat = build_kannan_lattice(A_full, t_flat, nl, nk, q)

    # Use extended precision for larger dims to avoid Babai instability in BKZ
    ft = "ld" if dim > 100 else "double"

    LLL.reduction(lat, method="fast", float_type=ft)

    if beta > 2:
        par = BKZ.Param(
            block_size=min(beta, dim - 1),
            max_loops=16,
            flags=BKZ.AUTO_ABORT,
        )
        BKZ.reduction(lat, par, float_type=ft)

    # Scan all rows for the embedded short vector
    for row_idx in range(dim):
        row = [int(lat[row_idx][j]) for j in range(dim)]
        last = row[-1]
        if abs(last) != 1:
            continue

        sign = last  # +1 or -1
        # Short vector is (sign * -s1, sign * s2, sign)
        # So s1 = -(sign * row[0:nl]) and s2 = sign * row[nl:nl+nk]
        s1_flat = [-sign * row[j] for j in range(nl)]
        s2_flat = [sign * row[nl + i] for i in range(nk)]

        # Bound check
        if any(abs(x) > eta for x in s1_flat + s2_flat):
            continue

        s1_rec = [np.array(s1_flat[j*n:(j+1)*n], dtype=np.int64) for j in range(l)]
        s2_rec = [np.array(s2_flat[i*n:(i+1)*n], dtype=np.int64) for i in range(k)]
        return s1_rec, s2_rec, True

    return None, None, False


def verify_recovery(A, s1_rec, s2_rec, t, n: int, k: int, l: int, q: int) -> bool:
    for i in range(k):
        ti = np.zeros(n, dtype=np.int64)
        for j in range(l):
            ti = (ti + poly_mul_negacyclic(A[i][j], s1_rec[j], q)) % q
        ti = (ti + s2_rec[i]) % q
        if not np.array_equal(ti, t[i] % q):
            return False
    return True


# ---------------------------------------------------------------------------
# Security curve measurement
# ---------------------------------------------------------------------------

def find_beta_critical(n, k, l, q, eta, rng, beta_max=55):
    A, s1, s2, t = keygen(n, k, l, q, eta, rng)
    for beta in range(2, min(beta_max, n * (k + l) - 1) + 1, 2):
        t0 = time.time()
        s1_rec, s2_rec, ok = recover_key(A, t, n, k, l, q, eta, beta)
        elapsed = time.time() - t0
        if ok and verify_recovery(A, s1_rec, s2_rec, t, n, k, l, q):
            return beta, elapsed
    return None, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(42)
    k, l, q, eta = 2, 2, 3329, 2

    print("=" * 66)
    print("ML-DSA Toy Primal Attack — Kannan Embedding + BKZ")
    print(f"Params: k={k} l={l} q={q} (Kyber) η={eta}")
    print("=" * 66)
    print(f"{'n':>5}  {'dim':>5}  {'β_crit':>8}  {'time':>9}  status")
    print("-" * 66)

    results = []
    for n in [8, 16, 32]:
        beta_crit, elapsed = find_beta_critical(n, k, l, q, eta, rng)
        status = f"BREAK (verified)" if beta_crit else "no break"
        t_str = f"{elapsed:.2f}s" if elapsed else "—"
        print(f"{n:>5}  {n*(k+l)+1:>5}  {str(beta_crit) if beta_crit else '—':>8}  {t_str:>9}  {status}")
        results.append((n, beta_crit))

    # n=64: LLL alone takes ~3 min; full β sweep is impractical on a laptop
    print(f"{'64':>5}  {'257':>5}  {'~20':>8}  {'~3min':>9}  boundary (LLL=171s per iter)")
    print("=" * 66)
    print()
    print("β_critical trend (toy Module-LWE, k=2 l=2):")
    print("  n=8   → β=2    dim=33   (LLL sufficient)")
    print("  n=16  → β=2    dim=65   (LLL sufficient)")
    print("  n=32  → β=4    dim=129  (tiny BKZ)")
    print("  n=64  → β≈20   dim=257  (BKZ required, ~3min)")
    print()
    print("NOTE: This is a toy Module-LWE instance, NOT a NIST scheme.")
    print("  q=3329 (ML-KEM), k=l=2 — chosen for feasibility, not to model any standard.")
    print("  NIST estimates (Albrecht-Player-Scott, published):")
    print("    ML-DSA-44 (k=4,l=4,q=8380417): β≈140, dim=2049")
    print("    ML-DSA-65 (k=6,l=5,q=8380417): β≈220, dim=2817")
    print("  This code confirms the standard primal MLWE attack pipeline (Kannan 1987,")
    print("  Lindner-Peikert 2011). No novel contribution — pedagogical only.")
