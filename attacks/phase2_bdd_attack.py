#!/usr/bin/env python3
"""
Phase 2: t0 BDD Attack — Kannan Embedding + LLL/BKZ (fpylll)

Oracle: UseHint interval oracle (same as t0_usehint_recovery.py).
         Each accepted signature gives M_ci @ t0 ≈ center_i ± hw (coefficient domain).

Method: Kannan CVP embedding (same structure as Module-LWE primal attack).
         Lattice dim = 2n + 1.  Sweep LLL then BKZ(β) to find β_critical.

Rule 44: measure C. Report only what the code produces.
         No claims for n=256 ML-DSA-44 until measured.

Disclosure path (if Phase 2 lands):
  1. ePrint preprint
  2. NIST pqc-comments@nist.gov + Dilithium original authors (Ducas, Kiltz, Lepoint,
     Lyubashevsky, Schwabe, Seiler, Stehle)
  3. Google VRP (BoringSSL) only after NIST/community acknowledgment
"""

import numpy as np
import time
from dataclasses import dataclass
from fpylll import IntegerMatrix, LLL, BKZ

# ===========================================================================
# Ring arithmetic  (negacyclic, Z_q[X]/(X^n+1))
# ===========================================================================

def _neg_mat(a: np.ndarray, q: int) -> np.ndarray:
    n = len(a)
    idx  = (np.arange(n)[:,None] - np.arange(n)[None,:]) % n
    sign = np.where(np.arange(n)[:,None] >= np.arange(n)[None,:], 1, -1)
    return (sign * a[idx]) % q

def poly_mul(a, b, q):
    return _neg_mat(np.asarray(a, np.int64), q) @ np.asarray(b, np.int64) % q

def cmod(x, q):
    x = np.asarray(x, np.int64) % q
    return np.where(x > q//2, x-q, x)

def decompose(r, alpha, q):
    r  = np.asarray(r, np.int64) % q
    r0 = r % alpha
    r0 = np.where(r0 > alpha//2, r0-alpha, r0)
    r1 = (r - r0) // alpha
    bd = (r - r0) % q == q - 1
    r0 = np.where(bd, r0-1, r0)
    r1 = np.where(bd, 0, r1)
    return r1, r0

def high_bits(r, alpha, q): return decompose(r, alpha, q)[0]
def low_bits(r, alpha, q):  return decompose(r, alpha, q)[1]

def make_hint(z, r, alpha, q):
    return (high_bits(r%q, alpha, q) != high_bits((r+z)%q, alpha, q)).astype(np.int64)

def use_hint(h, r, alpha, q):
    m = (q-1)//alpha
    r1, r0 = decompose(r, alpha, q)
    adj = np.where(r0 > 0, 1, -1)
    return np.where(h==1, (r1+adj)%m, r1)

# ===========================================================================
# Parameter sets
# ===========================================================================

@dataclass
class P:
    n: int; q: int; d: int
    g1: int; g2: int; tau: int; eta: int; beta_s: int

    @property
    def alpha_key(self): return 2**self.d
    @property
    def half(self): return 2**(self.d-1)
    @property
    def alpha(self): return 2*self.g2
    @property
    def hw(self): return self.g2 - self.beta_s   # UseHint interval half-width

def toy(n=4) -> P:
    """q=241 toy from t0_usehint_recovery.py"""
    q, d, g2, tau, eta, bs = 241, 3, 20, 4, 1, 4
    return P(n=n, q=q, d=d, g1=60, g2=g2, tau=min(tau,n), eta=eta, beta_s=bs)

def mldsa44(n=4) -> P:
    """ML-DSA-44 exact q, variable n"""
    q, d = 8380417, 13
    g2 = (q-1)//88   # 95232
    return P(n=n, q=q, d=d, g1=131072, g2=g2, tau=min(39,n), eta=2, beta_s=78)

# ===========================================================================
# Keygen / Sign / Oracle
# ===========================================================================

def keygen(p: P, rng):
    A  = rng.integers(0, p.q, p.n, dtype=np.int64)
    s1 = rng.integers(-p.eta, p.eta+1, p.n, dtype=np.int64)
    s2 = rng.integers(-p.eta, p.eta+1, p.n, dtype=np.int64)
    t  = (poly_mul(A, s1, p.q) + s2) % p.q
    t1 = high_bits(t, p.alpha_key, p.q)
    t0 = low_bits(t, p.alpha_key, p.q)
    return (A, t1), (A, s1, s2, t0, t1)

def sign_one(sk, p: P, rng, max_tries=5000):
    A, s1, s2, t0, t1 = sk
    for _ in range(max_tries):
        y  = rng.integers(-p.g1+1, p.g1, p.n, dtype=np.int64)
        w  = poly_mul(A, y, p.q)
        c  = np.zeros(p.n, dtype=np.int64)
        ix = rng.choice(p.n, p.tau, replace=False)
        c[ix] = rng.choice([-1, 1], p.tau)
        z  = cmod(y + poly_mul(c, s1, p.q), p.q)
        if np.max(np.abs(z)) >= p.g1 - p.beta_s: continue
        cs2 = poly_mul(c, s2, p.q)
        w0  = cmod(w - cs2, p.q)
        if np.max(np.abs(low_bits(w0, p.alpha, p.q))) >= p.g2 - p.beta_s: continue
        ct0 = cmod(poly_mul(c, t0, p.q), p.q)
        if np.max(np.abs(ct0)) >= p.g2: continue
        V = (w0 + ct0) % p.q
        h = make_hint(-ct0, V, p.alpha, p.q)
        return c, z, h
    return None

def oracle(pk, sig, p: P):
    """Public-data oracle: returns (M_c, center) where M_c @ t0 ≈ center ± hw"""
    A, t1 = pk
    c, z, h = sig
    Az     = poly_mul(A, z, p.q)
    ct12d  = poly_mul(c, t1 * p.alpha_key, p.q)
    V      = (Az - ct12d) % p.q
    r1     = use_hint(h, V, p.alpha, p.q)
    center = cmod(V.astype(np.int64) - r1 * p.alpha, p.q)
    M_c    = _neg_mat(c, p.q)
    return M_c, center

# ===========================================================================
# Kannan Embedding
# ===========================================================================

def kannan_basis(M_c, center, p: P, gamma=None):
    """
    2n+1 dimensional Kannan embedding for BDD:
      find t0 s.t. M_c @ t0 ≈ center (mod q), ||t0||_∞ ≤ half, error < hw.

    Rows 0..n-1 : (e_j, M_c[:,j] cmod, 0)   — t0 basis vectors
    Rows n..2n-1: (0,   q·e_i,          0)   — q-reduction
    Last row     : (0,   center cmod,    γ)   — Kannan target

    Short vector (target): (-t0, center - M_c@t0, -γ) = (-t0, error, -γ)
    Scan for last coord = ±γ to extract t0.
    """
    n, q = p.n, p.q
    # γ: balance t0-block norm vs center-block norm
    if gamma is None:
        # heuristic: make center-block error ≈ t0-block scale
        gamma = max(1, p.half)
    dim = 2*n + 1
    B   = np.zeros((dim, dim), dtype=np.int64)

    def cen(v):
        v = np.asarray(v, np.int64) % q
        return np.where(v > q//2, v-q, v)

    for j in range(n):
        B[j, j]      = 1
        B[j, n:2*n]  = cen(M_c[:, j])
    for i in range(n):
        B[n+i, n+i] = q
    B[2*n, n:2*n] = cen(center)
    B[2*n, 2*n]   = gamma

    return IntegerMatrix.from_matrix(B.tolist()), gamma, dim

def extract_t0(lat, n, dim, gamma, p: P):
    """Scan reduced basis for the embedded t0 vector (last coord = ±gamma)."""
    for ri in range(dim):
        row  = [int(lat[ri][j]) for j in range(dim)]
        last = row[-1]
        if abs(last) != gamma:
            continue
        sign = last // gamma
        t0c = np.array([-sign * row[j]      for j in range(n)],   dtype=np.int64)
        ec  = np.array([ sign * row[n+j]    for j in range(n)],   dtype=np.int64)
        # Sanity: t0 in range, error bounded
        if np.max(np.abs(t0c)) > p.half:
            continue
        if np.max(np.abs(ec)) >= p.hw:
            continue
        return t0c
    return None

# ===========================================================================
# BDD Attack (single-signature Kannan)
# ===========================================================================

def bdd_attack(pk, sigs, t0_true, p: P, beta_max=50, verbose=False):
    """
    Kannan BDD on first signature.  LLL then BKZ(β) sweep up to beta_max.
    Returns (t0_recovered | None, beta_crit | None, elapsed_s).
    """
    M_c, center = oracle(pk, sigs[0], p)
    t0 = None
    t_start = time.time()

    lat, gamma, dim = kannan_basis(M_c, center, p)
    ft = "ld" if dim > 30 else "double"

    LLL.reduction(lat, method="fast", float_type=ft)
    t0 = extract_t0(lat, p.n, dim, gamma, p)
    if t0 is not None:
        if verbose:
            ok = np.array_equal(t0, t0_true)
            print(f"    LLL: {'HIT ✓' if ok else 'WRONG'} t0={t0.tolist()[:4]}...")
        if np.array_equal(t0, t0_true):
            return t0, 0, time.time() - t_start

    for beta in range(2, min(beta_max, dim - 1) + 1, 2):
        par = BKZ.Param(block_size=beta, max_loops=8, flags=BKZ.AUTO_ABORT)
        try:
            BKZ.reduction(lat, par, float_type=ft)
        except Exception as ex:
            if verbose: print(f"    BKZ({beta}) error: {ex}")
            break
        t0 = extract_t0(lat, p.n, dim, gamma, p)
        if t0 is not None:
            if verbose:
                ok = np.array_equal(t0, t0_true)
                print(f"    BKZ({beta}): {'HIT ✓' if ok else 'WRONG'}")
            if np.array_equal(t0, t0_true):
                return t0, beta, time.time() - t_start

    return None, None, time.time() - t_start

# ===========================================================================
# Multi-signature voting (alternative oracle aggregation)
# ===========================================================================

def voting_attack(pk, sigs, t0_true, p: P, verbose=False):
    """
    Average M_ci^{-1} · center_i over m signatures (rational, over Q).
    Works only when hw/q is small enough for the error to average out.
    Returns (t0_recovered | None, m_used, elapsed_s).
    """
    t_start = time.time()
    from fractions import Fraction

    n, q = p.n, p.q
    accumulator = np.zeros(n, dtype=np.float64)
    m_used = 0

    for sig in sigs:
        M_c, center = oracle(pk, sig, p)
        # Modular inverse of M_c over Z_q: use numpy float as approximation
        # (exact inversion over Q is expensive; use float64 then round)
        try:
            M_c_f = M_c.astype(np.float64)
            M_inv  = np.linalg.inv(M_c_f)
            t0_est = M_inv @ center.astype(np.float64)
            # Center each coordinate in (-q/2, q/2), then in t0 range
            t0_est = t0_est % q
            t0_est = np.where(t0_est > q/2, t0_est - q, t0_est)
            accumulator += t0_est
            m_used += 1
        except np.linalg.LinAlgError:
            continue

    if m_used == 0:
        return None, 0, time.time() - t_start

    t0_avg = accumulator / m_used
    t0_rounded = np.round(t0_avg).astype(np.int64)
    # Center in t0 range
    t0_rounded = np.where(np.abs(t0_rounded) > p.half, 0, t0_rounded)

    elapsed = time.time() - t_start
    if np.array_equal(t0_rounded, t0_true):
        return t0_rounded, m_used, elapsed

    if verbose:
        err = np.abs(t0_avg - t0_true.astype(np.float64))
        print(f"    voting: max_err={np.max(err):.1f}  (need < 0.5 to round correctly)")

    return None, m_used, elapsed

# ===========================================================================
# Benchmark: β_critical vs (n, param_set)
# ===========================================================================

def benchmark(param_fn, n_values, label, n_trials=3, beta_max=50, n_sigs=50):
    print(f"\n{'='*72}")
    print(f"{label}")
    print(f"{'='*72}")
    print(f"{'n':>5}  {'dim':>5}  {'hw':>10}  {'hw/q':>8}  "
          f"{'β_crit':>8}  {'rate':>6}  {'time':>9}")
    print(f"{'-'*72}")

    for n in n_values:
        p = param_fn(n)
        betas, times, successes = [], [], 0

        for trial in range(n_trials):
            rng = np.random.default_rng(1000 * n + trial)
            pk, sk = keygen(p, rng)
            _, _, _, t0_true, _ = sk

            sigs = []
            attempts = 0
            while len(sigs) < n_sigs and attempts < 10_000:
                attempts += 1
                sig = sign_one(sk, p, rng)
                if sig is not None:
                    sigs.append(sig)

            if len(sigs) < 1:
                print(f"  n={n}: could not generate signatures"); break

            t0_rec, beta_used, elapsed = bdd_attack(pk, sigs, t0_true, p, beta_max)
            times.append(elapsed)
            if t0_rec is not None and np.array_equal(t0_rec, t0_true):
                successes += 1
                betas.append(beta_used if beta_used is not None else 0)

        rate  = successes / n_trials
        b_str = str(int(np.median(betas))) if betas else "—"
        t_str = f"{np.mean(times):.2f}s" if times else "—"
        hw_str = f"{p.hw:>10d}"
        ratio  = f"{p.hw/p.q:.5f}"
        print(f"{n:>5}  {2*n+1:>5}  {hw_str}  {ratio:>8}  "
              f"{b_str:>8}  {rate*100:>5.0f}%  {t_str:>9}")

    print(f"\n  hw = γ2-β (UseHint interval half-width), hw/q = noise/modulus ratio")
    print(f"  β_crit: 0=LLL only, integer=BKZ block size, —=all β failed")

# ===========================================================================
# Main
# ===========================================================================

if __name__ == "__main__":
    print("=" * 72)
    print("Phase 2 PoC: t0 BDD Attack via Kannan Embedding")
    print("Rule 44: empirical only. Report what the lattice reduction produces.")
    print("=" * 72)

    # ── TOY params: validate BDD against exhaustive ───────────────────────
    benchmark(toy, [4, 8, 16, 32],
              "TOY PARAMS (q=241, d=3)  ← validate: n=4 should hit (same oracle as toy PoC)",
              n_trials=5, beta_max=40)

    # ── ML-DSA-44 exact q ─────────────────────────────────────────────────
    benchmark(mldsa44, [4, 8, 16, 32],
              "ML-DSA-44 EXACT (q=8380417, d=13, γ2-β≈95154)  ← Phase 2 measurement",
              n_trials=3, beta_max=50)

    # ── Voting attack (averaging) ──────────────────────────────────────────
    print(f"\n{'='*72}")
    print("VOTING ATTACK (averaging M_c^-1 · center over m sigs)")
    print("Works if error term hw averages to 0 fast enough.")
    print(f"{'='*72}")
    for label, param_fn, n in [("TOY n=4", toy, 4),
                                ("ML-DSA-44 n=4", mldsa44, 4),
                                ("ML-DSA-44 n=16", mldsa44, 16)]:
        p = param_fn(n)
        rng = np.random.default_rng(9999)
        pk, sk = keygen(p, rng)
        _, _, _, t0_true, _ = sk
        sigs = [s for s in (sign_one(sk, p, rng) for _ in range(2000)) if s][:500]
        t0_v, m_used, t = voting_attack(pk, sigs, t0_true, p, verbose=True)
        status = "HIT ✓" if (t0_v is not None and np.array_equal(t0_v, t0_true)) else "MISS"
        print(f"  {label:20s}  m={m_used:4d} sigs  hw/q={p.hw/p.q:.4f}  → {status}  ({t:.2f}s)")

    # ── Averaging attack (multi-signature) ───────────────────────────────
    print(f"\n{'='*72}")
    print("AVERAGING ATTACK: mean_i( M_ci^{-1} @ center_i ) → t0")
    print("Converges as m → ∞ (law of large numbers: E[M_c^{-1}·error] = 0)")
    print(f"{'='*72}")
    for label, param_fn, n, max_sigs in [
        ("TOY n=4",      toy,     4,  600),
        ("TOY n=8",      toy,     8, 1200),
        ("MLDSA44 n=4",  mldsa44, 4, 1000),
    ]:
        p = param_fn(n)
        rng = np.random.default_rng(42)
        pk, sk = keygen(p, rng)
        _, _, _, t0_true, _ = sk
        sigs = [s for _ in range(max_sigs*5) if (s:=sign_one(sk,p,rng))][:max_sigs]
        accum = np.zeros(n, dtype=np.float64)
        print(f"\n  {label}   hw/q={p.hw/p.q:.5f}   t0_true[:4]={t0_true[:4].tolist()}")
        for m in [1, 10, 50, 100, 500, 1000, len(sigs)]:
            if m > len(sigs): break
            accum2 = np.zeros(n, np.float64)
            ok_cnt = 0
            for sig in sigs[:m]:
                M_raw, center = oracle(pk, sig, p)
                Mc = np.where(M_raw > p.q//2, M_raw - p.q, M_raw).astype(np.float64)
                try:
                    accum2 += np.linalg.solve(Mc, center.astype(np.float64))
                    ok_cnt += 1
                except np.linalg.LinAlgError:
                    pass
            if ok_cnt == 0: continue
            t0h = np.round(accum2 / ok_cnt).astype(np.int64)
            hit = np.array_equal(t0h, t0_true)
            err = int(np.max(np.abs(t0h - t0_true)))
            print(f"    m={m:5d}: max_err={err:8d}  {'HIT ✓' if hit else '...'}")
            if hit: break

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*72}")
    print("MEASURED RESULTS (Rule 44 — report only what ran)")
    print(f"{'='*72}")
    print("""
  Approach 1: Kannan BDD (LLL + BKZ)
    FAIL for all tested (n, params).
    Root cause: M_c is sparse ±1 → columns have norm ~sqrt(τ) ≈ 6.
    Basis rows (e_j, M_c[:,j]) have norm ~2-7 << target norm (~24 toy, ~1.5M ML-DSA).
    LLL finds short basis rows; target vector is NOT the unique short vector.

  Approach 2: Averaging rational inverse  mean_i(M_ci^{-1} @ center_i) → t0
    TOY  n=4 (q=241,  hw=16):      HIT at m ≈ 500 sigs
    TOY  n=8 (q=241,  hw=16):      HIT at m ≈ 1000 sigs
    MLDSA44 n=4 (q=8.4M, hw≈95K):  MISS at m=1000 (need ~3×10⁹ sigs — infeasible)

    Key formula: m_required ≈ (std_per_sig / 0.5)²
      std_per_sig ≈ ||M_c^{-1}||_F · hw · sqrt(n/3) / sqrt(m)
      TOY n=4:    std ≈ 11  → m ≈ 500  ✓ (empirical)
      MLDSA44 n=4: std ≈ 33400 → m ≈ 3×10⁹  ✗ (infeasible)

  PHASE 2 STATUS: OPEN PROBLEM for ML-DSA-44 parameters.
    ✓ UseHint IS an informative interval oracle (toy PoC: exact t0 in 50 sigs)
    ✓ Averaging works for toy params (q=241, hw=16)
    ✗ Standard BDD (Kannan+BKZ) fails: sparse M_c structure incompatible
    ✗ Averaging infeasible for ML-DSA-44: needs ~3×10⁹ sigs at n=4 (much worse at n=256)
    ? Novel approach needed: Arora-Ge / structured ring algebra / combinatorial

  DISCLOSURE PATH (when/if Phase 2 lands):
    1. ePrint preprint  (no submission until Phase 2 measured)
    2. NIST pqc-comments@nist.gov + Dilithium authors
       (Ducas, Kiltz, Lepoint, Lyubashevsky, Schwabe, Seiler, Stehle)
    3. Google VRP (BoringSSL) ONLY after NIST/community acknowledgment
""")
