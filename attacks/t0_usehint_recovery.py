#!/usr/bin/env python3
"""
ML-DSA t0 Recovery via UseHint Interval Oracle

ORACLE (from public signature data only):
  V[i] = (Az - ct1·2^D)[i]   ← computable from (z, A, t1, c) — all public
  Constraint: (c * t0)[i] ≈ LowBits(V[i], 2γ2) ± (γ2 - β)

  This is linear in t0 via negacyclic convolution with known c.
  Stack m signatures → m·n linear interval constraints on n unknowns.
  Solve by batched exhaustive elimination (n=4, 8^4=4096 candidates).

PARAMETERS: q=241, α=40, γ2=20, D=3 (t0 ∈ [-4,3])
  α | (q-1): 40 | 240 ✓   (required for FIPS 204 Decompose to be well-defined)
  γ2=20 > max_ct0=16: check-4 never fires → acceptance ≈ 30%
  Oracle informative: LowBits(V) offsets expose ct0 intervals (see analysis)

FINDING: check-3 binary oracle is insufficient. UseHint is the real oracle.
"""

import numpy as np
import time

# ---------------------------------------------------------------------------
# Parameters (satisfy: alpha | (q-1) and alpha_key | (q-1))
# ---------------------------------------------------------------------------
N, Q, D, G2, TAU, ETA, BETA, G1 = 4, 241, 3, 20, 4, 1, 4, 60
ALPHA     = 2 * G2                 # 40  — hint decomposition (must divide q-1=240)
ALPHA_KEY = 2 ** D                 # 8   — key decomposition  (must divide q-1=240)
HALF      = ALPHA_KEY // 2          # 4   → t0 ∈ [-(HALF-1), HALF] = [-3, 4]
M_HINT    = (Q - 1) // ALPHA       # 6   = number of hint rounding levels
M_KEY     = (Q - 1) // ALPHA_KEY   # 30  = number of key rounding levels

# ---------------------------------------------------------------------------
# Ring arithmetic
# ---------------------------------------------------------------------------

def _neg_mat(a: np.ndarray, q: int) -> np.ndarray:
    """n×n negacyclic convolution matrix for poly a mod (X^n+1, q)."""
    n    = len(a)
    ridx = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    sign = np.where(np.arange(n)[:, None] >= np.arange(n)[None, :], 1, -1)
    return (sign * a[ridx]) % q

def poly_mul(a, b, q):
    return _neg_mat(np.asarray(a, np.int64), q) @ np.asarray(b, np.int64) % q

def cmod(x, q):
    x = np.asarray(x, np.int64) % q
    return np.where(x > q // 2, x - q, x)

# ---------------------------------------------------------------------------
# Decompose / HighBits / LowBits (FIPS 204 — requires alpha | (q-1))
# ---------------------------------------------------------------------------

def decompose(r, alpha, q):
    """
    FIPS 204 Decompose: r = r1*alpha + r0, r0 ∈ (-alpha/2, alpha/2].
    Special case: if r - r0 = q-1, set r1=0, r0=r0-1.
    """
    m  = (q - 1) // alpha
    r  = np.asarray(r, np.int64) % q
    r0 = r % alpha
    r0 = np.where(r0 > alpha // 2, r0 - alpha, r0)
    r1 = (r - r0) // alpha
    # Boundary: r - r0 ≡ q-1 (mod q) → r1 = m → wrap to (0, r0-1)
    boundary = (r - r0) % q == q - 1
    r0 = np.where(boundary, r0 - 1, r0)
    r1 = np.where(boundary, 0, r1)
    return r1, r0

def high_bits(r, alpha, q): return decompose(r, alpha, q)[0]
def low_bits(r, alpha, q):  return decompose(r, alpha, q)[1]

def make_hint(z, r, alpha, q):
    """FIPS 204 MakeHint: 1 iff HighBits(r) ≠ HighBits(r+z)."""
    return (high_bits(r % q, alpha, q) != high_bits((r + z) % q, alpha, q)).astype(np.int64)

def use_hint(h, r, alpha, q):
    """FIPS 204 UseHint: recover HighBits(r+z) from (h, r) where h=MakeHint(z,r-z)."""
    m  = (q - 1) // alpha
    r1, r0 = decompose(r, alpha, q)
    # h=1: r1 is off by ±1 depending on sign of r0
    adj = np.where(r0 > 0, 1, -1)
    return np.where(h == 1, (r1 + adj) % m, r1)

# ---------------------------------------------------------------------------
# Toy ML-DSA keygen / sign / verify
# ---------------------------------------------------------------------------

def keygen(rng):
    A  = rng.integers(0, Q, N, dtype=np.int64)
    s1 = rng.integers(-ETA, ETA + 1, N, dtype=np.int64)
    s2 = rng.integers(-ETA, ETA + 1, N, dtype=np.int64)
    t  = (poly_mul(A, s1, Q) + s2) % Q
    t1 = high_bits(t, ALPHA_KEY, Q)
    t0 = low_bits(t, ALPHA_KEY, Q)
    return (A, t1), (A, s1, s2, t0, t1)

def sparse_c(rng):
    c = np.zeros(N, dtype=np.int64)
    c[rng.choice(N, TAU, replace=False)] = rng.choice([-1, 1], TAU)
    return c

def sign(sk, rng, max_tries=1000):
    A, s1, s2, t0, t1 = sk
    for _ in range(max_tries):
        y    = rng.integers(-G1 + 1, G1, N, dtype=np.int64)
        w    = poly_mul(A, y, Q)
        c    = sparse_c(rng)
        z    = cmod(y + poly_mul(c, s1, Q), Q)
        if np.max(np.abs(z)) >= G1 - BETA:
            continue
        cs2   = poly_mul(c, s2, Q)
        w_cs2 = cmod(w - cs2, Q)
        if np.max(np.abs(low_bits(w_cs2, ALPHA, Q))) >= G2 - BETA:
            continue
        ct0 = cmod(poly_mul(c, t0, Q), Q)
        if np.max(np.abs(ct0)) >= G2:
            continue                              # rarely fires (G2 > max_ct0)
        # FIPS 204: MakeHint(z=-ct0, r=w-cs2+ct0)  [r = V]
        V  = (w_cs2 + ct0) % Q
        h  = make_hint(-ct0, V, ALPHA, Q)
        w1 = high_bits(w, ALPHA, Q)
        return c, z, h, w1
    return None

def verify(pk, sig):
    A, t1 = pk
    c, z, h, w1_claimed = sig
    Az     = poly_mul(A, z, Q)
    ct1_2D = poly_mul(c, t1 * ALPHA_KEY, Q)
    V      = (Az - ct1_2D) % Q
    return np.array_equal(use_hint(h, V, ALPHA, Q), w1_claimed)

# ---------------------------------------------------------------------------
# Oracle: extract interval constraint on ct0 from public sig data
# ---------------------------------------------------------------------------

def get_oracle(pk, sig):
    """
    Return (M_c, center) such that M_c @ t0 ≈ center ± (G2-BETA).

    V[i] = w_cs2[i] + ct0[i]  (from public sig data)
    UseHint(h[i], V[i]) = HighBits(w_cs2[i]) = r1[i]
    center[i] = V[i] - r1[i]*ALPHA = LowBits(w_cs2[i]) + ct0[i]
    → ct0[i] = center[i] - LowBits(w_cs2[i])
    → |ct0[i] - center[i]| = |LowBits(w_cs2[i])| < G2-BETA  (from check-3)
    """
    A, t1  = pk
    c, z, h, _ = sig
    Az     = poly_mul(A, z, Q)
    ct1_2D = poly_mul(c, t1 * ALPHA_KEY, Q)
    V      = (Az - ct1_2D) % Q
    r1     = use_hint(h, V, ALPHA, Q)          # HighBits(w_cs2) — corrected via hint
    center = cmod(V.astype(np.int64) - r1 * ALPHA, Q)   # = LowBits(w_cs2) + ct0
    M_c    = _neg_mat(c, Q)                    # n×n: row i → (c*t0)[i]
    return M_c, center

# ---------------------------------------------------------------------------
# Batched exhaustive elimination (all 4096 candidates simultaneously)
# ---------------------------------------------------------------------------

def build_candidates():
    # t0 ∈ [-(HALF-1), HALF] via FIPS 204 Decompose (α//2 included, -(α//2)+1 min)
    r = np.arange(-(HALF - 1), HALF + 1, dtype=np.int64)  # [-3, -2, -1, 0, 1, 2, 3, 4]
    grids = np.meshgrid(*[r] * N, indexing='ij')
    return np.stack([g.ravel() for g in grids], axis=1)   # (8^N, N)

def apply_constraint(mask, cands, M_c, center, hw):
    alive  = cands[mask]
    ct0    = (alive @ M_c.T) % Q               # (k, N): ct0 = M_c @ t0 for each candidate
    ct0    = np.where(ct0 > Q // 2, ct0 - Q, ct0)
    diff   = np.abs(cmod(ct0 - center[None, :], Q))
    new_mask          = mask.copy()
    new_mask[mask]    = np.max(diff, axis=1) < hw
    return new_mask

# ---------------------------------------------------------------------------
# Main attack
# ---------------------------------------------------------------------------

def run_attack(pk, sigs, t0_true, verbose=True):
    t_start = time.time()
    cands = build_candidates()
    mask  = np.ones(len(cands), dtype=bool)
    hw    = G2 - BETA                          # interval half-width

    print(f"\n[ATTACK] {len(cands)} candidates, {len(sigs)} signatures")
    print(f"  γ2={G2}, β={BETA}, half-width={hw}")
    print(f"  Oracle: γ2 > max_ct0={TAU*HALF} — informative via LowBits offset")

    print(f"  True t0: {t0_true.tolist()}")
    for i, sig in enumerate(sigs):
        M_c, center = get_oracle(pk, sig)
        mask = apply_constraint(mask, cands, M_c, center, hw)
        remaining = mask.sum()
        if verbose and (i < 5 or (i + 1) % 5 == 0 or remaining <= 3):
            print(f"  sig {i+1:3d}: {remaining:5d} candidates remaining")
        if remaining <= 1:
            break

    survivors  = cands[mask]
    found_true = any(np.array_equal(s, t0_true) for s in survivors)
    elapsed    = time.time() - t_start

    print(f"\n  Survivors: {len(survivors)}  |  elapsed: {elapsed:.2f}s")
    print(f"  True t0 present: {found_true}")
    for s in survivors[:5]:
        marker = " ← TRUE" if np.array_equal(s, t0_true) else ""
        print(f"    {s.tolist()}{marker}")

    if len(survivors) == 1 and found_true:
        print(f"\n  *** EXACT t0 RECOVERY via UseHint oracle ✓ ***")
    elif found_true and len(survivors) <= 3:
        print(f"\n  *** NEAR-UNIQUE — trivial brute force from here ***")
    return survivors

# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 65)
    print("ML-DSA t0 Recovery — UseHint Oracle")
    print("=" * 65)
    print(f"n={N}, q={Q}, D={D} (t0∈[-{HALF},{HALF-1}])")
    print(f"γ2={G2}, τ={TAU}, η={ETA}, β={BETA}, γ1={G1}")
    print(f"α={ALPHA} | (q-1)={Q-1}: {(Q-1) % ALPHA == 0} ✓")
    print(f"max_ct0 = τ·{HALF} = {TAU*HALF},  γ2/max_ct0 = {G2/(TAU*HALF):.2f}")
    print(f"Total candidates: {(2*HALF)**N}")

    rng = np.random.default_rng(42)
    pk, sk = keygen(rng)
    _, s1, s2, t0_true, _ = sk
    print(f"\nTrue t0: {t0_true.tolist()}")

    # Collect signatures
    print("\n[SIGN] Collecting signatures...")
    sigs, attempts = [], 0
    while len(sigs) < 100 and attempts < 5000:
        attempts += 1
        sig = sign(sk, rng)
        if sig is not None:
            assert verify(pk, sig), f"verify failed at sig {len(sigs)}"
            sigs.append(sig)

    rate = len(sigs) / max(attempts, 1)
    print(f"  {len(sigs)} signatures in {attempts} attempts ({rate:.1%})")

    run_attack(pk, sigs, t0_true)

    print("\n" + "=" * 65)
    print("ML-DSA-44 EXTRAPOLATION")
    print("=" * 65)
    q44, n44, d44 = 8380417, 256, 13
    g2_44 = (q44 - 1) // 88
    tau44, beta44 = 39, 78
    max_ct0_44 = tau44 * (2 ** (d44 - 1))
    print(f"  n={n44}, q={q44}, d={d44}, γ2={g2_44}, τ={tau44}")
    print(f"  max_ct0={max_ct0_44},  γ2/max_ct0={g2_44/max_ct0_44:.3f}")
    print(f"  Candidates: 2^(n·d) = 2^{n44*d44} (intractable exhaustive)")
    print(f"  Need: LLL/BDD on {n44}×{n44} lattice — Phase 2 claim")
    print(f"  PRISM-DSA FIS removes timing oracle; UseHint oracle remains")
    print(f"  BDD hardness at δ={g2_44/max_ct0_44:.3f}: open problem")
