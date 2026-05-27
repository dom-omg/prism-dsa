#!/usr/bin/env python3
"""
Phase 5A — Nonce Reuse Key Recovery on ML-DSA
==============================================

Attack scenario: implementation omits `mu` (message hash) from rho' derivation.
Spec:  rho' = SHAKE256(key || rnd || mu)   [message-dependent, safe]
Bug:   rho' = SHAKE256(key)                [message-independent, catastrophic]

With the bug, ALL signatures share the same y = ExpandMask(rho', 0).
Given two valid signatures (z1, c1) and (z2, c2) on different messages:

  z1 = y + c1·s1  (mod R_q)
  z2 = y + c2·s1  (mod R_q)
  ──────────────────────────
  z1 - z2 = (c1 - c2)·s1   →   s1 = (z1 - z2) · inv(c1 - c2)  in R_q

Full secret key s1 recovered from 2 signatures.
This attack applies to: any implementation where rho' is constant across messages.
PRISM-DSA is NOT vulnerable: FIS derives rho' from key+rnd+mu each signing call.

CVE surface: implementations that use static/constant RNG seeds, or cache rho'
across calls without per-message freshness (e.g., buggy HSM firmware).
"""

import hashlib
import secrets
import struct
import time
from typing import Optional

# ─── Parameters ──────────────────────────────────────────────────────────────

# Toy parameters (instant, pedagogically clear)
Q_TOY = 241
N_TOY = 4
ETA_TOY = 2
TAU_TOY = 2
GAMMA1_TOY = 8
BETA_TOY = TAU_TOY * ETA_TOY  # 4

# ML-DSA-44 full parameters
Q = 8380417          # 2^23 - 2^13 + 1
N = 256
ETA = 2
TAU = 39
BETA = TAU * ETA     # 78
GAMMA1 = 1 << 17

# ─── Integer utilities ───────────────────────────────────────────────────────

def modinv(a: int, m: int) -> int:
    g, x, _ = _ext_gcd(a % m, m)
    if g != 1:
        raise ValueError(f"{a} not invertible mod {m}")
    return x % m

def _ext_gcd(a: int, b: int):
    if a == 0:
        return b, 0, 1
    g, x, y = _ext_gcd(b % a, a)
    return g, y - (b // a) * x, x

# ─── Polynomial ring R_q = Z_q[X]/(X^n+1) ───────────────────────────────────

def poly_add(a, b, q):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)): r[i] = (r[i] + a[i]) % q
    for i in range(len(b)): r[i] = (r[i] + b[i]) % q
    return r

def poly_sub(a, b, q):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)): r[i] = (r[i] + a[i]) % q
    for i in range(len(b)): r[i] = (r[i] - b[i]) % q
    return r

def poly_mul_neg(a, b, n, q):
    """Negacyclic polynomial multiplication in Z_q[X]/(X^n+1)."""
    r = [0] * n
    for i in range(n):
        for j in range(n):
            idx = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            r[idx] = (r[idx] + sign * a[i] * b[j]) % q
    return r

def poly_center(a, q):
    """Map coefficients from [0,q) to (-q/2, q/2]."""
    return [(x - q if x > q // 2 else x) for x in a]

# ─── Polynomial inversion via Extended Euclidean Algorithm ───────────────────
# Works for any n, q without NTT tables.

def _poly_degree(p):
    for i in range(len(p)-1, -1, -1):
        if p[i] != 0:
            return i
    return -1  # zero polynomial

def _poly_divmod(a, b, q):
    """Polynomial division over Z_q. Returns (quotient, remainder)."""
    a = list(a)
    b = list(b)
    deg_a = _poly_degree(a)
    deg_b = _poly_degree(b)
    if deg_b < 0:
        raise ZeroDivisionError("division by zero polynomial")
    if deg_a < deg_b:
        return [0], a

    quotient = [0] * (deg_a - deg_b + 1)
    lead_inv = modinv(b[deg_b], q)

    while True:
        deg_a = _poly_degree(a)
        if deg_a < deg_b:
            break
        factor = a[deg_a] * lead_inv % q
        shift = deg_a - deg_b
        quotient[shift] = factor
        for i in range(deg_b + 1):
            a[i + shift] = (a[i + shift] - factor * b[i]) % q

    return quotient, a

def _poly_mul_plain(a, b, q):
    """Polynomial multiplication over Z_q (no reduction mod X^n+1)."""
    if not a or not b:
        return [0]
    result = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            result[i + j] = (result[i + j] + ai * bj) % q
    return result

def _poly_sub_plain(a, b, q):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)): r[i] = (r[i] + a[i]) % q
    for i in range(len(b)): r[i] = (r[i] - b[i]) % q
    return r

def poly_inv_neg(f, n, q):
    """
    Invert f in R_q = Z_q[X]/(X^n+1) using Extended Euclidean Algorithm.
    Returns g such that f*g ≡ 1 (mod X^n+1, mod q), or None if not invertible.
    """
    # Modulus polynomial: X^n + 1
    modulus = [0] * (n + 1)
    modulus[0] = 1   # coeff of X^0
    modulus[n] = 1   # coeff of X^n

    # EEA: maintain r0*s0 ≡ f and r1*s1 ≡ f (mod modulus)
    # Initially: r0 = f, s0 = [1]; r1 = modulus, s1 = [0]
    r0 = list(f) + [0] * (n + 1 - len(f))
    r1 = list(modulus)
    s0 = [1] + [0] * n        # s0 * f ≡ r0 (mod modulus): trivially
    s1 = [0] * (n + 1)        # s1 * f ≡ r1 (mod modulus): trivially

    for _ in range(2 * (n + 1)):  # bounded iterations
        if _poly_degree(r1) < 0:
            break
        q_div, r_new = _poly_divmod(r0, r1, q)
        s_new = _poly_sub_plain(s0, _poly_mul_plain(q_div, s1, q), q)
        r0, r1 = r1, r_new
        s0, s1 = s1, s_new

    # r0 should now be the GCD
    deg = _poly_degree(r0)
    if deg > 0:
        return None  # GCD has degree > 0 → not invertible
    if r0[0] == 0:
        return None  # GCD is zero

    lead_inv = modinv(r0[0], q)
    # s0 * f ≡ r0 (mod modulus), so (s0 * lead_inv) * f ≡ 1 (mod X^n+1)
    result = [x * lead_inv % q for x in s0[:n]]
    return result

# ─── NTT for ML-DSA (n=256, q=8380417) ──────────────────────────────────────
# NTT-based inversion: O(n log n) vs O(n^2) EEA.
# For n=256 the EEA is fast enough (~65K ops), but NTT is cleaner for scale.

def _bitrev8(k: int) -> int:
    return int(format(k & 0xFF, '08b')[::-1], 2)

# ML-DSA ZETAS: ζ = 1753 is a primitive 512th root of unity mod Q (ζ^256 ≡ -1 mod Q).
# ZETAS[k] = ζ^{bitrev_8(k)} mod Q — direct powers, NOT psi = ζ^{(q-1)/512}.
_ZETAS_256 = [pow(1753, _bitrev8(k), Q) for k in range(N)]
_N_INV_256 = modinv(N, Q)

def ntt256(a):
    """In-place NTT for ML-DSA (length 256). Matches the ML-DSA spec butterfly."""
    a = [x % Q for x in a]
    k = 0
    length = N >> 1
    while length >= 1:
        start = 0
        while start < N:
            k += 1
            zeta = _ZETAS_256[k]
            for j in range(start, start + length):
                t = zeta * a[j + length] % Q
                a[j + length] = (a[j] - t) % Q
                a[j] = (a[j] + t) % Q
            start += length << 1
        length >>= 1
    return a

def intt256(a):
    """In-place INTT for ML-DSA (length 256). Uses negated zetas in reverse."""
    a = [x % Q for x in a]
    k = N
    length = 1
    while length < N:
        start = 0
        while start < N:
            k -= 1
            zeta_neg = (Q - _ZETAS_256[k]) % Q
            for j in range(start, start + length):
                t = a[j]
                a[j] = (t + a[j + length]) % Q
                a[j + length] = zeta_neg * (t - a[j + length]) % Q  # t - b, not b - t
            start += length << 1
        length <<= 1
    for j in range(N):
        a[j] = a[j] * _N_INV_256 % Q
    return a

def poly_inv_ntt(f, n=N, q=Q):
    """
    Invert f in R_q = Z_q[X]/(X^n+1) using NTT (n=256 only).
    Returns None if f is not invertible (some NTT eval is 0 mod q).
    """
    f_ntt = ntt256(list(f))
    if any(x == 0 for x in f_ntt):
        return None  # not invertible: some root of X^n+1 is also a root of f
    inv_ntt = [pow(x, q - 2, q) for x in f_ntt]  # Fermat: x^{-1} = x^{q-2}
    return intt256(inv_ntt)

# ─── Minimal signing oracle ───────────────────────────────────────────────────
# Models only the nonce reuse scenario.
# "Buggy" oracle: rho' = H(key) — message-independent.
# "Correct" oracle: rho' = H(key || mu) — message-dependent.

import numpy as np

_RNG = np.random.default_rng(0)

def shake128_sample(seed: bytes, n: int, q: int, bound: int) -> list:
    """Sample a uniform polynomial from seed via SHAKE-128."""
    h = hashlib.shake_128(seed).digest(n * 3)
    coeffs = []
    i = 0
    while len(coeffs) < n:
        v = struct.unpack_from('<H', h, i % len(h))[0] % (2 * bound + 1) - bound
        coeffs.append(v % q)
        i += 2
    return coeffs[:n]

def sample_secret_poly(n: int, eta: int, q: int, seed: bytes) -> list:
    """Sample secret poly with coeffs in {-eta,...,eta}."""
    h = hashlib.shake_256(seed).digest(n * 2)
    coeffs = []
    for i in range(n):
        v = h[i] % (2 * eta + 1) - eta
        coeffs.append(v % q)
    return coeffs

def sample_challenge(n: int, tau: int, q: int, seed: bytes) -> list:
    """Sample sparse ternary challenge polynomial with exactly tau non-zero ±1 entries."""
    import random
    rng = random.Random(seed)
    c = [0] * n
    positions = rng.sample(range(n), tau)
    for p in positions:
        c[p] = rng.choice([1, q - 1])  # ±1 mod q
    return c

def expand_mask(rho_prime: bytes, nonce: int, gamma1: int, n: int, q: int) -> list:
    """Sample y coefficients uniform in [-gamma1+1, gamma1]."""
    seed = rho_prime + struct.pack('<H', nonce)
    h = hashlib.shake_256(seed).digest(n * 4 + 8)
    coeffs = []
    for i in range(n):
        raw = struct.unpack_from('<I', h, i * 4)[0]
        v = raw % (2 * gamma1) - gamma1 + 1
        coeffs.append(v % q)
    return coeffs

def buggy_sign(key: bytes, message: bytes, s1: list, n: int, q: int, gamma1: int, beta: int, tau: int) -> tuple:
    """
    VULNERABLE signing oracle: rho' = H(key) — OMITS message hash from rho'.
    Same y is used for ALL messages with the same key.
    Returns (z, c) where z = y + c*s1 in R_q.
    """
    rho_prime = hashlib.shake_256(key).digest(64)  # BUG: no mu

    mu = hashlib.shake_256(message).digest(32)

    y = expand_mask(rho_prime, 0, gamma1, n, q)

    # c = H(mu || w1) — we simplify w1 as HighBits(y) for this demo
    c = sample_challenge(n, tau, q, mu + bytes([y[0] % 256]))

    # z = y + c*s1 mod (X^n+1, q)
    cs1 = poly_mul_neg(c, s1, n, q)
    z = poly_add(y, cs1, q)

    # Check norm (simplified — skip hint check for clarity)
    max_z = max(abs(x if x <= q//2 else x - q) for x in z)
    if max_z >= gamma1 - beta:
        # Would retry in real impl; for demo, return anyway (overflow is detectable)
        pass

    return z, c, y  # y exposed ONLY for verification — attacker does NOT see y

# ─── TOY ATTACK (n=4, q=241) ─────────────────────────────────────────────────

def run_toy_attack():
    """
    Full key recovery from 2 signatures on toy ML-DSA.
    n=4, q=241, eta=2. Runs instantly (~ms).
    """
    print("═" * 62)
    print("TOY ATTACK: n=4, q=241")
    print("═" * 62)

    n, q, eta, tau, gamma1, beta = N_TOY, Q_TOY, ETA_TOY, TAU_TOY, GAMMA1_TOY, BETA_TOY

    # Keygen
    key = b"secret-key-toy"
    s1 = [1, -1, 2, 0]   # secret key: small coeffs in {-eta,...,eta}
    s1_q = [x % q for x in s1]
    print(f"\nSecret s1 = {s1}  (target)")

    # Sign two different messages with BUGGY oracle
    m1, m2 = b"message one", b"message two"
    z1, c1, y1 = buggy_sign(key, m1, s1_q, n, q, gamma1, beta, tau)
    z2, c2, y2 = buggy_sign(key, m2, s1_q, n, q, gamma1, beta, tau)

    print(f"\nSig 1 (m='{m1.decode()}'): z1 = {z1}")
    print(f"           challenge c1 = {c1}")
    print(f"Sig 2 (m='{m2.decode()}'): z2 = {z2}")
    print(f"           challenge c2 = {c2}")

    # Verify y reuse (attacker does NOT see y — this is just for verification)
    assert y1 == y2, "y must be equal (nonce reuse confirmed)"
    print(f"\n✓ y reuse confirmed (hidden from attacker): y = {y1}")

    # ATTACK: z1 - z2 = (c1 - c2) * s1  →  s1 = (z1-z2) * inv(c1-c2)
    diff_z = poly_sub(z1, z2, q)
    diff_c = poly_sub(c1, c2, q)

    print(f"\nz1 - z2 = {diff_z}")
    print(f"c1 - c2 = {diff_c}")

    c_inv = poly_inv_neg(diff_c, n, q)
    if c_inv is None:
        print("ERROR: (c1-c2) not invertible in R_q — need different messages")
        return False

    s1_recovered = poly_mul_neg(diff_z, c_inv, n, q)
    # Center to recover signed representation
    s1_rec_signed = [x if x <= q//2 else x - q for x in s1_recovered]

    print(f"\n→ s1 recovered = {s1_rec_signed}")
    print(f"→ s1 actual    = {s1}")

    match = (s1_rec_signed == s1)
    print(f"\n{'✓ FULL KEY RECOVERY — ATTACK SUCCEEDED' if match else '✗ RECOVERY FAILED'}")

    if match:
        print(f"  s1 extracted in 2 signatures from buggy oracle")
        print(f"  No brute force. No lattice reduction. Pure algebra.")

    return match

# ─── FULL ML-DSA-44 ATTACK (n=256, q=8380417) ────────────────────────────────

def run_full_attack(verbose: bool = True):
    """
    Full key recovery from 2 signatures on ML-DSA-44 parameters.
    n=256, q=8380417. Uses NTT for O(n log n) inversion.
    """
    print("\n" + "═" * 62)
    print("FULL ATTACK: n=256, q=8380417 (ML-DSA-44 parameters)")
    print("═" * 62)

    n, q, eta, tau, gamma1, beta = N, Q, ETA, TAU, GAMMA1, BETA

    # Keygen: random secret key with small coefficients
    key = secrets.token_bytes(32)
    s1 = sample_secret_poly(n, eta, q, key + b"s1")
    s1_signed = [x if x <= q//2 else x - q for x in s1]
    print(f"\nKey generated. ||s1||_inf = {max(abs(x) for x in s1_signed)}")
    print(f"  Expected: ||s1||_inf ≤ {eta}  {'✓' if max(abs(x) for x in s1_signed) <= eta else '✗'}")

    # Two messages
    m1 = b"classified document alpha 2026-05-27"
    m2 = b"classified document beta  2026-05-27"

    t_start = time.perf_counter()

    # Sign with BUGGY oracle
    z1, c1, y1 = buggy_sign(key, m1, s1, n, q, gamma1, beta, tau)
    z2, c2, y2 = buggy_sign(key, m2, s1, n, q, gamma1, beta, tau)

    t_sign = time.perf_counter() - t_start

    assert y1 == y2, "y must be equal"
    print(f"\nTwo signatures collected in {t_sign*1000:.1f}ms")
    print(f"  Nonce reuse confirmed: y[0]={y1[0]}, y[1]={y1[1]}")

    # ATTACK
    t_atk = time.perf_counter()

    diff_z = poly_sub(z1, z2, q)
    diff_c = poly_sub(c1, c2, q)

    # NTT-based inversion: O(n log n)
    c_inv = poly_inv_ntt(diff_c, n, q)
    if c_inv is None:
        print("ERROR: (c1-c2) not invertible — retry with different messages")
        return False

    s1_recovered_ntt = poly_mul_neg(diff_z, c_inv, n, q)

    t_atk_end = time.perf_counter() - t_atk

    # Verify recovery
    s1_rec_signed = [x if x <= q//2 else x - q for x in s1_recovered_ntt]
    s1_orig_signed = [x if x <= q//2 else x - q for x in s1]

    match = (s1_rec_signed == s1_orig_signed)

    if verbose:
        print(f"\nAttack computation: {t_atk_end*1000:.1f}ms")
        print(f"Recovered s1[:8]:  {s1_rec_signed[:8]}")
        print(f"Actual    s1[:8]:  {s1_orig_signed[:8]}")
        print(f"||s1_recovered||_inf = {max(abs(x) for x in s1_rec_signed)}")

    print(f"\n{'✓ FULL KEY RECOVERY — ATTACK SUCCEEDED' if match else '✗ RECOVERY FAILED'}")
    if match:
        print(f"  2 signatures → full 256-coefficient secret key s1")
        print(f"  NTT inversion: O(n log n) = O({n} × {n.bit_length()})")
        print(f"  Total wall time: {(t_sign + t_atk_end)*1000:.1f}ms")
        print()
        print(f"  IMPACT: Any ML-DSA implementation where rho' is")
        print(f"  derived without per-message randomness is fully broken")
        print(f"  by an attacker who can request 2 signatures.")

    return match

# ─── Invertibility analysis ───────────────────────────────────────────────────

def analyze_invertibility(trials: int = 1000):
    """
    How often is (c1-c2) invertible in R_q?
    c1, c2 are sparse ternary polys with tau non-zero entries.
    (c1-c2) has coefficients in {-2,-1,0,1,2}.
    """
    print("\n" + "═" * 62)
    print("INVERTIBILITY ANALYSIS: P((c1-c2) invertible in R_q)")
    print("═" * 62)
    print(f"  Testing {trials} random pairs (c1, c2)...")

    invertible = 0
    for i in range(trials):
        seed1 = i.to_bytes(4, 'little')
        seed2 = (i + 10000).to_bytes(4, 'little')
        c1 = sample_challenge(N_TOY, TAU_TOY, Q_TOY, seed1)
        c2 = sample_challenge(N_TOY, TAU_TOY, Q_TOY, seed2)
        diff = poly_sub(c1, c2, Q_TOY)
        if any(x != 0 for x in diff):  # c1 ≠ c2
            inv = poly_inv_neg(diff, N_TOY, Q_TOY)
            if inv is not None:
                invertible += 1

    p = invertible / trials
    print(f"  P(invertible) = {p:.3f}  ({invertible}/{trials})")
    print(f"  Expected: ≈ 1 - 1/q ≈ {1 - 1/Q_TOY:.3f}")
    print(f"  Each pair of distinct messages → {p*100:.0f}% chance of instant key recovery")

# ─── Main ─────────────────────────────────────────────────────────────────────

def self_test():
    """Verify EEA and NTT inversion are correct before running attacks."""
    # Test EEA: f * inv(f) ≡ 1 (mod X^4+1, mod 241)
    n, q = N_TOY, Q_TOY
    f = [3, 1, 2, 0]  # degree-2 poly
    f_q = [x % q for x in f]
    inv_f = poly_inv_neg(f_q, n, q)
    assert inv_f is not None, "EEA: f should be invertible"
    product = poly_mul_neg(f_q, inv_f, n, q)
    assert product[0] == 1 and all(x == 0 for x in product[1:]), \
        f"EEA self-test failed: f * inv(f) = {product}, expected [1,0,0,0]"

    # Test NTT: INTT(NTT(f)) = f
    f256 = list(range(1, N + 1))  # simple test polynomial for n=256
    f256_ntt = ntt256(f256[:])
    f256_rec = intt256(f256_ntt[:])
    assert all(a % Q == b % Q for a, b in zip(f256, f256_rec)), \
        "NTT self-test failed: INTT(NTT(f)) ≠ f"

    # Test NTT inversion: f * inv(f) ≡ 1 in R_q (n=256)
    f256_test = [0] * N
    f256_test[0] = 3; f256_test[1] = 1; f256_test[2] = 2  # small poly
    inv_f256 = poly_inv_ntt(f256_test)
    assert inv_f256 is not None, "NTT inv: test poly should be invertible"
    product256 = poly_mul_neg(f256_test, inv_f256, N, Q)
    assert product256[0] == 1 and all(x % Q == 0 for x in product256[1:]), \
        f"NTT inv self-test failed: f * inv(f)[0] = {product256[0]}"

    print("Self-tests passed: EEA ✓  NTT round-trip ✓  NTT inversion ✓")

if __name__ == "__main__":
    print("PRISM-DSA Phase 5A — Nonce Reuse Key Recovery Attack")
    print("Target: ML-DSA implementations with faulty rho' derivation")
    print()

    self_test()
    print()

    ok1 = run_toy_attack()
    ok2 = run_full_attack(verbose=True)
    analyze_invertibility(trials=200)

    print("\n" + "═" * 62)
    print("SUMMARY")
    print("═" * 62)
    print(f"  Toy  (n=4,   q=241):        {'PASS ✓' if ok1 else 'FAIL ✗'}")
    print(f"  Full (n=256, q=8380417):    {'PASS ✓' if ok2 else 'FAIL ✗'}")
    print()
    print("  Root cause:  rho' = H(key) without mu → y reuse across messages")
    print("  Math:        z1 - z2 = (c1-c2)·s1 → NTT inversion → s1")
    print("  Complexity:  O(n log n) NTT, 2 oracle queries, milliseconds")
    print("  Fix:         rho' = H(key || rnd || mu) [FIPS 204 §5.2, correct]")
    print("  PRISM-DSA:   NOT vulnerable (FIS + per-call rho' derivation)")
