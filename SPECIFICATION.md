# PRISM-DSA: Formal Specification

**Post-quantum Ring-based Ideal Signature Mechanism with Fixed-Iteration Signing**

Version 1.0 — Draft for Review  
Status: Research prototype. Not for production use.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Mathematical Preliminaries](#2-mathematical-preliminaries)
3. [Parameter Sets](#3-parameter-sets)
4. [Key Generation](#4-key-generation)
5. [Signing — Fixed-Iteration Signing (FIS)](#5-signing--fixed-iteration-signing-fis)
6. [Verification](#6-verification)
7. [Encoding and Packing](#7-encoding-and-packing)
8. [Security Levels](#8-security-levels)
9. [Differences from ML-DSA (FIPS 204)](#9-differences-from-ml-dsa-fips-204)
10. [Auxiliary Algorithms](#10-auxiliary-algorithms)

---

## 1. Overview

### 1.1 Purpose

PRISM-DSA is a post-quantum digital signature scheme based on the hardness of the
Module Short Integer Solution (Module-SIS) and Module Learning With Errors (Module-LWE)
problems over the cyclotomic ring R_q = Z_q[X]/(X^N + 1). It is designed to provide
the same security guarantees as ML-DSA (FIPS 204) while additionally providing a
**deterministic, fixed-time signing path** that eliminates the variable-iteration
rejection-sampling loop present in the original Fiat-Shamir with Aborts (FSwA) paradigm.

### 1.2 Relationship to ML-DSA (FIPS 204)

PRISM-DSA is a strict superset of ML-DSA in the following sense:

- **Same ring**: R_q = Z_q[X]/(X^256 + 1), q = 8380417
- **Same parameter structure**: (K, L, η, τ, β, γ₁, γ₂, ω) per security level
- **Same key generation**: KeyGen is identical
- **Same verification**: Verify is identical; PRISM-DSA signatures are verified with the
  ML-DSA verifier (for matching parameter sets)
- **Different signing**: ML-DSA uses an unbounded rejection loop; PRISM-DSA uses the
  Fixed-Iteration Signing (FIS) construction (Section 5)
- **Same signature format**: (c̃, z, h) structure and bit-packing are unchanged
- **Same sizes**: PK, SK, and signature byte lengths are identical to ML-DSA-44/65/87

### 1.3 The Fixed-Iteration Signing (FIS) Construction

Standard FSwA signatures loop until a valid signature is found. This loop:
1. Leaks timing information correlated with the secret key and message
2. Creates potential DoS vectors (adversaries can probe for worst-case inputs)
3. Complicates formal verification (unbounded loops require specialized proof techniques)

FIS replaces the unbounded loop with a **fixed number of attempts** (FIS_SLOTS = 64).
In each slot, a full signing attempt is computed unconditionally. The first valid slot is
selected using constant-time conditional moves (cmov) from the `subtle` crate. The
signing time is exactly `FIS_SLOTS × (one full attempt)` regardless of the key, message,
or randomness.

Failure probability (all slots rejected) is ≤ 2^{-27} for all parameter sets.

### 1.4 Notation

| Symbol | Meaning |
|--------|---------|
| R_q | Polynomial ring Z_q[X]/(X^N + 1) |
| N | Ring degree = 256 |
| q | Prime modulus = 8380417 |
| ζ | Primitive 512th root of unity mod q: ζ = 1753 |
| ‖·‖_∞ | ℓ_∞ norm (maximum absolute coefficient) |
| ‖·‖_1 | ℓ_1 norm (sum of absolute values) |
| ⊙ | Pointwise (NTT-domain) multiplication |
| H | SHAKE-256 hash function |
| H₅₁₂ | SHAKE-256 with 512-bit output |
| S_η | Set of polynomials with coefficients in {-η,...,η} |
| B_γ | Set of polynomials with ‖·‖_∞ < γ |

---

## 2. Mathematical Preliminaries

### 2.1 The Cyclotomic Ring

All arithmetic takes place in the polynomial ring

    R = Z[X]/(X^N + 1),    N = 256

with coefficients reduced modulo the prime q = 8380417. The quotient ring is

    R_q = Z_q[X]/(X^N + 1)

The polynomial X^256 + 1 is the 512th cyclotomic polynomial, which is irreducible over Z
and splits into N degree-1 factors over Z_q (since 512 | (q-1), verified below). This
enables the Number Theoretic Transform.

**Verification**: q - 1 = 8380416 = 2^9 × 3 × 5 × 7 × 157. Since 512 = 2^9 divides
q-1, there exist primitive 512th roots of unity mod q. The value ζ = 1753 is a primitive
512th root of unity mod q (precomputed and verified in the reference implementation).

### 2.2 Number Theoretic Transform (NTT)

The NTT is a 256-point negacyclic transform over Z_q. Negacyclicity arises from the
X^N + 1 quotient: the ring product in R_q corresponds to a "negacyclic convolution"
in the coefficient domain.

**Forward NTT**: Given f ∈ R_q with coefficients (f₀, ..., f_{N-1}), the NTT outputs
ĝ = (f̂₀, ..., f̂_{N-1}) where

    f̂_k = Σ_{j=0}^{N-1} f_j · ζ^{(2·BitRev(k)+1)·j}    (mod q)

where BitRev(k) is the bit-reversal of k in log₂(N) = 8 bits.

In the Cooley-Tukey butterfly form (as implemented), the forward NTT operates in-place
with bit-reversed output. The twiddle factors (ZETAS array) are precomputed as

    ZETAS[k] = ζ^{BitRev8(k)} · MONT   (mod q),    k = 0, ..., N-1

where MONT = 2^32 mod q = 4193792 is the Montgomery constant.

**Inverse NTT**: The inverse uses the Gentleman-Sande butterfly. After `invntt_tomont`,
coefficients are in Montgomery domain:

    a[i] ← NTT^{-1}(a)[i] · MONT   (mod q)

The normalization factor is n^{-1} mod q = 8347681 (since 256 · 8347681 ≡ 1 mod q).
Implemented as the Montgomery constant F = n^{-1} · MONT mod q = 41978.

**NTT multiplication**: For f, g ∈ R_q, the ring product h = f · g mod (X^N + 1) is
computed as

    h = NTT^{-1}( NTT(f) ⊙ NTT(g) )

where ⊙ denotes pointwise multiplication. This uses O(N log N) operations instead of
O(N^2) for schoolbook multiplication.

**Montgomery reduction**: For a ∈ [-q·2^31, q·2^31]:

    MontRed(a) = (a - ((a · q^{-1} mod 2^32) · q)) >> 32

where q^{-1} mod 2^32 = 58728449. Output satisfies |MontRed(a)| < q.

**Reduce32**: Maps a ∈ Z to [-6283008, 6283008]:

    Reduce32(a) = a - round(a / q) · q = a - ((a + 2^22) >> 23) · q

**caddq**: Maps a ∈ [-q, q) to [0, q-1):

    caddq(a) = a + ((a >> 31) & q)

### 2.3 Power2Round

For a ∈ Z_q, the Power2Round decomposition produces (a₁, a₀) such that

    a = a₁ · 2^D + a₀

with a₀ ∈ (-2^{D-1}, 2^{D-1}] and D = 13.

**Algorithm** (coefficient-wise):

    a₁ ← ⌈a / 2^D⌉ = (a + 2^{D-1} - 1) >> D
    a₀ ← a - a₁ · 2^D

Note: a₁ ∈ [0, ⌈q / 2^D⌉] = [0, 1023] (10 bits), a₀ ∈ (-4096, 4096] (13 bits signed).

### 2.4 Decompose

For a ∈ Z_q and rounding parameter α = 2γ₂, the Decompose function produces (a₁, a₀) such that

    a ≡ a₁ · α + a₀   (mod q)

with a₀ ∈ (-γ₂, γ₂]. Special case: if a₁ = (q-1)/α, then a₁ ← 0, a₀ ← a₀ - 1.

**Algorithm** (coefficient-wise):

    a ← a mod q   (map to [0, q-1])
    a₀ ← a mod α
    if a₀ > γ₂: a₀ ← a₀ - α   (center in (-γ₂, γ₂])
    a₁ ← (a - a₀) / α
    if a₁ = (q-1)/α: a₁ ← 0; a₀ ← a₀ - 1   (boundary case)

**Implementation note**: The boundary check `a₁ == (q-1)/α` is performed using the
constant-time identity `((a₁ - top) | (top - a₁)) >> 31` where `top = (q-1)/α`.

HighBits(a, α) = a₁ from Decompose(a, α)
LowBits(a, α) = a₀ from Decompose(a, α)

### 2.5 MakeHint and UseHint

These functions enable the verifier to recover the high bits of A·y - c·s₂ + c·t₀
using only the public key (which stores only t₁ = HighBits(A·s₁+s₂, 2^D)).

**MakeHint** (coefficient-wise, a₀ = low bits, a₁ = high bits of w):

    MakeHint(a₀, a₁, γ₂) = 1   if a₀ > γ₂ ∨ a₀ < -γ₂ ∨ (a₀ = -γ₂ ∧ a₁ ≠ 0)
                             0   otherwise

Interpretation: the hint bit records whether the low bits of (w₀ - c·s₂ + c·t₀)
overflow into the high bits.

**UseHint** (coefficient-wise, h ∈ {0,1}, a ∈ [0, q-1]):

    Let (a₁, a₀) = Decompose(a, γ₂), m = (q-1)/(2γ₂)
    UseHint(0, a, γ₂) = a₁
    UseHint(1, a, γ₂) = if a₀ > 0:
                             if a₁ = m-1: 0   else: a₁ + 1
                         else:
                             if a₁ = 0: m-1   else: a₁ - 1

Correctness: UseHint(h, w'_approx, γ₂) recovers w₁ = HighBits(w, 2γ₂) from the
approximation w'_approx = A·z - c·t₁·2^D when the hint was computed correctly during signing.

---

## 3. Parameter Sets

Three parameter sets are defined, targeting NIST security levels 1, 3, and 5.

### 3.1 Shared Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| N | 256 | Ring degree; polynomial ring R_q = Z_q[X]/(X^{256}+1) |
| q | 8380417 | NTT-friendly prime; q = 2^23 - 2^13 + 1 |
| D | 13 | Power-of-two rounding bits for public key compression |
| SEED_BYTES | 32 | Length of seed ξ and derived seeds ρ, key |
| CRH_BYTES | 64 | Length of collision-resistant hash outputs (tr, μ, ρ') |
| TR_BYTES | 64 | Length of public key hash tr |
| RND_BYTES | 32 | Length of per-signing randomness rnd |
| FIS_SLOTS | 64 | Number of fixed signing iterations |

### 3.2 PRISM-128 Parameters (≈ NIST Level 1)

| Parameter | Value | Description |
|-----------|-------|-------------|
| K | 4 | Number of rows in module matrix A |
| L | 4 | Number of columns in module matrix A |
| η | 2 | Secret key bound: s₁, s₂ ∈ S_η^L × S_η^K |
| τ | 39 | Challenge weight: ‖c‖_0 = τ, c_i ∈ {-1, 0, 1} |
| β | 78 | = τ · η; bound on ‖c·s₁‖_∞ and ‖c·s₂‖_∞ |
| γ₁ | 131072 = 2^17 | Nonce bound for mask vector y |
| γ₂ | 95232 = (q-1)/88 | Rounding parameter |
| ω | 80 | Maximum hint weight (‖h‖_1 ≤ ω) |
| λ | 128 bits | Challenge hash output length; LAMBDA_BYTES = 32 |

**Derived byte sizes** (PRISM-128):

| Component | Formula | Bytes |
|-----------|---------|-------|
| PK | 32 + K·320 | 1312 |
| SK | 64 + 64 + L·96 + K·96 + K·416 | 2560 |
| Signature | 32 + L·576 + ω + K | 2420 |

### 3.3 PRISM-192 Parameters (≈ NIST Level 3)

| Parameter | Value | Description |
|-----------|-------|-------------|
| K | 6 | Rows of A |
| L | 5 | Columns of A |
| η | 4 | Secret key bound |
| τ | 49 | Challenge weight |
| β | 196 | = τ · η |
| γ₁ | 524288 = 2^19 | Nonce bound |
| γ₂ | 261888 = (q-1)/32 | Rounding parameter |
| ω | 55 | Maximum hint weight |
| λ | 192 bits | LAMBDA_BYTES = 48 |

**Derived byte sizes** (PRISM-192):

| Component | Formula | Bytes |
|-----------|---------|-------|
| PK | 32 + K·320 | 1952 |
| SK | 64 + 64 + L·128 + K·128 + K·416 | 4032 |
| Signature | 48 + L·640 + ω + K | 3309 |

### 3.4 PRISM-256 Parameters (≈ NIST Level 5)

| Parameter | Value | Description |
|-----------|-------|-------------|
| K | 8 | Rows of A |
| L | 7 | Columns of A |
| η | 2 | Secret key bound |
| τ | 60 | Challenge weight |
| β | 120 | = τ · η |
| γ₁ | 524288 = 2^19 | Nonce bound |
| γ₂ | 261888 = (q-1)/32 | Rounding parameter |
| ω | 75 | Maximum hint weight |
| λ | 256 bits | LAMBDA_BYTES = 64 |

**Derived byte sizes** (PRISM-256):

| Component | Formula | Bytes |
|-----------|---------|-------|
| PK | 32 + K·320 | 2592 |
| SK | 64 + 64 + L·96 + K·96 + K·416 | 4896 |
| Signature | 64 + L·640 + ω + K | 4627 |

---

## 4. Key Generation

Key generation is identical to ML-DSA (FIPS 204), Algorithm 1.

### 4.1 Algorithm: KeyGen(ξ)

**Input**: Optional seed ξ ∈ {0,1}^{256} (if absent, sampled uniformly at random)  
**Output**: Public key pk = (ρ, t₁), Secret key sk = (ρ, key, tr, s₁, s₂, t₀)

```
KeyGen(ξ):
  1. If ξ not provided: ξ ←$ {0,1}^256
  2. (ρ ‖ ρ' ‖ key) ← H₅₁₂(ξ ‖ [K, L])         // 32+64+32 = 128 bytes
  3. A ← ExpandA(ρ)                               // A ∈ R_q^{K×L}
  4. s₁ ← ExpandS(ρ', 0, η)                      // s₁ ∈ S_η^L
  5. s₂ ← ExpandS(ρ', L, η)                      // s₂ ∈ S_η^K
  6. t̂ ← NTT(A) ⊙ NTT(s₁) + NTT(s₂)            // in NTT domain
  7. t ← NTT^{-1}(t̂)                             // t = A·s₁ + s₂ ∈ R_q^K
  8. (t₁, t₀) ← Power2Round(t, 2^D)             // coefficient-wise
  9. pk ← Encode(ρ, t₁)                          // 32 + K·320 bytes
  10. tr ← H₅₁₂(pk)                              // 64-byte public key hash
  11. sk ← Encode(ρ, key, tr, s₁, s₂, t₀)
  12. return (pk, sk)
```

**Notes**:
- Step 2 uses SHAKE-256 with output 128 bytes. The domain separator `[K, L]` distinguishes
  key types across parameter sets.
- Step 3: ExpandA generates A[i][j] = SHAKE256(ρ ‖ j ‖ i) for i ∈ [0,K), j ∈ [0,L),
  using rejection sampling to draw uniform elements of Z_q.
- Steps 4–5: ExpandS uses SHAKE256(ρ' ‖ nonce) with coefficient sampling using the
  η=2 (3-bit) or η=4 (4-bit) technique described in Section 10.
- Step 8 is applied coefficient-wise across all K polynomials.
- tr binds the public key to all subsequent signatures via the message digest μ.

### 4.2 ExpandA: Matrix Expansion

```
ExpandA(ρ):
  for i = 0 to K-1:
    for j = 0 to L-1:
      XOF ← SHAKE256(ρ ‖ [j, i])
      A[i][j] ← RejUniform(XOF)   // draw N uniform coefficients in [0, q-1]
  return A
```

**RejUniform**: Draw bytes 3 at a time; interpret as 23-bit value;
accept if < q (acceptance probability = q/2^23 ≈ 0.9998).

### 4.3 ExpandS: Secret Vector Sampling

For η = 2 (PRISM-128, PRISM-256):
```
ExpandS<η=2>(ρ', nonce, N):
  XOF ← SHAKE256(ρ' ‖ nonce_le16)
  for i = 0 to N/2 - 1:
    b ← XOF.read(1)
    b₀ ← b & 0x0F; b₁ ← b >> 4
    if b₀ < 15: t₀ = b₀ - floor(205·b₀/1024)·5; coeffs[2i] = 2 - t₀
    if b₁ < 15: t₁ = b₁ - floor(205·b₁/1024)·5; coeffs[2i+1] = 2 - t₁
```
(The expression `b - (205*b >> 10)*5` computes `b mod 5` using integer arithmetic,
a standard branchless modular reduction for small values.)

For η = 4 (PRISM-192):
```
ExpandS<η=4>(ρ', nonce, N):
  XOF ← SHAKE256(ρ' ‖ nonce_le16)
  for i = 0 to N/2 - 1:
    b ← XOF.read(1)
    b₀ ← b & 0x0F; b₁ ← b >> 4
    if b₀ < 9: coeffs[2i] = 4 - b₀
    if b₁ < 9: coeffs[2i+1] = 4 - b₁
```

---

## 5. Signing — Fixed-Iteration Signing (FIS)

FIS is the core innovation of PRISM-DSA. It replaces the unbounded rejection loop
of ML-DSA with a fixed `FIS_SLOTS = 64` iterations.

### 5.1 Algorithm: FIS-Sign(sk, m, ctx)

**Input**: sk (secret key bytes), m (message), ctx ∈ {0,1}^* with |ctx| ≤ 255  
**Output**: σ = (c̃, z, h) or FAILURE

**Precondition**: |ctx| ≤ 255

```
FIS-Sign(sk, m, ctx):
  1. Parse sk → (ρ, key, tr, s₁, s₂, t₀)
  2. A ← ExpandA(ρ); NTT(A)                     // precompute, constant across slots
  3. ŝ₁ ← NTT(s₁); ŝ₂ ← NTT(s₂); t̂₀ ← NTT(t₀)
  4. μ ← H₅₁₂(tr ‖ [0x00, |ctx|] ‖ ctx ‖ m)    // 64-byte message digest
  5. rnd ←$ {0,1}^256                             // fresh per-signature randomness
  6. ρ' ← H₅₁₂(key ‖ rnd ‖ μ)                   // 64-byte signing randomness
  7. result_bytes ← 0^{SIG_BYTES}                 // output buffer
  8. found ← 0                                    // flag: first valid slot found?

  // FIS loop: exactly FIS_SLOTS iterations, no early exit
  9. for slot = 0 to FIS_SLOTS - 1:
       a. nonce ← slot × L
       b. y ← ExpandMask(ρ', nonce, γ₁)          // y ∈ (-γ₁, γ₁]^L
       c. ŷ ← NTT(y)
       d. w ← NTT^{-1}(A ⊙ ŷ)                   // w = A·y ∈ R_q^K
       e. (w₁, w₀) ← Decompose(w, 2γ₂)           // coefficient-wise
       f. c̃ ← H_{2λ}(μ ‖ Encode(w₁))             // challenge hash
       g. c ← SampleInBall(c̃, τ)                 // sparse ternary challenge
       h. ĉ ← NTT(c)
       i. cs₁ ← NTT^{-1}(ĉ ⊙ ŝ₁)               // cs₁ = c·s₁
       j. z ← y + cs₁                             // response
       k. valid ← (‖z‖_∞ < γ₁ - β)              // CT norm check
       l. cs₂ ← NTT^{-1}(ĉ ⊙ ŝ₂)               // cs₂ = c·s₂
       m. w₀' ← w₀ - cs₂                         // adjusted low bits
       n. valid ← valid & (‖w₀'‖_∞ < γ₂ - β)   // CT norm check
       o. ct₀ ← NTT^{-1}(ĉ ⊙ t̂₀)               // ct₀ = c·t₀
       p. valid ← valid & (‖ct₀‖_∞ < γ₂)        // CT norm check
       q. h ← MakeHint(w₀' + ct₀, w₁, γ₂)       // hint polynomial
       r. valid ← valid & (‖h‖_1 ≤ ω)           // hint weight check
       s. candidate ← Encode(c̃, z, h)
       t. use_this ← CT_AND(valid, CT_NOT(found)) // first valid slot only
       u. result_bytes ← CT_SELECT(use_this, candidate, result_bytes)
       v. found ← found | valid                   // latch: never resets to 0

  10. if found = 0: return FAILURE
  11. return result_bytes
```

### 5.2 Constant-Time Selection

Step 9t–9v implement constant-time conditional selection using the `subtle` crate:

```
use_this = Choice::from(valid & (found ^ 1))
for each byte i:
  result_bytes[i] = u8::conditional_select(result_bytes[i], candidate[i], use_this)
found |= valid
```

**Semantics**: `use_this` is 1 iff this slot is valid AND no previous slot was valid.
The `conditional_select` is a cmov (conditional move) operation with no data-dependent
branching on the `choice` bit.

### 5.3 Norm Checks

All norm checks in the FIS loop use the **constant-time** function `check_norm`:

```
check_norm(v, bound):
  ok ← true
  for i = 0 to N-1:
    r ← Reduce32(v[i])    // maps to approximately [-q/2, q/2]
    ok ← ok & (|r| < bound)
  return ok
```

This iterates all N = 256 coefficients regardless of intermediate results. No early exit.

The vector version iterates all polynomials in the vector:
```
norm_check_ct(vec, bound):
  ok ← 1
  for each polynomial p in vec:
    ok ← ok & check_norm(p, bound)
  return ok
```

### 5.4 Hint Weight Check

After computing the hint polynomial:

```
h_weight ← Σᵢ Σⱼ h[i][j]   // sum all 1-bits across K polynomials
valid ← valid & ((OMEGA.wrapping_sub(h_weight) >> (usize::BITS - 1)) ^ 1)
```

The expression `(OMEGA.wrapping_sub(h_weight) >> (BITS-1)) ^ 1` evaluates to 1 iff
`h_weight ≤ OMEGA`, using wrapping subtraction and sign-bit extraction to avoid a branch.

### 5.5 ExpandMask

For γ₁ = 2^17 (PRISM-128), each coefficient requires 18 bits:
```
ExpandMask<γ₁=2^17>(ρ', nonce):
  for i = 0 to L-1:
    XOF ← SHAKE256(ρ' ‖ (nonce+i)_le16)
    buf ← XOF.read(N × 18 / 8)   // 576 bytes
    for j = 0 to N/4 - 1:
      [z₀,z₁,z₂,z₃] from 9 bytes of buf, each 18 bits
      y[i][4j+k] = γ₁ - zₖ
```

For γ₁ = 2^19 (PRISM-192, PRISM-256), each coefficient requires 20 bits:
```
ExpandMask<γ₁=2^19>(ρ', nonce):
  for i = 0 to L-1:
    XOF ← SHAKE256(ρ' ‖ (nonce+i)_le16)
    buf ← XOF.read(N × 20 / 8)   // 640 bytes
    for j = 0 to N/2 - 1:
      [z₀,z₁] from 5 bytes of buf, each 20 bits
      y[i][2j+k] = γ₁ - zₖ
```

---

## 6. Verification

Verification is identical to ML-DSA (FIPS 204) Algorithm 3. PRISM-DSA signatures
output by FIS are verified by the standard ML-DSA verifier.

### 6.1 Algorithm: Verify(pk, m, ctx, σ)

**Input**: pk (public key bytes), m (message), ctx (context), σ (signature bytes)  
**Output**: Accept or Reject

```
Verify(pk, m, ctx, σ):
  1. Parse pk → (ρ, t₁)
  2. Parse σ → (c̃, z, h)
  3. if |ctx| > 255: return Reject
  4. if ‖z‖_∞ ≥ γ₁ - β: return Reject             // Check 1: bound on z
  5. if ‖h‖_1 > ω: return Reject                   // Check 2: hint weight
  6. tr ← H₅₁₂(pk)
  7. μ ← H₅₁₂(tr ‖ [0x00, |ctx|] ‖ ctx ‖ m)
  8. c ← SampleInBall(c̃, τ)
  9. A ← ExpandA(ρ); NTT(A)
  10. ẑ ← NTT(z); t̂₁ ← NTT(t₁)
  11. Az ← NTT^{-1}(A ⊙ ẑ)                         // A·z
  12. ct₁ ← NTT^{-1}(NTT(c) ⊙ t̂₁)                 // c·t₁
  13. w'_approx[i][j] ← (Az[i][j] - ct₁[i][j] · 2^D) mod q   // in [0, q-1]
  14. w₁' ← UseHint(h, w'_approx, γ₂)              // recover high bits
  15. c̃' ← H_{2λ}(μ ‖ Encode(w₁'))
  16. if CT_EQ(c̃, c̃'): return Accept else: return Reject
```

**Implementation note** (Step 13): The computation `Az[j] - ct₁[j] · 2^D` requires
64-bit arithmetic to avoid overflow (ct₁[j] can be near q ≈ 2^23, so `ct₁[j] · 2^D`
can reach 2^36). Implemented in i64 with `rem_euclid` reduction to [0, q-1].

**Step 16**: Uses `subtle::ConstantTimeEq` for constant-time byte comparison,
preventing timing oracles on the challenge hash comparison.

---

## 7. Encoding and Packing

All packed representations are defined such that `Unpack(Pack(x)) = x` for all valid inputs.

### 7.1 t₁ Packing: 10 bits per coefficient

t₁ coefficients ∈ [0, ⌈q/2^D⌉ - 1] = [0, 1023]. Packed 4 coefficients in 5 bytes:

```
// Pack coefficients [a₀, a₁, a₂, a₃], each 10 bits, into 5 bytes:
out[5i]   =  a₀ & 0xFF
out[5i+1] = (a₀ >> 8) | (a₁ << 2)    // 2 bits of a₀, 6 bits of a₁
out[5i+2] = (a₁ >> 6) | (a₂ << 4)
out[5i+3] = (a₂ >> 4) | (a₃ << 6)
out[5i+4] =  a₃ >> 2
```

Total per polynomial: 256 × 10 / 8 = 320 bytes.

### 7.2 t₀ Packing: 13 bits per coefficient (signed)

t₀ coefficients ∈ (-2^{D-1}, 2^{D-1}] = (-4096, 4096]. Stored as unsigned:
`u = 2^{D-1} - t₀ ∈ [0, 8191]`. Packed 8 coefficients in 13 bytes.

```
// u[k] = 2^12 - t₀[8i+k], for k = 0..7
out[13i]    =  u[0] & 0xFF
out[13i+1]  = (u[0] >> 8) | (u[1] << 5)
out[13i+2]  =  u[1] >> 3
out[13i+3]  = (u[1] >> 11) | (u[2] << 2)
out[13i+4]  = (u[2] >> 6) | (u[3] << 7)
out[13i+5]  =  u[3] >> 1
out[13i+6]  = (u[3] >> 9) | (u[4] << 4)
out[13i+7]  =  u[4] >> 4
out[13i+8]  = (u[4] >> 12) | (u[5] << 1)
out[13i+9]  = (u[5] >> 7) | (u[6] << 6)
out[13i+10] =  u[6] >> 2
out[13i+11] = (u[6] >> 10) | (u[7] << 3)
out[13i+12] =  u[7] >> 5
```

Total: 256 × 13 / 8 = 416 bytes.

### 7.3 η=2 Packing: 3 bits per coefficient

Coefficients ∈ {-2,-1,0,1,2}, stored as `2 - coeff ∈ {0,1,2,3,4}`. Packed 8 in 3 bytes:

```
// a[k] = 2 - coeff[8i+k], for k=0..7, each 3 bits
out[3i]   = a[0] | (a[1] << 3) | (a[2] << 6)
out[3i+1] = (a[2] >> 2) | (a[3] << 1) | (a[4] << 4) | (a[5] << 7)
out[3i+2] = (a[5] >> 1) | (a[6] << 2) | (a[7] << 5)
```

Total: 256 × 3 / 8 = 96 bytes.

### 7.4 η=4 Packing: 4 bits per coefficient

Coefficients ∈ {-4,...,4}, stored as `4 - coeff ∈ {0,...,8}`. Packed 2 in 1 byte:

```
out[i] = (4 - coeff[2i]) | ((4 - coeff[2i+1]) << 4)
```

Total: 256 × 4 / 8 = 128 bytes.

### 7.5 z Packing (γ₁ = 2^17): 18 bits per coefficient

z coefficients ∈ (-γ₁, γ₁], stored as `γ₁ - z ∈ [0, 2γ₁)`. Packed 4 in 9 bytes:

```
// v[k] = γ₁ - z[4i+k], each 18 bits
out[9i]   =  v[0] & 0xFF
out[9i+1] =  v[0] >> 8
out[9i+2] = (v[0] >> 16) | (v[1] << 2)   // 2 bits of v[0], 6 of v[1]
out[9i+3] =  v[1] >> 6
out[9i+4] = (v[1] >> 14) | (v[2] << 4)
out[9i+5] =  v[2] >> 4
out[9i+6] = (v[2] >> 12) | (v[3] << 6)
out[9i+7] =  v[3] >> 2
out[9i+8] =  v[3] >> 10
```

Total: 256 × 18 / 8 = 576 bytes.

### 7.6 z Packing (γ₁ = 2^19): 20 bits per coefficient

Packed 4 in 10 bytes (2 coefficients per 5-byte group):

```
// v[k] = γ₁ - z[4i+k], each 20 bits
out[10i]   =  v[0] & 0xFF
out[10i+1] =  v[0] >> 8
out[10i+2] = (v[0] >> 16) | (v[1] << 4)  // 4 bits of v[0], 4 of v[1]
out[10i+3] =  v[1] >> 4
out[10i+4] =  v[1] >> 12
out[10i+5] =  v[2] & 0xFF
out[10i+6] =  v[2] >> 8
out[10i+7] = (v[2] >> 16) | (v[3] << 4)
out[10i+8] =  v[3] >> 4
out[10i+9] =  v[3] >> 12
```

Total: 256 × 20 / 8 = 640 bytes.

### 7.7 w₁ Packing (γ₂ = (q-1)/88): 6 bits per coefficient

w₁ coefficients ∈ [0, (q-1)/(2γ₂) - 1] = [0, 43]. Packed 4 in 3 bytes:

```
out[3i]   = coeff[4i]       | (coeff[4i+1] << 6)
out[3i+1] = (coeff[4i+1] >> 2) | (coeff[4i+2] << 4)
out[3i+2] = (coeff[4i+2] >> 4) | (coeff[4i+3] << 2)
```

Total: 192 bytes.

### 7.8 w₁ Packing (γ₂ = (q-1)/32): 4 bits per coefficient

w₁ coefficients ∈ [0, (q-1)/(2γ₂) - 1] = [0, 15]. Packed 2 in 1 byte:

```
out[i] = coeff[2i] | (coeff[2i+1] << 4)
```

Total: 128 bytes.

### 7.9 Hint Packing

The hint vector h ∈ {0,1}^{K×N} is sparse (‖h‖_1 ≤ ω). Stored in ω + K bytes:

**Format**: First ω bytes store the positions of all 1-bits, sorted within each
polynomial and then concatenated. The K offset bytes `out[ω+i]` store the running
count of 1-bits through polynomial i (i.e., the exclusive end index for polynomial i's
1-bit positions in the positions array).

```
// Packing:
idx ← 0
for i = 0 to K-1:
  for j = 0 to N-1:
    if h[i][j] = 1 and idx < ω:
      out[idx] ← j
      idx ← idx + 1
  out[ω + i] ← idx   // cumulative count

// Unpacking with validation:
k ← 0
for i = 0 to K-1:
  end ← buf[ω + i]
  if end < k or end > ω: return INVALID
  for j = k to end-1:
    if j > k and buf[j] ≤ buf[j-1]: return INVALID   // must be strictly increasing
    h[i][buf[j]] ← 1
  k ← end
```

Strictly increasing positions guarantee canonical encoding and uniqueness.

### 7.10 Signature Layout

```
σ = c̃ ‖ z[0] ‖ z[1] ‖ ... ‖ z[L-1] ‖ hint_positions ‖ hint_offsets
```

| Field | Offset | Size (PRISM-128) |
|-------|--------|-----------------|
| c̃ | 0 | 32 bytes (λ=128) |
| z[0..L) | 32 | L × 576 bytes |
| hint positions | 32 + 2304 | 80 bytes (ω=80) |
| hint offsets | 2416 | 4 bytes (K=4) |
| **Total** | | **2420 bytes** |

---

## 8. Security Levels

### 8.1 Claimed Security

| Variant | Classical | Quantum | NIST Level |
|---------|-----------|---------|------------|
| PRISM-128 | 128 bits | 64 bits | Level 1 |
| PRISM-192 | 192 bits | 96 bits | Level 3 |
| PRISM-256 | 256 bits | 128 bits | Level 5 |

Security is defined against:
- **EUF-CMA**: Existential unforgeability under chosen message attack
- **Strong unforgeability**: Infeasibility of generating new valid (m, σ) pairs for
  seen m with new σ (ML-DSA provides this under Module-SIS)

### 8.2 Hardness Assumptions

**Module-SIS_{q,K+L,β}**: Given A ∈ R_q^{K×L} uniform at random, find (z₁, z₂) ∈ R_q^{K+L}
with ‖(z₁ ‖ z₂)‖_∞ ≤ β and A·z₁ + z₂ = 0 (equivalently, find a short vector in a
module lattice with the L-structure). This underlies unforgeability.

**Module-LWE_{q,K,η}**: Given (A, A·s + e) where s ∈ S_η^L and e ∈ S_η^K, distinguish
from (A, u) with u uniform. This underlies key secrecy (hiding s₁, s₂ in the public key).

Both assumptions are believed hard for polynomial-time quantum algorithms when the
parameters satisfy the conditions specified in [FIPS 204].

### 8.3 FIS Security Notes

FIS does not weaken ML-DSA's security guarantees:

1. **Signatures are valid ML-DSA signatures**: The FIS acceptance conditions (norm checks,
   hint weight bound) are identical to ML-DSA. Any valid PRISM-DSA signature passes
   ML-DSA verification.

2. **Output distribution is identical**: Each accepted slot produces a signature drawn
   from the same distribution as a single successful ML-DSA signing attempt
   (proof sketch in `docs/security-proof.md`).

3. **Timing uniformity**: Signing time is deterministically FIS_SLOTS × (one iteration),
   independent of key and message. Formal statement and caveat in Section 9.

### 8.4 Failure Probability

| Variant | p_accept | P(failure) = (1-p_accept)^64 |
|---------|----------|-------------------------------|
| PRISM-128 | ≈ 0.22 | ≈ 2^{-27} ≈ 7.5 × 10^{-9} |
| PRISM-192 | ≈ 0.17 | ≈ 2^{-22} ≈ 2.4 × 10^{-7} |
| PRISM-256 | ≈ 0.17 | ≈ 2^{-22} ≈ 2.4 × 10^{-7} |

For applications requiring lower failure probability, FIS_SLOTS can be increased:
- FIS_SLOTS = 100: P(failure, 128) ≈ 2^{-42}
- FIS_SLOTS = 200: P(failure, 128) ≈ 2^{-83}

---

## 9. Differences from ML-DSA (FIPS 204)

### 9.1 What Changes

| Aspect | ML-DSA | PRISM-DSA |
|--------|--------|-----------|
| Signing loop | Unbounded rejection loop | Exactly FIS_SLOTS = 64 iterations |
| Timing | Variable; ~4–8 iterations avg | Fixed; always FIS_SLOTS iterations |
| Abort condition | `continue` on rejection | CT selection; no control-flow abort |
| Output selection | First successful iteration | CT_SELECT(first_valid) |
| Failure mode | Never fails in practice | Fails with probability ≤ 2^{-27} |
| Formal CT proof | Not claimed (unbounded loop) | Timing uniformity claimed (caveat: NTT) |

### 9.2 What Does Not Change

| Aspect | Status |
|--------|--------|
| Ring and modulus | R_q = Z_q[X]/(X^256+1), q = 8380417 |
| Parameter sets | (K, L, η, τ, β, γ₁, γ₂, ω) identical to ML-DSA-44/65/87 |
| KeyGen algorithm | Identical to ML-DSA |
| Verify algorithm | Identical to ML-DSA |
| Signature format | (c̃, z, h) with identical packing |
| Signature sizes | Identical to ML-DSA-44/65/87 |
| Public key sizes | Identical |
| Hash functions | SHAKE-256 everywhere (same as ML-DSA) |
| ExpandA | Identical |
| ExpandS | Identical |
| SampleInBall | Identical |
| Acceptance conditions | Identical norm bounds |

### 9.3 ML-DSA Compatibility

PRISM-DSA signatures ARE valid ML-DSA signatures (for matching K,L parameters),
and are verified by the ML-DSA verifier. The converse is also true: ML-DSA signatures
can be verified by the PRISM-DSA verifier. The only interoperability concern is that
a PRISM-DSA signer will sometimes fail (2^{-27}) where an ML-DSA signer would not.

### 9.4 Open Deviations

The following items deviate from a production-ready implementation:

1. **NTT not fully CT**: `montgomery_reduce` and `reduce32` contain arithmetic that may
   branch on modern processors due to the `>> 23` shift and potential carry operations.
   A hardware-verified constant-time NTT is future work.

2. **No Zeroize on drop**: SecretKey structs should implement `Zeroize` to clear secret
   material from heap on drop. TODO in implementation.

3. **No formal QROM proof**: The output distribution claim (Section 8.3, item 2) is
   argued informally. Formal verification via EasyCrypt is planned.

4. **SampleInBall not CT**: The Fisher-Yates position sampling uses rejection (variable
   number of XOF reads), but this is acceptable since c̃ is a public hash output.

---

## 10. Auxiliary Algorithms

### 10.1 SampleInBall

Samples a sparse ternary polynomial c ∈ R_q with exactly τ nonzero coefficients,
each ±1. The input c̃ is a public challenge hash.

```
SampleInBall(c̃, τ):
  XOF ← SHAKE256(c̃)
  sign_bytes ← XOF.read(8)           // 64 sign bits
  signs ← u64::from_le_bytes(sign_bytes)
  c ← [0; N]
  for i = N-τ to N-1:               // Fisher-Yates over last τ positions
    j ← rejection_sample(XOF, i)   // j ∈ {0, ..., i} (variable-time; input is public)
    c[i] ← c[j]
    c[j] ← 1 - 2·(signs & 1)       // ±1 based on sign bit
    signs >>= 1
  return c
```

### 10.2 Hash Functions

All hash functions use SHAKE-256 (SHA-3 extendable output function):

- `H₅₁₂(inputs...)`: SHAKE-256 with 512-bit (64-byte) output
- `H_{2λ}(inputs...)`: SHAKE-256 with 2λ-bit output (λ = 128, 192, or 256 bits)

Multiple input segments are concatenated before hashing (domain separation via
context bytes where applicable).

### 10.3 Serialized Key Structure

**Public Key** (PRISM-128 example, 1312 bytes):
```
pk = ρ[0..32) ‖ PackT1(t₁[0]) ‖ PackT1(t₁[1]) ‖ ... ‖ PackT1(t₁[K-1])
```

**Secret Key** (PRISM-128 example, 2560 bytes):
```
sk = ρ[0..32) ‖ key[0..32) ‖ tr[0..64)
   ‖ PackEta2(s₁[0]) ‖ ... ‖ PackEta2(s₁[L-1])
   ‖ PackEta2(s₂[0]) ‖ ... ‖ PackEta2(s₂[K-1])
   ‖ PackT0(t₀[0]) ‖ ... ‖ PackT0(t₀[K-1])
```

For PRISM-192 (η=4): PackEta2 is replaced by PackEta4 (128 bytes/polynomial).

---

*End of PRISM-DSA Specification v1.0*
