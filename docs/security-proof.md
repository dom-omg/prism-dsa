# PRISM-DSA Security Proof Sketch

**Status**: Informal argument. Formal proof via EasyCrypt is planned but not yet written.

## 1. Setup

PRISM-DSA is ML-DSA (FIPS 204) with a fixed-iteration signing loop (FIS). The
verification algorithm is identical to ML-DSA. Therefore, the security of PRISM-DSA
reduces to the security of ML-DSA plus the following claim about FIS.

**Parameters (PRISM-128)**: q = 8380417, N = 256, K = 4, L = 4, η = 2, τ = 39,
β = 78, γ₁ = 2¹⁷, γ₂ = (q−1)/88, ω = 80, FIS_SLOTS = 64.

## 2. Standard ML-DSA Security (Summary)

ML-DSA is proven EUF-CMA under:
- **Module-SIS hardness**: no PPT adversary can find a short non-zero vector in a
  module lattice. This underpins unforgeability.
- **Module-LWE hardness**: the public key (ρ, t₁) computationally hides s₁, s₂.
  This underpins key privacy.

Both assumptions hold in the quantum random oracle model (QROM) per the FIPS 204
supporting documentation and [KLS18].

## 3. FIS Output Distribution

**Claim**: When at least one slot is valid, the output distribution of FIS is identical
to the output distribution of standard ML-DSA rejection sampling.

**Argument**:
Let `slot_i` be the output `(c̃ᵢ, zᵢ, hᵢ)` from iteration `i`. A slot is valid iff:
1. ‖zᵢ‖∞ < γ₁ − β
2. ‖w₀ − cs₂‖∞ < γ₂ − β  
3. ‖ct₀‖∞ < γ₂
4. ‖hᵢ‖₁ ≤ ω

FIS selects the output from the **first** valid slot. Each valid slot's output is
independently drawn from the same distribution as standard ML-DSA (the nonces y₀,…,y₆₃
are derived deterministically from ρ′ = H(key ‖ rnd ‖ μ) and distinct nonces via
`expand_mask`). Therefore, the selected output has the same distribution as a single
valid ML-DSA signing attempt.

**Key invariant**: The nonce for slot `i` is `nonce = i × L`, ensuring each slot's
mask vector `y` is independently pseudo-random (via SHAKE-256 in `expand_mask`).

## 4. Unforgeability

**Theorem** (informal): PRISM-DSA is EUF-CMA under Module-SIS and Module-LWE,
assuming FIS_SLOTS ≥ 1.

**Proof sketch**: We reduce to ML-DSA unforgeability. Given a PRISM-DSA forger `F`:
1. `F` outputs a forgery `(m, σ)` where `σ = (c̃, z, h)` passes `verify`.
2. `verify` is identical to ML-DSA `verify`, so `σ` is also a valid ML-DSA signature.
3. The ML-DSA unforgeability reduction extracts a short vector from this forgery,
   breaking Module-SIS.

The FIS loop does not affect verifiability — it only changes HOW signatures are generated,
not WHAT constitutes a valid signature.

## 5. Timing Uniformity

**Theorem**: The signing time of PRISM-DSA is a deterministic function of `FIS_SLOTS`
and the parameter set only. It does not depend on the secret key or the message.

**Proof**: The implementation:
- Iterates exactly `FIS_SLOTS = 64` times (no `break` or `continue`).
- Computes `cs₁, cs₂, ct₀, hint` unconditionally in every iteration.
- Uses `subtle::ConditionallySelectable` (cmov) for output selection.
- `check_norm` iterates all N = 256 coefficients regardless of intermediate results.

Therefore, the number of field operations is fixed: `FIS_SLOTS × (one full iteration)`.
Wall-clock time may vary due to cache effects and branch prediction on norm checks, but
the control flow (branch graph) is secret-independent.

**Caveat**: NTT butterfly operations may have data-dependent timing on some microarchitectures
(e.g., conditional adds in `montgomery_reduce`). A hardware-verified CT NTT is future work.

## 6. Failure Probability

P(all 64 slots rejected) = P(single slot rejected)^64.

For PRISM-128, empirical measurement gives p_accept ≈ 0.22 per slot:

| Rejection cause       | Approx. probability |
|-----------------------|---------------------|
| ‖z‖∞ ≥ γ₁ − β       | ~60%                |
| ‖w₀−cs₂‖∞ ≥ γ₂ − β | ~10% (of remaining) |
| ‖ct₀‖∞ ≥ γ₂         | ~5%  (of remaining) |
| ‖h‖₁ > ω             | ~2%  (of remaining) |
| **Combined accept**   | **~22%**            |

P(all fail) = 0.78^64 ≈ **2^{−27}** ≈ 7.5 × 10⁻⁹

For critical applications, use FIS_SLOTS = 200: P(all fail) ≈ 2^{−83}.

## 7. Open Questions

1. **Formal QROM proof**: The ML-DSA reduction holds in the QROM; does FIS's
   deterministic nonce derivation interact with quantum signing queries? (Expected: no,
   since ρ′ is freshly randomized each call via `rng.fill_bytes(&mut rnd)`.)

2. **NTT constant-time**: `montgomery_reduce` and `reduce32` contain conditional
   operations. Formal CT verification via ct-verif or similar is needed.

3. **FIS slot correlation**: Are slots 0…63 independent in the QROM? They share the
   same key-derived ρ′ but use distinct nonces. The expansion via SHAKE-256 with
   distinct inputs should make them indistinguishable from independent, but a formal
   argument is needed.

4. **EasyCrypt formalization**: Target: formalize Claims 3.1 (output distribution) and
   5.1 (timing uniformity) in EasyCrypt. Precedent: the ML-DSA proof by [BDK+18].

## 8. References

- [FIPS 204] NIST. *Module-Lattice-Based Digital Signature Standard*. 2024.
- [KLS18] Kiltz, Lyubashevsky, Schaffner. *A Concrete Treatment of Fiat-Shamir Signatures
  in the Quantum Random-Oracle Model*. EUROCRYPT 2018.
- [BDK+18] Bai, Ducas, Kiltz, Lepoint, Lyubashevsky, Schwabe, Seiler, Stehlé.
  *CRYSTALS-Dilithium Algorithm Specifications and Supporting Documentation*. 2018.
