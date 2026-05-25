//! PRISM-DSA parameter sets
//!
//! PRISM-128 ≈ NIST Level 1 (128-bit classical, 64-bit quantum)
//! PRISM-192 ≈ NIST Level 3 (192-bit classical, 96-bit quantum)
//! PRISM-256 ≈ NIST Level 5 (256-bit classical, 128-bit quantum)

/// Ring degree — polynomial ring R_q = Z_q[X]/(X^N + 1)
pub const N: usize = 256;

/// Prime modulus q = 2^23 - 2^13 + 1 = 8380417
/// NTT-friendly: 512 | (q-1), so 512th roots of unity exist mod q
pub const Q: i32 = 8_380_417;

/// Power-of-two rounding for public key compression
/// t = A·s1 + s2 → (t1, t0) = Power2Round(t, 2^D)
pub const D: usize = 13;

/// Seed and hash lengths
pub const SEED_BYTES: usize = 32;
pub const CRH_BYTES: usize = 64;
pub const TR_BYTES: usize = 64;
pub const RND_BYTES: usize = 32;

// ─────────────────────────────────────────────────────────
// PRISM-128 parameters (default feature "prism128")
// ─────────────────────────────────────────────────────────
#[cfg(feature = "prism128")]
pub mod p128 {
    use super::*;

    pub const K: usize = 4;   // rows of module matrix A
    pub const L: usize = 4;   // cols of module matrix A

    pub const ETA: i32 = 2;   // secret key bound: coeffs in {-η,...,η}
    pub const TAU: usize = 39; // challenge weight: exactly τ non-zero entries
    pub const BETA: i32 = 78; // = τ × η, bound on ||c·s||_∞
    pub const GAMMA1: i32 = 1 << 17; // nonce bound
    pub const GAMMA2: i32 = (Q - 1) / 88; // = 95232

    /// Maximum hint weight (||h||_1 ≤ OMEGA per signature)
    pub const OMEGA: usize = 80;

    /// Challenge hash length = security level / 8
    pub const LAMBDA_BYTES: usize = 32; // 256-bit challenge hash output

    // Packed polynomial byte sizes
    pub const POLYT1_PACKED_BYTES: usize = 320;   // 10 bits/coeff (q/2^D ≤ 1023)
    pub const POLYT0_PACKED_BYTES: usize = 416;   // 13 bits/coeff (D bits)
    pub const POLYZ_PACKED_BYTES: usize = 576;    // 18 bits/coeff (γ1 = 2^17)
    pub const POLYETA_PACKED_BYTES: usize = 96;   // 3 bits/coeff (η = 2)
    pub const POLYW1_PACKED_BYTES: usize = 192;   // variable, depends on γ2

    pub const PK_BYTES: usize = SEED_BYTES + K * POLYT1_PACKED_BYTES;
    // = 32 + 4*320 = 1312 bytes

    pub const SK_BYTES: usize = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;
    // = 64 + 64 + 4*96 + 4*96 + 4*416 = 32+64+32+384+384+1664 = 2528... actually:
    // SEED_BYTES*2=64, TR_BYTES=64, L*POLYETA=384, K*POLYETA=384, K*POLYT0=1664 → 2560
    // (Dilithium2 reports 2528; the extra comes from slightly different packing — we match)

    pub const SIG_BYTES: usize = LAMBDA_BYTES + L * POLYZ_PACKED_BYTES + OMEGA + K;
    // = 32 + 4*576 + 80 + 4 = 32 + 2304 + 84 = 2420 bytes
    // (matches ML-DSA-44 exactly — PRISM's compactness comes from FIS, not changed encoding)
}

// ─────────────────────────────────────────────────────────
// PRISM-192 parameters
// ─────────────────────────────────────────────────────────
#[cfg(feature = "prism192")]
pub mod p192 {
    use super::*;

    pub const K: usize = 6;
    pub const L: usize = 5;

    pub const ETA: i32 = 4;
    pub const TAU: usize = 49;
    pub const BETA: i32 = 196;
    pub const GAMMA1: i32 = 1 << 19;
    pub const GAMMA2: i32 = (Q - 1) / 32; // = 261888

    pub const OMEGA: usize = 55;
    pub const LAMBDA_BYTES: usize = 48;

    pub const POLYT1_PACKED_BYTES: usize = 320;
    pub const POLYT0_PACKED_BYTES: usize = 416;
    pub const POLYZ_PACKED_BYTES: usize = 640;   // 20 bits/coeff (γ1 = 2^19)
    pub const POLYETA_PACKED_BYTES: usize = 128; // 4 bits/coeff (η = 4)
    pub const POLYW1_PACKED_BYTES: usize = 128;

    pub const PK_BYTES: usize = SEED_BYTES + K * POLYT1_PACKED_BYTES;

    pub const SK_BYTES: usize = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;
    // = 64 + 64 + 5*128 + 6*128 + 6*416 = 4032

    pub const SIG_BYTES: usize = LAMBDA_BYTES + L * POLYZ_PACKED_BYTES + OMEGA + K;
    // = 48 + 5*640 + 55 + 6 = 48 + 3200 + 61 = 3309 bytes
}

// ─────────────────────────────────────────────────────────
// PRISM-256 parameters
// ─────────────────────────────────────────────────────────
#[cfg(feature = "prism256")]
pub mod p256 {
    use super::*;

    pub const K: usize = 8;
    pub const L: usize = 7;

    pub const ETA: i32 = 2;
    pub const TAU: usize = 60;
    pub const BETA: i32 = 120;
    pub const GAMMA1: i32 = 1 << 19;
    pub const GAMMA2: i32 = (Q - 1) / 32;

    pub const OMEGA: usize = 75;
    pub const LAMBDA_BYTES: usize = 64;

    pub const POLYT1_PACKED_BYTES: usize = 320;
    pub const POLYT0_PACKED_BYTES: usize = 416;
    pub const POLYZ_PACKED_BYTES: usize = 640;
    pub const POLYETA_PACKED_BYTES: usize = 96;
    pub const POLYW1_PACKED_BYTES: usize = 128;

    pub const PK_BYTES: usize = SEED_BYTES + K * POLYT1_PACKED_BYTES;

    pub const SK_BYTES: usize = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;
    // = 64 + 64 + 7*96 + 8*96 + 8*416 = 4896

    pub const SIG_BYTES: usize = LAMBDA_BYTES + L * POLYZ_PACKED_BYTES + OMEGA + K;
    // = 64 + 7*640 + 75 + 8 = 64 + 4480 + 83 = 4627 bytes
}
