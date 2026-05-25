//! # PRISM-DSA
//!
//! **Post-quantum Ring-based Ideal Signature Mechanism with Fixed-Iteration Signing**
//!
//! ## Overview
//!
//! PRISM-DSA is a post-quantum digital signature scheme based on the hardness of
//! Module-SIS and Module-LWE over the cyclotomic ring R_q = Z_q[X]/(X^256 + 1).
//!
//! ## Core Innovation: Fixed-Iteration Signing (FIS)
//!
//! Standard Fiat-Shamir with Aborts (FSwA) signatures use an unbounded rejection loop.
//! This creates:
//! - **Timing side-channels**: the number of iterations correlates with the secret key
//! - **DoS vectors**: adversaries can probe for worst-case signing inputs
//! - **Formal verification difficulty**: unbounded loops require special treatment in proof systems
//!
//! PRISM-DSA solves this with **Fixed-Iteration Signing**:
//! - Always runs exactly `FIS_SLOTS = 64` signing attempts
//! - Uses constant-time conditional selection to choose the first valid slot
//! - Failure probability: `(p_reject)^64 ≈ 2^{-27}` for PRISM-128
//! - **Timing uniformity**: signing time is `FIS_SLOTS × (one attempt)`, always
//!
//! ## Security
//!
//! PRISM-128 targets NIST Level 1 (128-bit classical security, 64-bit quantum security).
//!
//! Hardness assumptions:
//! 1. Module-SIS_{q, k+l, β}: find short vector in a module lattice (unforgeability)
//! 2. Module-LWE_{q, k, η}: distinguish noisy module product from random (key secrecy)
//!
//! Both assumptions are standard NIST PQC assumptions (same as ML-DSA).
//!
//! **NOT PROVEN**: No formal reduction has been written. FIS security is conjectured
//! based on the proof sketch in `sign.rs`. Formal verification via EasyCrypt is planned.
//!
//! ## API
//!
//! ```rust,ignore
//! use prism_dsa::{keygen128, sign128, verify128};
//! use rand::rngs::OsRng;
//!
//! // Key generation
//! let (pk, sk) = keygen128(&mut OsRng).unwrap();
//!
//! // Signing
//! let msg = b"Hello, post-quantum world";
//! let ctx = b"my-application-v1";
//! let sig = sign128(&sk.bytes, msg, ctx, &mut OsRng).unwrap();
//!
//! // Verification
//! verify128(&pk.bytes, msg, ctx, &sig.bytes).unwrap();
//! ```
//!
//! ## Sizes
//!
//! | Variant    | PK (bytes) | SK (bytes) | Sig (bytes) |
//! |------------|------------|------------|-------------|
//! | PRISM-128  | 1312       | 4032       | 2420        |
//! | PRISM-192  | 1952       | 4032       | 3309        |
//! | PRISM-256  | 2592       | 4896       | 4627        |
//!
//! ## Implementation Notes
//!
//! - No `unsafe` code
//! - Constant-time signing (FIS + CT selection with `subtle` crate)
//! - Verification is constant-time in the comparison step
//! - No secret-dependent branches in core signing path (TODO: audit NTT)
//! - `zeroize` for secret key cleanup on drop (TODO: implement Zeroize for SecretKey)

pub mod error;
pub mod params;
pub mod poly;
pub mod ntt;
pub mod reduce;
pub mod sample;
pub mod hash;
pub mod packing;
pub mod keygen;
pub mod sign;
pub mod verify;

// Re-export the main API
pub use keygen::{keygen128, keygen128_from_seed, PublicKey128, SecretKey128};
pub use sign::{sign128, Signature128};
pub use verify::verify128;

#[cfg(feature = "prism192")]
pub use keygen::{keygen192, keygen192_from_seed, PublicKey192, SecretKey192};
#[cfg(feature = "prism192")]
pub use sign::{sign192, Signature192};
#[cfg(feature = "prism192")]
pub use verify::verify192;

#[cfg(feature = "prism256")]
pub use keygen::{keygen256, keygen256_from_seed, PublicKey256, SecretKey256};
#[cfg(feature = "prism256")]
pub use sign::{sign256, Signature256};
#[cfg(feature = "prism256")]
pub use verify::verify256;

pub use error::{KeyGenError, SignError, VerifyError};

/// Convenience wrapper: generate keys, sign, verify in one call.
/// Mainly useful for testing.
#[cfg(test)]
pub fn sign_and_verify_128(message: &[u8]) -> bool {
    use rand::rngs::OsRng;
    let (pk, sk) = keygen128(&mut OsRng).unwrap();
    let sig = sign128(&sk.bytes, message, b"test", &mut OsRng).unwrap();
    verify128(&pk.bytes, message, b"test", &sig.bytes).is_ok()
}
