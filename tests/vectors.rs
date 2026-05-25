//! Correctness tests for PRISM-DSA
//!
//! Coverage:
//! - PRISM-128 sign + verify round-trip (10 tests)
//! - KAT determinism: same inputs → same signature bytes
//! - Differential: wrong message/context/key all fail
//! - PRISM-192 sign + verify round-trip (feature-gated)
//! - PRISM-256 sign + verify round-trip (feature-gated)
//! - FIS completeness: 50 consecutive signs succeed

use prism_dsa::{keygen128_from_seed, sign128, verify128, VerifyError};
use rand_core::{RngCore, CryptoRng};

/// Deterministic zero RNG for KAT tests
struct ZeroRng;
impl RngCore for ZeroRng {
    fn next_u32(&mut self) -> u32 { 0 }
    fn next_u64(&mut self) -> u64 { 0 }
    fn fill_bytes(&mut self, dest: &mut [u8]) { dest.fill(0); }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        dest.fill(0); Ok(())
    }
}
impl CryptoRng for ZeroRng {}

// ─────────────────────────────────────────────────────────
// PRISM-128 tests
// ─────────────────────────────────────────────────────────

#[test]
fn sign_verify_roundtrip() {
    let seed = [0x42u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    let message = b"PRISM-DSA test message";
    let context = b"test-context-v1";
    let sig = sign128(&sk.bytes, message, context, &mut ZeroRng).unwrap();
    assert!(verify128(&pk.bytes, message, context, &sig.bytes).is_ok());
}

#[test]
fn wrong_message_fails() {
    let seed = [0x01u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    let sig = sign128(&sk.bytes, b"correct", b"ctx", &mut ZeroRng).unwrap();
    assert_eq!(verify128(&pk.bytes, b"wrong", b"ctx", &sig.bytes).unwrap_err(), VerifyError::Forgery);
}

#[test]
fn wrong_context_fails() {
    let seed = [0x02u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    let sig = sign128(&sk.bytes, b"msg", b"ctx-v1", &mut ZeroRng).unwrap();
    assert_eq!(verify128(&pk.bytes, b"msg", b"ctx-v2", &sig.bytes).unwrap_err(), VerifyError::Forgery);
}

#[test]
fn keygen_deterministic() {
    let seed = [0xABu8; 32];
    let (pk1, sk1) = keygen128_from_seed(&seed).unwrap();
    let (pk2, sk2) = keygen128_from_seed(&seed).unwrap();
    assert_eq!(pk1.bytes, pk2.bytes);
    assert_eq!(sk1.bytes, sk2.bytes);
}

#[test]
fn flipped_ctilde_fails() {
    let seed = [0x03u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    let mut sig = sign128(&sk.bytes, b"test", b"ctx", &mut ZeroRng).unwrap();
    sig.bytes[0] ^= 0x01;
    assert!(verify128(&pk.bytes, b"test", b"ctx", &sig.bytes).is_err());
}

#[test]
fn context_too_long() {
    let (_, sk) = keygen128_from_seed(&[0x04u8; 32]).unwrap();
    let ctx = vec![0u8; 256];
    assert!(sign128(&sk.bytes, b"msg", &ctx, &mut ZeroRng).is_err());
}

#[test]
fn multiple_messages() {
    let seed = [0x05u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    for i in 0u8..10 {
        let msg = vec![i; 100];
        let sig = sign128(&sk.bytes, &msg, b"app", &mut ZeroRng).unwrap();
        assert!(verify128(&pk.bytes, &msg, b"app", &sig.bytes).is_ok());
    }
}

#[test]
fn empty_message() {
    let (pk, sk) = keygen128_from_seed(&[0x06u8; 32]).unwrap();
    let sig = sign128(&sk.bytes, b"", b"", &mut ZeroRng).unwrap();
    assert!(verify128(&pk.bytes, b"", b"", &sig.bytes).is_ok());
}

#[test]
fn large_message() {
    let (pk, sk) = keygen128_from_seed(&[0x07u8; 32]).unwrap();
    let msg = vec![0xFEu8; 10_000];
    let sig = sign128(&sk.bytes, &msg, b"", &mut ZeroRng).unwrap();
    assert!(verify128(&pk.bytes, &msg, b"", &sig.bytes).is_ok());
}

#[test]
fn wrong_key_fails() {
    let (_, sk1) = keygen128_from_seed(&[0x10u8; 32]).unwrap();
    let (pk2, _) = keygen128_from_seed(&[0x11u8; 32]).unwrap();
    let sig = sign128(&sk1.bytes, b"test", b"", &mut ZeroRng).unwrap();
    assert!(verify128(&pk2.bytes, b"test", b"", &sig.bytes).is_err());
}

/// KAT: same inputs must produce identical signature bytes (determinism with ZeroRng)
#[test]
fn kat_deterministic_signature() {
    let seed = [0x42u8; 32];
    let (pk, sk) = keygen128_from_seed(&seed).unwrap();
    let msg = b"KAT test message";
    let ctx = b"kat-v1";
    let sig1 = sign128(&sk.bytes, msg, ctx, &mut ZeroRng).unwrap();
    let sig2 = sign128(&sk.bytes, msg, ctx, &mut ZeroRng).unwrap();
    assert_eq!(sig1.bytes, sig2.bytes, "signing must be deterministic with ZeroRng");
    assert!(verify128(&pk.bytes, msg, ctx, &sig1.bytes).is_ok());
}

/// Differential: many seeds × many messages, all must round-trip
#[test]
fn differential_sign_verify() {
    let seeds: [[u8; 32]; 5] = [[1;32], [2;32], [3;32], [100;32], [255;32]];
    let messages: &[&[u8]] = &[b"msg1", b"", &[0xFFu8; 100]];
    for seed in &seeds {
        let (pk, sk) = keygen128_from_seed(seed).unwrap();
        for msg in messages {
            let sig = sign128(&sk.bytes, msg, b"diff-test", &mut ZeroRng).unwrap();
            assert!(verify128(&pk.bytes, msg, b"diff-test", &sig.bytes).is_ok());
        }
    }
}

/// FIS completeness: 50 consecutive sign calls must all succeed (P(fail) ≈ 2^{-27} each)
#[test]
fn fis_never_fails_in_practice() {
    use rand::rngs::OsRng;
    let (pk, sk) = keygen128_from_seed(&[0x99u8; 32]).unwrap();
    for i in 0u8..50 {
        let msg = [i; 32];
        let sig = sign128(&sk.bytes, &msg, b"fis-test", &mut OsRng).unwrap();
        assert!(verify128(&pk.bytes, &msg, b"fis-test", &sig.bytes).is_ok());
    }
}

// ─────────────────────────────────────────────────────────
// PRISM-192 tests
// ─────────────────────────────────────────────────────────

#[cfg(feature = "prism192")]
mod prism192_tests {
    use prism_dsa::{keygen192_from_seed, sign192, verify192, VerifyError};
    use super::ZeroRng;

    #[test]
    fn sign_verify_roundtrip_192() {
        let (pk, sk) = keygen192_from_seed(&[0x42u8; 32]).unwrap();
        let sig = sign192(&sk.bytes, b"test message 192", b"ctx", &mut ZeroRng).unwrap();
        assert!(verify192(&pk.bytes, b"test message 192", b"ctx", &sig.bytes).is_ok());
    }

    #[test]
    fn wrong_message_fails_192() {
        let (pk, sk) = keygen192_from_seed(&[0x10u8; 32]).unwrap();
        let sig = sign192(&sk.bytes, b"correct", b"ctx", &mut ZeroRng).unwrap();
        assert_eq!(verify192(&pk.bytes, b"wrong", b"ctx", &sig.bytes).unwrap_err(), VerifyError::Forgery);
    }

    #[test]
    fn wrong_key_fails_192() {
        let (_, sk1) = keygen192_from_seed(&[0x20u8; 32]).unwrap();
        let (pk2, _) = keygen192_from_seed(&[0x21u8; 32]).unwrap();
        let sig = sign192(&sk1.bytes, b"msg", b"", &mut ZeroRng).unwrap();
        assert!(verify192(&pk2.bytes, b"msg", b"", &sig.bytes).is_err());
    }

    #[test]
    fn multiple_messages_192() {
        let (pk, sk) = keygen192_from_seed(&[0x30u8; 32]).unwrap();
        for i in 0u8..5 {
            let msg = [i; 64];
            let sig = sign192(&sk.bytes, &msg, b"app192", &mut ZeroRng).unwrap();
            assert!(verify192(&pk.bytes, &msg, b"app192", &sig.bytes).is_ok());
        }
    }

    #[test]
    fn kat_deterministic_192() {
        let (pk, sk) = keygen192_from_seed(&[0x42u8; 32]).unwrap();
        let sig1 = sign192(&sk.bytes, b"kat", b"v1", &mut ZeroRng).unwrap();
        let sig2 = sign192(&sk.bytes, b"kat", b"v1", &mut ZeroRng).unwrap();
        assert_eq!(sig1.bytes, sig2.bytes);
        assert!(verify192(&pk.bytes, b"kat", b"v1", &sig1.bytes).is_ok());
    }
}

// ─────────────────────────────────────────────────────────
// PRISM-256 tests
// ─────────────────────────────────────────────────────────

#[cfg(feature = "prism256")]
mod prism256_tests {
    use prism_dsa::{keygen256_from_seed, sign256, verify256, VerifyError};
    use super::ZeroRng;

    #[test]
    fn sign_verify_roundtrip_256() {
        let (pk, sk) = keygen256_from_seed(&[0x42u8; 32]).unwrap();
        let sig = sign256(&sk.bytes, b"test message 256", b"ctx", &mut ZeroRng).unwrap();
        assert!(verify256(&pk.bytes, b"test message 256", b"ctx", &sig.bytes).is_ok());
    }

    #[test]
    fn wrong_message_fails_256() {
        let (pk, sk) = keygen256_from_seed(&[0x10u8; 32]).unwrap();
        let sig = sign256(&sk.bytes, b"correct", b"ctx", &mut ZeroRng).unwrap();
        assert_eq!(verify256(&pk.bytes, b"wrong", b"ctx", &sig.bytes).unwrap_err(), VerifyError::Forgery);
    }

    #[test]
    fn wrong_key_fails_256() {
        let (_, sk1) = keygen256_from_seed(&[0x20u8; 32]).unwrap();
        let (pk2, _) = keygen256_from_seed(&[0x21u8; 32]).unwrap();
        let sig = sign256(&sk1.bytes, b"msg", b"", &mut ZeroRng).unwrap();
        assert!(verify256(&pk2.bytes, b"msg", b"", &sig.bytes).is_err());
    }

    #[test]
    fn multiple_messages_256() {
        let (pk, sk) = keygen256_from_seed(&[0x30u8; 32]).unwrap();
        for i in 0u8..5 {
            let msg = [i; 64];
            let sig = sign256(&sk.bytes, &msg, b"app256", &mut ZeroRng).unwrap();
            assert!(verify256(&pk.bytes, &msg, b"app256", &sig.bytes).is_ok());
        }
    }

    #[test]
    fn kat_deterministic_256() {
        let (pk, sk) = keygen256_from_seed(&[0x42u8; 32]).unwrap();
        let sig1 = sign256(&sk.bytes, b"kat", b"v1", &mut ZeroRng).unwrap();
        let sig2 = sign256(&sk.bytes, b"kat", b"v1", &mut ZeroRng).unwrap();
        assert_eq!(sig1.bytes, sig2.bytes);
        assert!(verify256(&pk.bytes, b"kat", b"v1", &sig1.bytes).is_ok());
    }
}
