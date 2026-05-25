//! Known Answer Tests for PRISM-DSA
//!
//! Tests against pre-generated KAT vectors in tests/kat/prism128_kat.txt
//! Generated with ZeroRng (rnd = 0^32) — fully deterministic.
//!
//! To verify against NIST ML-DSA-44 vectors:
//!   Place sigVer.rsp from NIST ACVP at tests/kat/ml_dsa44_sigVer.rsp
//!   Format: pk=<hex> message=<hex> context=<hex> signature=<hex> testPassed=true

use prism_dsa::{keygen128_from_seed, sign128, verify128};
use rand_core::{RngCore, CryptoRng};
use std::collections::HashMap;

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

fn from_hex(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}

fn parse_kat_file(path: &str) -> Vec<HashMap<String, String>> {
    let content = std::fs::read_to_string(path).expect("KAT file missing");
    let mut records = Vec::new();
    let mut current = HashMap::new();

    for line in content.lines() {
        // Trim only the right side so "message = " (trailing space) stays parseable.
        let line = line.trim_end();
        let line_trimmed = line.trim();
        if line_trimmed.is_empty() || line_trimmed.starts_with('#') {
            if current.contains_key("seed") {
                records.push(std::mem::take(&mut current));
            }
            continue;
        }
        // Split on " = " first; if that fails, try " =" to handle "key =" with empty value.
        if let Some((k, v)) = line_trimmed.split_once(" = ") {
            current.insert(k.trim().to_string(), v.trim().to_string());
        } else if let Some(k) = line_trimmed.strip_suffix(" =") {
            current.insert(k.trim().to_string(), String::new());
        }
    }
    if current.contains_key("seed") {
        records.push(current);
    }
    records
}

/// Load and verify all KAT vectors from the pre-generated file.
/// For each vector: re-derive pk from seed, sign with ZeroRng, check sig matches stored sig.
#[test]
fn kat_prism128_from_file() {
    let kat_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/kat/prism128_kat.txt");
    let records = parse_kat_file(kat_path);
    assert!(!records.is_empty(), "KAT file must have at least one record");

    for (i, rec) in records.iter().enumerate() {
        let seed_hex = rec.get("seed").unwrap();
        let msg_hex = rec.get("message").unwrap();
        let ctx_hex = rec.get("context").unwrap();
        let sig_hex = rec.get("sig").unwrap();

        let seed: [u8; 32] = from_hex(seed_hex).try_into()
            .unwrap_or_else(|_| panic!("bad seed at record {}", i));
        let msg = from_hex(msg_hex);
        let ctx = from_hex(ctx_hex);
        let expected_sig = from_hex(sig_hex);

        // Re-derive key
        let (pk, sk) = keygen128_from_seed(&seed).unwrap();

        // Sign with ZeroRng — must match stored vector
        let sig = sign128(&sk.bytes, &msg, &ctx, &mut ZeroRng).unwrap();
        assert_eq!(
            sig.bytes.as_ref(),
            expected_sig.as_slice(),
            "signature mismatch at KAT record {} (seed={}...)",
            i, &seed_hex[..8]
        );

        // Verify stored sig against re-derived pk
        let sig_arr: [u8; prism_dsa::params::p128::SIG_BYTES] = expected_sig
            .try_into()
            .unwrap_or_else(|_| panic!("wrong sig length at record {}", i));
        assert!(
            verify128(&pk.bytes, &msg, &ctx, &sig_arr).is_ok(),
            "stored sig failed verification at record {}", i
        );
    }

    eprintln!("KAT PRISM-128: {} vectors verified ✓", records.len());
}

/// Self-consistency: sign with ZeroRng, then verify — all seeds must round-trip.
#[test]
fn kat_self_consistency_exhaustive() {
    let test_cases: &[([u8; 32], &[u8], &[u8])] = &[
        ([0x00; 32], b"",                    b""),
        ([0x01; 32], b"test",                b"ctx"),
        ([0x42; 32], b"KAT message",         b"kat-v1"),
        ([0xFF; 32], b"The quick brown fox", b"prism128"),
        ([0xAB; 32], b"PRISM-DSA prototype", b"research"),
        ([0x7F; 32], b"edge case",           b""),
        ([0x80; 32], b"",                    b"empty-msg"),
    ];

    for (seed, msg, ctx) in test_cases {
        let (pk, sk) = keygen128_from_seed(seed).unwrap();
        let sig = sign128(&sk.bytes, msg, ctx, &mut ZeroRng).unwrap();
        assert!(
            verify128(&pk.bytes, msg, ctx, &sig.bytes).is_ok(),
            "self-consistency failed for seed {:02x}...", seed[0]
        );

        // Wrong message must fail
        let bad_msg = b"WRONG_MESSAGE_DO_NOT_VERIFY";
        if *msg != bad_msg.as_slice() {
            assert!(
                verify128(&pk.bytes, bad_msg, ctx, &sig.bytes).is_err(),
                "wrong message should fail for seed {:02x}...", seed[0]
            );
        }
    }
}

/// Signature size must match the declared constant.
#[test]
fn kat_signature_size_invariant() {
    use prism_dsa::params::p128;
    let (_, sk) = keygen128_from_seed(&[0x42; 32]).unwrap();
    let sig = sign128(&sk.bytes, b"size check", b"sz", &mut ZeroRng).unwrap();
    assert_eq!(sig.bytes.len(), p128::SIG_BYTES,
        "signature size must equal p128::SIG_BYTES = {}", p128::SIG_BYTES);
}
