// KAT vector generator for PRISM-DSA
// Run: cargo run --features prism128 --bin gen_kat
// Outputs deterministic test vectors using ZeroRng (rnd = [0u8; 32])

use rand_core::{RngCore, CryptoRng};
use prism_dsa::{keygen128_from_seed, sign128, verify128};

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

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{:02x}", x)).collect()
}

fn main() {
    println!("# PRISM-DSA PRISM-128 Known Answer Tests");
    println!("# Generated with ZeroRng (rnd = 0^32)");
    println!("# Format: count / seed / message / context / pk / sk_hash / sig / verified");
    println!();

    let test_cases = [
        ([0x00u8; 32], b"".as_slice(),          b"".as_slice()),
        ([0x01u8; 32], b"test",                  b"ctx"),
        ([0x42u8; 32], b"KAT message",           b"kat-v1"),
        ([0xFFu8; 32], b"The quick brown fox",   b"prism128"),
        ([0xABu8; 32], b"PRISM-DSA prototype",   b"research"),
    ];

    for (i, (seed, msg, ctx)) in test_cases.iter().enumerate() {
        let (pk, sk) = keygen128_from_seed(seed).unwrap();
        let sig = sign128(&sk.bytes, msg, ctx, &mut ZeroRng).unwrap();
        let verified = verify128(&pk.bytes, msg, ctx, &sig.bytes).is_ok();

        // Use SHA3-256 of sk as compact representation (sk is 4032 bytes)
        // For the KAT we store the full pk and sig
        println!("count = {}", i);
        println!("seed = {}", hex(seed));
        println!("message = {}", hex(msg));
        println!("context = {}", hex(ctx));
        println!("pk = {}", hex(&pk.bytes));
        println!("sig = {}", hex(&sig.bytes));
        println!("verified = {}", verified);
        println!();
    }
}
