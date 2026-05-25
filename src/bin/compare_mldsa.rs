use std::time::{Duration, Instant};
use std::hint::black_box;

use ml_dsa::{MlDsa44, MlDsa65, MlDsa87, Keypair, Signer, Verifier, Generate};
use prism_dsa::{keygen128_from_seed, sign128, verify128};
use rand::rngs::OsRng;

const WARMUP: usize = 5;
const ITERS: usize = 100;

fn time_fn<F: FnMut() -> ()>(mut f: F) -> Duration {
    for _ in 0..WARMUP { f(); }
    let t = Instant::now();
    for _ in 0..ITERS { f(); }
    t.elapsed() / ITERS as u32
}

fn us(d: Duration) -> f64 { d.as_secs_f64() * 1_000_000.0 }

fn main() {
    println!("\n╔═══════════════════════════════════════════════════════════════════╗");
    println!("║  PRISM-DSA vs ML-DSA (FIPS 204) — Performance Comparison          ║");
    println!("║  {} warmup + {} measured iterations, ARM64 Linux                ║", WARMUP, ITERS);
    println!("╚═══════════════════════════════════════════════════════════════════╝\n");

    println!("{:<30} {:>12} {:>12} {:>12}", "Algorithm", "keygen (µs)", "sign (µs)", "verify (µs)");
    println!("{}", "─".repeat(70));

    // ── PRISM-128 ─────────────────────────────────────────────────────────────
    let (pk128, sk128) = keygen128_from_seed(&[0x42u8; 32]).unwrap();
    let msg = b"benchmark message";
    let sig128 = sign128(&sk128.bytes, msg, b"bench", &mut OsRng).unwrap();

    let t_kg = time_fn(|| { black_box(prism_dsa::keygen128(&mut OsRng).unwrap()); });
    let t_sg = time_fn(|| { black_box(sign128(&sk128.bytes, msg, b"bench", &mut OsRng).unwrap()); });
    let t_vf = time_fn(|| { black_box(verify128(&pk128.bytes, msg, b"bench", &sig128.bytes).unwrap()); });
    println!("{:<30} {:>12.1} {:>12.1} {:>12.1}",
        "PRISM-128 (CT, 64-slot FIS)", us(t_kg), us(t_sg), us(t_vf));

    // ── ML-DSA-44 (NIST Level 2, K=4 L=4) ───────────────────────────────────
    let sk44 = ml_dsa::SigningKey::<MlDsa44>::generate();
    let vk44 = sk44.verifying_key();
    let sig44 = sk44.sign(msg);

    let t_kg = time_fn(|| { black_box(ml_dsa::SigningKey::<MlDsa44>::generate()); });
    let t_sg = time_fn(|| { black_box(sk44.sign(msg)); });
    let t_vf = time_fn(|| { black_box(vk44.verify(msg, &sig44).unwrap()); });
    println!("{:<30} {:>12.1} {:>12.1} {:>12.1}",
        "ML-DSA-44 (non-CT, ~1.4 iter)", us(t_kg), us(t_sg), us(t_vf));

    // ── ML-DSA-65 (NIST Level 3, K=6 L=5) ───────────────────────────────────
    let sk65 = ml_dsa::SigningKey::<MlDsa65>::generate();
    let vk65 = sk65.verifying_key();
    let sig65 = sk65.sign(msg);

    let t_kg = time_fn(|| { black_box(ml_dsa::SigningKey::<MlDsa65>::generate()); });
    let t_sg = time_fn(|| { black_box(sk65.sign(msg)); });
    let t_vf = time_fn(|| { black_box(vk65.verify(msg, &sig65).unwrap()); });
    println!("{:<30} {:>12.1} {:>12.1} {:>12.1}",
        "ML-DSA-65 (non-CT, ~1.4 iter)", us(t_kg), us(t_sg), us(t_vf));

    // ── ML-DSA-87 (NIST Level 5, K=8 L=7) ───────────────────────────────────
    let sk87 = ml_dsa::SigningKey::<MlDsa87>::generate();
    let vk87 = sk87.verifying_key();
    let sig87 = sk87.sign(msg);

    let t_kg = time_fn(|| { black_box(ml_dsa::SigningKey::<MlDsa87>::generate()); });
    let t_sg = time_fn(|| { black_box(sk87.sign(msg)); });
    let t_vf = time_fn(|| { black_box(vk87.verify(msg, &sig87).unwrap()); });
    println!("{:<30} {:>12.1} {:>12.1} {:>12.1}",
        "ML-DSA-87 (non-CT, ~1.4 iter)", us(t_kg), us(t_sg), us(t_vf));

    println!("{}", "─".repeat(70));

    // ── Ratio summary ─────────────────────────────────────────────────────────
    let prism128_sign = {
        for _ in 0..WARMUP { let _ = sign128(&sk128.bytes, msg, b"bench", &mut OsRng); }
        let t = Instant::now();
        for _ in 0..ITERS { let _ = sign128(&sk128.bytes, msg, b"bench", &mut OsRng); }
        t.elapsed() / ITERS as u32
    };
    let mldsa65_sign = {
        for _ in 0..WARMUP { let _ = sk65.sign(msg); }
        let t = Instant::now();
        for _ in 0..ITERS { let _ = sk65.sign(msg); }
        t.elapsed() / ITERS as u32
    };

    println!("\nPRISM-128 sign:   {:.3}ms (CT timing guarantee, exactly 64 FIS iterations)",
        prism128_sign.as_secs_f64() * 1000.0);
    println!("ML-DSA-65 sign:   {:.3}ms (non-CT, variable rejection sampling)",
        mldsa65_sign.as_secs_f64() * 1000.0);
    println!("Overhead factor:  {:.1}x  (cost of timing-invariant signing)",
        prism128_sign.as_secs_f64() / mldsa65_sign.as_secs_f64());
    println!("\nNote: PRISM-128 ≈ NIST Level 1 (same K=4, L=4 as ML-DSA-44)");
    println!("      Both implementations are pure Rust, no hand-optimized assembly.");
}
