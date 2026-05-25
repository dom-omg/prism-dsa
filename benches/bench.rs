use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::rngs::OsRng;
use prism_dsa::{keygen128, keygen128_from_seed, sign128, verify128};

// ── PRISM-DSA benchmarks ──────────────────────────────────────────────────────

fn bench_keygen128(c: &mut Criterion) {
    c.bench_function("prism128/keygen", |b| {
        b.iter(|| black_box(keygen128(&mut OsRng).unwrap()))
    });
}

fn bench_sign128(c: &mut Criterion) {
    let (_, sk) = keygen128_from_seed(&[0x42u8; 32]).unwrap();
    let msg = b"benchmark message";
    c.bench_function("prism128/sign_fis64", |b| {
        b.iter(|| black_box(sign128(&sk.bytes, msg, b"bench", &mut OsRng).unwrap()))
    });
}

fn bench_verify128(c: &mut Criterion) {
    let (pk, sk) = keygen128_from_seed(&[0x42u8; 32]).unwrap();
    let msg = b"benchmark message";
    let sig = sign128(&sk.bytes, msg, b"bench", &mut OsRng).unwrap();
    c.bench_function("prism128/verify", |b| {
        b.iter(|| black_box(verify128(&pk.bytes, msg, b"bench", &sig.bytes).unwrap()))
    });
}

fn bench_sign_verify128(c: &mut Criterion) {
    let (pk, sk) = keygen128(&mut OsRng).unwrap();
    let msg = b"sign+verify benchmark";
    c.bench_function("prism128/sign+verify", |b| {
        b.iter(|| {
            let sig = sign128(&sk.bytes, msg, b"bench", &mut OsRng).unwrap();
            black_box(verify128(&pk.bytes, msg, b"bench", &sig.bytes).unwrap())
        })
    });
}

#[cfg(feature = "prism192")]
fn bench_sign192(c: &mut Criterion) {
    use prism_dsa::{keygen192_from_seed, sign192, verify192};
    let (pk, sk) = keygen192_from_seed(&[0x42u8; 32]).unwrap();
    let msg = b"benchmark message 192";
    let sig = sign192(&sk.bytes, msg, b"bench", &mut OsRng).unwrap();
    c.bench_function("prism192/sign_fis64", |b| {
        b.iter(|| black_box(sign192(&sk.bytes, msg, b"bench", &mut OsRng).unwrap()))
    });
    c.bench_function("prism192/verify", |b| {
        b.iter(|| black_box(verify192(&pk.bytes, msg, b"bench", &sig.bytes).unwrap()))
    });
}

#[cfg(not(feature = "prism192"))]
fn bench_sign192(_c: &mut Criterion) {}

#[cfg(feature = "prism256")]
fn bench_sign256(c: &mut Criterion) {
    use prism_dsa::{keygen256_from_seed, sign256, verify256};
    let (pk, sk) = keygen256_from_seed(&[0x42u8; 32]).unwrap();
    let msg = b"benchmark message 256";
    let sig = sign256(&sk.bytes, msg, b"bench", &mut OsRng).unwrap();
    c.bench_function("prism256/sign_fis64", |b| {
        b.iter(|| black_box(sign256(&sk.bytes, msg, b"bench", &mut OsRng).unwrap()))
    });
    c.bench_function("prism256/verify", |b| {
        b.iter(|| black_box(verify256(&pk.bytes, msg, b"bench", &sig.bytes).unwrap()))
    });
}

#[cfg(not(feature = "prism256"))]
fn bench_sign256(_c: &mut Criterion) {}

// ── ML-DSA comparison benchmarks (FIPS 204 reference) ────────────────────────

fn bench_mldsa44(c: &mut Criterion) {
    use ml_dsa::{MlDsa44, Keypair, Signer, Verifier, Generate};
    let sk = ml_dsa::SigningKey::<MlDsa44>::generate();
    let vk = sk.verifying_key();
    let msg = b"benchmark message";
    let sig = sk.sign(msg);

    c.bench_function("mldsa44/keygen", |b| {
        b.iter(|| black_box(ml_dsa::SigningKey::<MlDsa44>::generate()))
    });
    c.bench_function("mldsa44/sign", |b| {
        b.iter(|| black_box(sk.sign(msg)))
    });
    c.bench_function("mldsa44/verify", |b| {
        b.iter(|| black_box(vk.verify(msg, &sig).unwrap()))
    });
}

fn bench_mldsa65(c: &mut Criterion) {
    use ml_dsa::{MlDsa65, Keypair, Signer, Verifier, Generate};
    let sk = ml_dsa::SigningKey::<MlDsa65>::generate();
    let vk = sk.verifying_key();
    let msg = b"benchmark message";
    let sig = sk.sign(msg);

    c.bench_function("mldsa65/keygen", |b| {
        b.iter(|| black_box(ml_dsa::SigningKey::<MlDsa65>::generate()))
    });
    c.bench_function("mldsa65/sign", |b| {
        b.iter(|| black_box(sk.sign(msg)))
    });
    c.bench_function("mldsa65/verify", |b| {
        b.iter(|| black_box(vk.verify(msg, &sig).unwrap()))
    });
}

fn bench_mldsa87(c: &mut Criterion) {
    use ml_dsa::{MlDsa87, Keypair, Signer, Verifier, Generate};
    let sk = ml_dsa::SigningKey::<MlDsa87>::generate();
    let vk = sk.verifying_key();
    let msg = b"benchmark message";
    let sig = sk.sign(msg);

    c.bench_function("mldsa87/keygen", |b| {
        b.iter(|| black_box(ml_dsa::SigningKey::<MlDsa87>::generate()))
    });
    c.bench_function("mldsa87/sign", |b| {
        b.iter(|| black_box(sk.sign(msg)))
    });
    c.bench_function("mldsa87/verify", |b| {
        b.iter(|| black_box(vk.verify(msg, &sig).unwrap()))
    });
}

criterion_group!(
    prism_benches,
    bench_keygen128,
    bench_sign128,
    bench_verify128,
    bench_sign_verify128,
    bench_sign192,
    bench_sign256
);
criterion_group!(
    mldsa_benches,
    bench_mldsa44,
    bench_mldsa65,
    bench_mldsa87
);
criterion_main!(prism_benches, mldsa_benches);
