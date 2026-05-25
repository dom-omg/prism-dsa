//! PRISM-DSA Key Generation
//!
//! Security: under Module-LWE and Module-SIS assumptions.
//! Same key generation algorithm as ML-DSA; keys are compatible.

use rand_core::{RngCore, CryptoRng};

use crate::params::*;
use crate::poly::Poly;
use crate::ntt::{ntt, invntt_tomont, matrix_vector_product};
use crate::sample::{expand_a, expand_s};
use crate::hash::shake256;
use crate::packing::*;
use crate::error::KeyGenError;

/// PRISM-128 public key (1312 bytes)
/// pk = (ρ, t1) where:
///   ρ: public seed for matrix A
///   t1: high bits of A·s1 + s2, packed in 10 bits per coefficient
pub struct PublicKey128 {
    pub bytes: [u8; p128::PK_BYTES],
}

/// PRISM-128 secret key (2528 bytes)
pub struct SecretKey128 {
    pub bytes: [u8; p128::SK_BYTES],
}

/// Generate a PRISM-128 key pair from OS randomness.
pub fn keygen128<R: RngCore + CryptoRng>(rng: &mut R) -> Result<(PublicKey128, SecretKey128), KeyGenError> {
    let mut xi = [0u8; SEED_BYTES];
    rng.try_fill_bytes(&mut xi).map_err(|_| KeyGenError::RngFailure)?;
    keygen128_from_seed(&xi)
}

/// Generate a PRISM-128 key pair from an explicit seed ξ.
/// Fully deterministic. For testing; production code should use keygen128().
pub fn keygen128_from_seed(xi: &[u8; SEED_BYTES]) -> Result<(PublicKey128, SecretKey128), KeyGenError> {
    use p128::*;

    // Derive ρ, ρ', key from seed
    let mut seed_expanded = [0u8; 2 * SEED_BYTES + CRH_BYTES];
    let mode_bytes = [K as u8, L as u8];
    shake256(&[xi as &[u8], &mode_bytes], &mut seed_expanded);

    let rho: [u8; SEED_BYTES] = seed_expanded[..SEED_BYTES].try_into().unwrap();
    let rho_prime: [u8; CRH_BYTES] = seed_expanded[SEED_BYTES..SEED_BYTES + CRH_BYTES].try_into().unwrap();
    let key: [u8; SEED_BYTES] = seed_expanded[SEED_BYTES + CRH_BYTES..].try_into().unwrap();

    // Expand matrix A (in coefficient domain; NTT done at signing)
    let mut a_mat = expand_a::<K, L>(&rho);

    // Sample secret vectors s1 ∈ S_η^l, s2 ∈ S_η^k
    let s1_coeffs = expand_s::<L>(&rho_prime, 0, ETA);
    let s2_coeffs = expand_s::<K>(&rho_prime, L as u16, ETA);

    // Compute t = A·s1 + s2
    // First NTT-transform A and s1
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() {
            ntt(p);
        }
    }

    let mut s1_ntt = s1_coeffs;
    for p in s1_ntt.iter_mut() { ntt(p); }

    // t_raw[i] = sum_j A[i][j] * s1[j] (in NTT domain)
    let mut t_raw = matrix_vector_product::<K, L>(&a_mat, &s1_ntt);

    // Back to coefficient domain and add s2
    for p in t_raw.iter_mut() {
        invntt_tomont(p);
    }
    for i in 0..K {
        for j in 0..N {
            t_raw[i].coeffs[j] += s2_coeffs[i].coeffs[j];
        }
        t_raw[i].caddq();
    }

    // Power2Round: t = t1·2^D + t0
    let mut t1_polys = [Poly::ZERO; K];
    let mut t0_polys = [Poly::ZERO; K];
    for i in 0..K {
        let (a1, a0) = t_raw[i].power2round();
        t1_polys[i] = a1;
        t0_polys[i] = a0;
    }

    // Pack public key: (ρ, t1)
    let mut pk_bytes = [0u8; PK_BYTES];
    pk_bytes[..SEED_BYTES].copy_from_slice(&rho);
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        pack_t1(&mut buf, &t1_polys[i]);
        pk_bytes[start..start + POLYT1_PACKED_BYTES].copy_from_slice(&buf);
    }

    // tr = SHAKE256(pk)
    let mut tr = [0u8; TR_BYTES];
    shake256(&[&pk_bytes as &[u8]], &mut tr);

    // Pack secret key: (ρ, key, tr, s1, s2, t0)
    let sk_size = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;

    let mut sk_bytes = vec![0u8; sk_size];
    let mut off = 0;

    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&rho);
    off += SEED_BYTES;
    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&key);
    off += SEED_BYTES;
    sk_bytes[off..off + TR_BYTES].copy_from_slice(&tr);
    off += TR_BYTES;

    for p in &s1_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta2(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &s2_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta2(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &t0_polys {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        pack_t0(&mut buf, p);
        sk_bytes[off..off + POLYT0_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYT0_PACKED_BYTES;
    }

    let mut sk_arr = [0u8; SK_BYTES];
    let actual_size = sk_size.min(SK_BYTES);
    sk_arr[..actual_size].copy_from_slice(&sk_bytes[..actual_size]);

    Ok((
        PublicKey128 { bytes: pk_bytes },
        SecretKey128 { bytes: sk_arr },
    ))
}

// ─────────────────────────────────────────────────────────
// PRISM-192 key generation
// ─────────────────────────────────────────────────────────

/// PRISM-192 public key
#[cfg(feature = "prism192")]
pub struct PublicKey192 {
    pub bytes: [u8; p192::PK_BYTES],
}

/// PRISM-192 secret key
#[cfg(feature = "prism192")]
pub struct SecretKey192 {
    pub bytes: [u8; p192::SK_BYTES],
}

/// Generate a PRISM-192 key pair from OS randomness.
#[cfg(feature = "prism192")]
pub fn keygen192<R: RngCore + CryptoRng>(rng: &mut R) -> Result<(PublicKey192, SecretKey192), KeyGenError> {
    let mut xi = [0u8; SEED_BYTES];
    rng.try_fill_bytes(&mut xi).map_err(|_| KeyGenError::RngFailure)?;
    keygen192_from_seed(&xi)
}

/// Generate a PRISM-192 key pair from an explicit seed ξ.
/// Fully deterministic. For testing; production code should use keygen192().
#[cfg(feature = "prism192")]
pub fn keygen192_from_seed(xi: &[u8; SEED_BYTES]) -> Result<(PublicKey192, SecretKey192), KeyGenError> {
    use p192::*;
    use crate::packing::pack_eta4;

    // Derive ρ, ρ', key from seed
    let mut seed_expanded = [0u8; 2 * SEED_BYTES + CRH_BYTES];
    let mode_bytes = [K as u8, L as u8];
    shake256(&[xi as &[u8], &mode_bytes], &mut seed_expanded);

    let rho: [u8; SEED_BYTES] = seed_expanded[..SEED_BYTES].try_into().unwrap();
    let rho_prime: [u8; CRH_BYTES] = seed_expanded[SEED_BYTES..SEED_BYTES + CRH_BYTES].try_into().unwrap();
    let key: [u8; SEED_BYTES] = seed_expanded[SEED_BYTES + CRH_BYTES..].try_into().unwrap();

    // Expand matrix A (in coefficient domain; NTT done at signing)
    let mut a_mat = expand_a::<K, L>(&rho);

    // Sample secret vectors s1 ∈ S_η^l, s2 ∈ S_η^k
    let s1_coeffs = expand_s::<L>(&rho_prime, 0, ETA);
    let s2_coeffs = expand_s::<K>(&rho_prime, L as u16, ETA);

    // Compute t = A·s1 + s2
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() {
            ntt(p);
        }
    }

    let mut s1_ntt = s1_coeffs;
    for p in s1_ntt.iter_mut() { ntt(p); }

    let mut t_raw = matrix_vector_product::<K, L>(&a_mat, &s1_ntt);

    for p in t_raw.iter_mut() {
        invntt_tomont(p);
    }
    for i in 0..K {
        for j in 0..N {
            t_raw[i].coeffs[j] += s2_coeffs[i].coeffs[j];
        }
        t_raw[i].caddq();
    }

    // Power2Round: t = t1·2^D + t0
    let mut t1_polys = [Poly::ZERO; K];
    let mut t0_polys = [Poly::ZERO; K];
    for i in 0..K {
        let (a1, a0) = t_raw[i].power2round();
        t1_polys[i] = a1;
        t0_polys[i] = a0;
    }

    // Pack public key: (ρ, t1)
    let mut pk_bytes = [0u8; PK_BYTES];
    pk_bytes[..SEED_BYTES].copy_from_slice(&rho);
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        pack_t1(&mut buf, &t1_polys[i]);
        pk_bytes[start..start + POLYT1_PACKED_BYTES].copy_from_slice(&buf);
    }

    // tr = SHAKE256(pk)
    let mut tr = [0u8; TR_BYTES];
    shake256(&[&pk_bytes as &[u8]], &mut tr);

    // Pack secret key: (ρ, key, tr, s1, s2, t0)
    let sk_size = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;

    let mut sk_bytes = vec![0u8; sk_size];
    let mut off = 0;

    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&rho);
    off += SEED_BYTES;
    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&key);
    off += SEED_BYTES;
    sk_bytes[off..off + TR_BYTES].copy_from_slice(&tr);
    off += TR_BYTES;

    for p in &s1_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta4(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &s2_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta4(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &t0_polys {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        pack_t0(&mut buf, p);
        sk_bytes[off..off + POLYT0_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYT0_PACKED_BYTES;
    }

    let mut sk_arr = [0u8; SK_BYTES];
    let actual_size = sk_size.min(SK_BYTES);
    sk_arr[..actual_size].copy_from_slice(&sk_bytes[..actual_size]);

    Ok((
        PublicKey192 { bytes: pk_bytes },
        SecretKey192 { bytes: sk_arr },
    ))
}

// ─────────────────────────────────────────────────────────
// PRISM-256 key generation
// ─────────────────────────────────────────────────────────

/// PRISM-256 public key
#[cfg(feature = "prism256")]
pub struct PublicKey256 {
    pub bytes: [u8; p256::PK_BYTES],
}

/// PRISM-256 secret key
#[cfg(feature = "prism256")]
pub struct SecretKey256 {
    pub bytes: [u8; p256::SK_BYTES],
}

/// Generate a PRISM-256 key pair from OS randomness.
#[cfg(feature = "prism256")]
pub fn keygen256<R: RngCore + CryptoRng>(rng: &mut R) -> Result<(PublicKey256, SecretKey256), KeyGenError> {
    let mut xi = [0u8; SEED_BYTES];
    rng.try_fill_bytes(&mut xi).map_err(|_| KeyGenError::RngFailure)?;
    keygen256_from_seed(&xi)
}

/// Generate a PRISM-256 key pair from an explicit seed ξ.
/// Fully deterministic. For testing; production code should use keygen256().
#[cfg(feature = "prism256")]
pub fn keygen256_from_seed(xi: &[u8; SEED_BYTES]) -> Result<(PublicKey256, SecretKey256), KeyGenError> {
    use p256::*;

    // Derive ρ, ρ', key from seed
    let mut seed_expanded = [0u8; 2 * SEED_BYTES + CRH_BYTES];
    let mode_bytes = [K as u8, L as u8];
    shake256(&[xi as &[u8], &mode_bytes], &mut seed_expanded);

    let rho: [u8; SEED_BYTES] = seed_expanded[..SEED_BYTES].try_into().unwrap();
    let rho_prime: [u8; CRH_BYTES] = seed_expanded[SEED_BYTES..SEED_BYTES + CRH_BYTES].try_into().unwrap();
    let key: [u8; SEED_BYTES] = seed_expanded[SEED_BYTES + CRH_BYTES..].try_into().unwrap();

    // Expand matrix A (in coefficient domain; NTT done at signing)
    let mut a_mat = expand_a::<K, L>(&rho);

    // Sample secret vectors s1 ∈ S_η^l, s2 ∈ S_η^k
    let s1_coeffs = expand_s::<L>(&rho_prime, 0, ETA);
    let s2_coeffs = expand_s::<K>(&rho_prime, L as u16, ETA);

    // Compute t = A·s1 + s2
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() {
            ntt(p);
        }
    }

    let mut s1_ntt = s1_coeffs;
    for p in s1_ntt.iter_mut() { ntt(p); }

    let mut t_raw = matrix_vector_product::<K, L>(&a_mat, &s1_ntt);

    for p in t_raw.iter_mut() {
        invntt_tomont(p);
    }
    for i in 0..K {
        for j in 0..N {
            t_raw[i].coeffs[j] += s2_coeffs[i].coeffs[j];
        }
        t_raw[i].caddq();
    }

    // Power2Round: t = t1·2^D + t0
    let mut t1_polys = [Poly::ZERO; K];
    let mut t0_polys = [Poly::ZERO; K];
    for i in 0..K {
        let (a1, a0) = t_raw[i].power2round();
        t1_polys[i] = a1;
        t0_polys[i] = a0;
    }

    // Pack public key: (ρ, t1)
    let mut pk_bytes = [0u8; PK_BYTES];
    pk_bytes[..SEED_BYTES].copy_from_slice(&rho);
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        pack_t1(&mut buf, &t1_polys[i]);
        pk_bytes[start..start + POLYT1_PACKED_BYTES].copy_from_slice(&buf);
    }

    // tr = SHAKE256(pk)
    let mut tr = [0u8; TR_BYTES];
    shake256(&[&pk_bytes as &[u8]], &mut tr);

    // Pack secret key: (ρ, key, tr, s1, s2, t0)
    let sk_size = 2 * SEED_BYTES + TR_BYTES
        + L * POLYETA_PACKED_BYTES
        + K * POLYETA_PACKED_BYTES
        + K * POLYT0_PACKED_BYTES;

    let mut sk_bytes = vec![0u8; sk_size];
    let mut off = 0;

    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&rho);
    off += SEED_BYTES;
    sk_bytes[off..off + SEED_BYTES].copy_from_slice(&key);
    off += SEED_BYTES;
    sk_bytes[off..off + TR_BYTES].copy_from_slice(&tr);
    off += TR_BYTES;

    for p in &s1_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta2(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &s2_coeffs {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        pack_eta2(&mut buf, p);
        sk_bytes[off..off + POLYETA_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYETA_PACKED_BYTES;
    }
    for p in &t0_polys {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        pack_t0(&mut buf, p);
        sk_bytes[off..off + POLYT0_PACKED_BYTES].copy_from_slice(&buf);
        off += POLYT0_PACKED_BYTES;
    }

    let mut sk_arr = [0u8; SK_BYTES];
    let actual_size = sk_size.min(SK_BYTES);
    sk_arr[..actual_size].copy_from_slice(&sk_bytes[..actual_size]);

    Ok((
        PublicKey256 { bytes: pk_bytes },
        SecretKey256 { bytes: sk_arr },
    ))
}
