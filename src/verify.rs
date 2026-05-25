//! PRISM-DSA Signature Verification
//!
//! Verification is identical to ML-DSA-44 verification.
//! FIS only affects the signing algorithm; verification is unchanged.

use crate::params::*;
use crate::poly::Poly;
use crate::ntt::{ntt, invntt_tomont, matrix_vector_product};
use crate::sample::{expand_a, sample_in_ball};
use crate::hash::shake256;
use crate::packing::*;
use crate::error::VerifyError;

/// Verify a PRISM-128 signature.
pub fn verify128(
    pk_bytes: &[u8; p128::PK_BYTES],
    message: &[u8],
    context: &[u8],
    sig_bytes: &[u8; p128::SIG_BYTES],
) -> Result<(), VerifyError> {
    use p128::*;

    if context.len() > 255 {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Unpack public key (ρ, t1) ─────────────────────────────────
    let rho: [u8; SEED_BYTES] = pk_bytes[..SEED_BYTES].try_into().unwrap();
    let mut t1 = [Poly::ZERO; K];
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        buf.copy_from_slice(&pk_bytes[start..start + POLYT1_PACKED_BYTES]);
        unpack_t1(&mut t1[i], &buf);
    }

    // ── Unpack signature (c̃, z, h) ──────────────────────────────
    let ctilde: [u8; LAMBDA_BYTES] = sig_bytes[..LAMBDA_BYTES].try_into().unwrap();
    let mut sig_off = LAMBDA_BYTES;

    let mut z = [Poly::ZERO; L];
    for i in 0..L {
        let mut buf = [0u8; POLYZ_PACKED_BYTES];
        buf.copy_from_slice(&sig_bytes[sig_off..sig_off + POLYZ_PACKED_BYTES]);
        unpack_z_gamma1_17(&mut z[i], &buf, GAMMA1);
        sig_off += POLYZ_PACKED_BYTES;
    }

    let hint_bytes = &sig_bytes[sig_off..sig_off + OMEGA + K];
    let mut h = [Poly::ZERO; K];
    if !unpack_hint::<K>(&mut h, hint_bytes, OMEGA) {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Validity checks ───────────────────────────────────────────
    // Check 1: ||z||_∞ < γ1 - β
    for p in &z {
        if !p.check_norm(GAMMA1 - BETA) {
            return Err(VerifyError::InvalidSignature);
        }
    }

    // Check 2: hint weight ≤ ω
    let h_weight: usize = h.iter().map(|p| p.popcount()).sum();
    if h_weight > OMEGA {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Compute μ = H(tr || prefix || message) ────────────────────
    let mut tr = [0u8; TR_BYTES];
    shake256(&[pk_bytes as &[u8]], &mut tr);
    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    // ── c = SampleInBall(c̃), NTT(c) ─────────────────────────────
    let mut c = sample_in_ball(&ctilde, TAU);
    ntt(&mut c);

    // ── Expand A and NTT-transform t1 ────────────────────────────
    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() { ntt(p); }
    }

    let mut t1_ntt = t1;
    for p in t1_ntt.iter_mut() { ntt(p); }

    // ── Compute w'_approx = Az - ct1·2^D ────────────────────────
    let mut z_ntt = z;
    for p in z_ntt.iter_mut() { ntt(p); }

    // Az
    let mut az = matrix_vector_product::<K, L>(&a_mat, &z_ntt);
    for p in az.iter_mut() { invntt_tomont(p); }

    // ct1
    let mut ct1 = [Poly::ZERO; K];
    for i in 0..K {
        ct1[i] = c.pointwise_montgomery(&t1_ntt[i]);
    }
    for p in ct1.iter_mut() { invntt_tomont(p); }

    // w'_approx = Az - ct1 * 2^D mod q, map to [0, q-1]
    // Note: ct1[j] can be ~q, so (ct1[j] << D) overflows i32.
    // Compute in i64 then reduce mod q.
    let mut w_approx = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            let val = (az[i].coeffs[j] as i64)
                - ((ct1[i].coeffs[j] as i64) << D);
            // Reduce to [0, q-1]
            let reduced = val.rem_euclid(Q as i64) as i32;
            w_approx[i].coeffs[j] = reduced;
        }
    }

    // ── w1_recovered = UseHint(h, w'_approx) ─────────────────────
    let mut w1_recovered = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            w1_recovered[i].coeffs[j] = crate::poly::use_hint_coeff(
                h[i].coeffs[j],
                w_approx[i].coeffs[j],
                GAMMA2,
            );
        }
    }

    // ── Recompute challenge hash and compare ─────────────────────
    let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
    for i in 0..K {
        let mut buf = [0u8; POLYW1_PACKED_BYTES];
        pack_w1_gamma2_88(&mut buf, &w1_recovered[i]);
        w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES]
            .copy_from_slice(&buf);
    }

    let mut ctilde_check = [0u8; LAMBDA_BYTES];
    shake256(&[&mu as &[u8], &w1_packed], &mut ctilde_check);

    // Constant-time comparison
    use subtle::ConstantTimeEq;
    if ctilde.as_ref().ct_eq(ctilde_check.as_ref()).unwrap_u8() == 0 {
        return Err(VerifyError::Forgery);
    }

    Ok(())
}

/// Verify a PRISM-192 signature.
#[cfg(feature = "prism192")]
pub fn verify192(
    pk_bytes: &[u8; p192::PK_BYTES],
    message: &[u8],
    context: &[u8],
    sig_bytes: &[u8; p192::SIG_BYTES],
) -> Result<(), VerifyError> {
    use p192::*;

    if context.len() > 255 {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Unpack public key (ρ, t1) ─────────────────────────────────
    let rho: [u8; SEED_BYTES] = pk_bytes[..SEED_BYTES].try_into().unwrap();
    let mut t1 = [Poly::ZERO; K];
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        buf.copy_from_slice(&pk_bytes[start..start + POLYT1_PACKED_BYTES]);
        unpack_t1(&mut t1[i], &buf);
    }

    // ── Unpack signature (c̃, z, h) ──────────────────────────────
    let ctilde: [u8; LAMBDA_BYTES] = sig_bytes[..LAMBDA_BYTES].try_into().unwrap();
    let mut sig_off = LAMBDA_BYTES;

    let mut z = [Poly::ZERO; L];
    for i in 0..L {
        let mut buf = [0u8; POLYZ_PACKED_20];
        buf.copy_from_slice(&sig_bytes[sig_off..sig_off + POLYZ_PACKED_20]);
        unpack_z_gamma1_19(&mut z[i], &buf, GAMMA1);
        sig_off += POLYZ_PACKED_20;
    }

    let hint_bytes = &sig_bytes[sig_off..sig_off + OMEGA + K];
    let mut h = [Poly::ZERO; K];
    if !unpack_hint::<K>(&mut h, hint_bytes, OMEGA) {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Validity checks ───────────────────────────────────────────
    // Check 1: ||z||_∞ < γ1 - β
    for p in &z {
        if !p.check_norm(GAMMA1 - BETA) {
            return Err(VerifyError::InvalidSignature);
        }
    }

    // Check 2: hint weight ≤ ω
    let h_weight: usize = h.iter().map(|p| p.popcount()).sum();
    if h_weight > OMEGA {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Compute μ = H(tr || prefix || message) ────────────────────
    let mut tr = [0u8; TR_BYTES];
    shake256(&[pk_bytes as &[u8]], &mut tr);
    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    // ── c = SampleInBall(c̃), NTT(c) ─────────────────────────────
    let mut c = sample_in_ball(&ctilde, TAU);
    ntt(&mut c);

    // ── Expand A and NTT-transform t1 ────────────────────────────
    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() { ntt(p); }
    }

    let mut t1_ntt = t1;
    for p in t1_ntt.iter_mut() { ntt(p); }

    // ── Compute w'_approx = Az - ct1·2^D ────────────────────────
    let mut z_ntt = z;
    for p in z_ntt.iter_mut() { ntt(p); }

    // Az
    let mut az = matrix_vector_product::<K, L>(&a_mat, &z_ntt);
    for p in az.iter_mut() { invntt_tomont(p); }

    // ct1
    let mut ct1 = [Poly::ZERO; K];
    for i in 0..K {
        ct1[i] = c.pointwise_montgomery(&t1_ntt[i]);
    }
    for p in ct1.iter_mut() { invntt_tomont(p); }

    // w'_approx = Az - ct1 * 2^D mod q, map to [0, q-1]
    let mut w_approx = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            let val = (az[i].coeffs[j] as i64)
                - ((ct1[i].coeffs[j] as i64) << D);
            let reduced = val.rem_euclid(Q as i64) as i32;
            w_approx[i].coeffs[j] = reduced;
        }
    }

    // ── w1_recovered = UseHint(h, w'_approx) ─────────────────────
    let mut w1_recovered = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            w1_recovered[i].coeffs[j] = crate::poly::use_hint_coeff(
                h[i].coeffs[j],
                w_approx[i].coeffs[j],
                GAMMA2,
            );
        }
    }

    // ── Recompute challenge hash and compare ─────────────────────
    let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
    for i in 0..K {
        let mut buf = [0u8; POLYW1_PACKED_32];
        pack_w1_gamma2_32(&mut buf, &w1_recovered[i]);
        w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES]
            .copy_from_slice(&buf);
    }

    let mut ctilde_check = [0u8; LAMBDA_BYTES];
    shake256(&[&mu as &[u8], &w1_packed], &mut ctilde_check);

    // Constant-time comparison
    use subtle::ConstantTimeEq;
    if ctilde.as_ref().ct_eq(ctilde_check.as_ref()).unwrap_u8() == 0 {
        return Err(VerifyError::Forgery);
    }

    Ok(())
}

/// Verify a PRISM-256 signature.
#[cfg(feature = "prism256")]
pub fn verify256(
    pk_bytes: &[u8; p256::PK_BYTES],
    message: &[u8],
    context: &[u8],
    sig_bytes: &[u8; p256::SIG_BYTES],
) -> Result<(), VerifyError> {
    use p256::*;

    if context.len() > 255 {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Unpack public key (ρ, t1) ─────────────────────────────────
    let rho: [u8; SEED_BYTES] = pk_bytes[..SEED_BYTES].try_into().unwrap();
    let mut t1 = [Poly::ZERO; K];
    for i in 0..K {
        let start = SEED_BYTES + i * POLYT1_PACKED_BYTES;
        let mut buf = [0u8; POLYT1_PACKED_BYTES];
        buf.copy_from_slice(&pk_bytes[start..start + POLYT1_PACKED_BYTES]);
        unpack_t1(&mut t1[i], &buf);
    }

    // ── Unpack signature (c̃, z, h) ──────────────────────────────
    let ctilde: [u8; LAMBDA_BYTES] = sig_bytes[..LAMBDA_BYTES].try_into().unwrap();
    let mut sig_off = LAMBDA_BYTES;

    let mut z = [Poly::ZERO; L];
    for i in 0..L {
        let mut buf = [0u8; POLYZ_PACKED_20];
        buf.copy_from_slice(&sig_bytes[sig_off..sig_off + POLYZ_PACKED_20]);
        unpack_z_gamma1_19(&mut z[i], &buf, GAMMA1);
        sig_off += POLYZ_PACKED_20;
    }

    let hint_bytes = &sig_bytes[sig_off..sig_off + OMEGA + K];
    let mut h = [Poly::ZERO; K];
    if !unpack_hint::<K>(&mut h, hint_bytes, OMEGA) {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Validity checks ───────────────────────────────────────────
    // Check 1: ||z||_∞ < γ1 - β
    for p in &z {
        if !p.check_norm(GAMMA1 - BETA) {
            return Err(VerifyError::InvalidSignature);
        }
    }

    // Check 2: hint weight ≤ ω
    let h_weight: usize = h.iter().map(|p| p.popcount()).sum();
    if h_weight > OMEGA {
        return Err(VerifyError::InvalidSignature);
    }

    // ── Compute μ = H(tr || prefix || message) ────────────────────
    let mut tr = [0u8; TR_BYTES];
    shake256(&[pk_bytes as &[u8]], &mut tr);
    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    // ── c = SampleInBall(c̃), NTT(c) ─────────────────────────────
    let mut c = sample_in_ball(&ctilde, TAU);
    ntt(&mut c);

    // ── Expand A and NTT-transform t1 ────────────────────────────
    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() {
        for p in row.iter_mut() { ntt(p); }
    }

    let mut t1_ntt = t1;
    for p in t1_ntt.iter_mut() { ntt(p); }

    // ── Compute w'_approx = Az - ct1·2^D ────────────────────────
    let mut z_ntt = z;
    for p in z_ntt.iter_mut() { ntt(p); }

    // Az
    let mut az = matrix_vector_product::<K, L>(&a_mat, &z_ntt);
    for p in az.iter_mut() { invntt_tomont(p); }

    // ct1
    let mut ct1 = [Poly::ZERO; K];
    for i in 0..K {
        ct1[i] = c.pointwise_montgomery(&t1_ntt[i]);
    }
    for p in ct1.iter_mut() { invntt_tomont(p); }

    // w'_approx = Az - ct1 * 2^D mod q, map to [0, q-1]
    let mut w_approx = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            let val = (az[i].coeffs[j] as i64)
                - ((ct1[i].coeffs[j] as i64) << D);
            let reduced = val.rem_euclid(Q as i64) as i32;
            w_approx[i].coeffs[j] = reduced;
        }
    }

    // ── w1_recovered = UseHint(h, w'_approx) ─────────────────────
    let mut w1_recovered = [Poly::ZERO; K];
    for i in 0..K {
        for j in 0..N {
            w1_recovered[i].coeffs[j] = crate::poly::use_hint_coeff(
                h[i].coeffs[j],
                w_approx[i].coeffs[j],
                GAMMA2,
            );
        }
    }

    // ── Recompute challenge hash and compare ─────────────────────
    let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
    for i in 0..K {
        let mut buf = [0u8; POLYW1_PACKED_32];
        pack_w1_gamma2_32(&mut buf, &w1_recovered[i]);
        w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES]
            .copy_from_slice(&buf);
    }

    let mut ctilde_check = [0u8; LAMBDA_BYTES];
    shake256(&[&mu as &[u8], &w1_packed], &mut ctilde_check);

    // Constant-time comparison
    use subtle::ConstantTimeEq;
    if ctilde.as_ref().ct_eq(ctilde_check.as_ref()).unwrap_u8() == 0 {
        return Err(VerifyError::Forgery);
    }

    Ok(())
}
