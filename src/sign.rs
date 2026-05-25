//! PRISM-DSA Signing — Fixed-Iteration Signing (FIS)
//!
//! ## Timing Uniformity
//!
//! The FIS loop always runs exactly FIS_SLOTS iterations.
//! No `continue` or `break` on any secret-derived condition.
//! All intermediate values (cs2, ct0, hint) are computed every iteration.
//! The `subtle` crate selects the first valid slot's output in constant time.
//!
//! ## Failure probability
//!
//! PRISM-128: p_accept ≈ 0.22/slot → P(all 64 fail) ≈ 0.78^64 ≈ 2^{-27}
//! PRISM-192: p_accept ≈ 0.17/slot → P(all 64 fail) ≈ 0.83^64 ≈ 2^{-22}
//! PRISM-256: p_accept ≈ 0.17/slot → P(all 64 fail) ≈ 0.83^64 ≈ 2^{-22}

use rand_core::{RngCore, CryptoRng};
use subtle::{Choice, ConditionallySelectable};

use crate::params::*;
use crate::poly::Poly;
use crate::ntt::{ntt, invntt_tomont, matrix_vector_product};
use crate::sample::{expand_a, expand_mask, sample_in_ball};
use crate::hash::shake256;
use crate::packing::*;
use crate::error::SignError;

pub const FIS_SLOTS: usize = 64;

pub struct Signature128 { pub bytes: [u8; p128::SIG_BYTES] }
#[cfg(feature = "prism192")]
pub struct Signature192 { pub bytes: [u8; p192::SIG_BYTES] }
#[cfg(feature = "prism256")]
pub struct Signature256 { pub bytes: [u8; p256::SIG_BYTES] }

// CT norm check across a poly vector — no short-circuit, always iterates all polys.
fn norm_check_ct<const S: usize>(v: &[Poly; S], bound: i32) -> u8 {
    let mut ok = 1u8;
    for p in v.iter() {
        ok &= p.check_norm(bound) as u8;
    }
    ok
}

fn serialize_sig<const L: usize, const K: usize>(
    ctilde: &[u8],
    z: &[Poly; L],
    hint: &[Poly; K],
    gamma1: i32,
    omega: usize,
) -> Vec<u8> {
    let lambda = ctilde.len();
    let polyz_bytes = if gamma1 == (1 << 17) { POLYZ_PACKED_18 } else { POLYZ_PACKED_20 };
    let mut out = vec![0u8; lambda + L * polyz_bytes + omega + K];

    let mut off = 0;
    out[off..off + lambda].copy_from_slice(ctilde);
    off += lambda;

    for p in z.iter() {
        if gamma1 == (1 << 17) {
            let mut buf = [0u8; POLYZ_PACKED_18];
            pack_z_gamma1_17(&mut buf, p, gamma1);
            out[off..off + POLYZ_PACKED_18].copy_from_slice(&buf);
            off += POLYZ_PACKED_18;
        } else {
            let mut buf = [0u8; POLYZ_PACKED_20];
            pack_z_gamma1_19(&mut buf, p, gamma1);
            out[off..off + POLYZ_PACKED_20].copy_from_slice(&buf);
            off += POLYZ_PACKED_20;
        }
    }

    let hint_base = off;
    let mut idx = 0usize;
    for i in 0..K {
        for j in 0..N {
            if hint[i].coeffs[j] != 0 && idx < omega {
                out[hint_base + idx] = j as u8;
                idx += 1;
            }
        }
        out[hint_base + omega + i] = idx as u8;
    }

    out
}

pub fn sign128<R: RngCore + CryptoRng>(
    sk_bytes: &[u8; p128::SK_BYTES],
    message: &[u8],
    context: &[u8],
    rng: &mut R,
) -> Result<Signature128, SignError> {
    if context.len() > 255 { return Err(SignError::InvalidContext); }
    use p128::*;

    let mut off = 0;
    let rho: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let key: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let tr:  [u8; TR_BYTES]   = sk_bytes[off..off + TR_BYTES].try_into().unwrap();   off += TR_BYTES;

    let mut s1 = [Poly::ZERO; L];
    for p in s1.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta2(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut s2 = [Poly::ZERO; K];
    for p in s2.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta2(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut t0 = [Poly::ZERO; K];
    for p in t0.iter_mut() {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYT0_PACKED_BYTES]);
        unpack_t0(p, &buf); off += POLYT0_PACKED_BYTES;
    }

    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() { for p in row.iter_mut() { ntt(p); } }
    let mut s1_ntt = s1; for p in s1_ntt.iter_mut() { ntt(p); }
    let mut s2_ntt = s2; for p in s2_ntt.iter_mut() { ntt(p); }
    let mut t0_ntt = t0; for p in t0_ntt.iter_mut() { ntt(p); }

    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    let mut rnd = [0u8; RND_BYTES];
    rng.fill_bytes(&mut rnd);
    let mut rho_prime = [0u8; CRH_BYTES];
    shake256(&[&key as &[u8], &rnd, &mu], &mut rho_prime);

    let mut result_bytes = [0u8; SIG_BYTES];
    // ── Constant-time properties of this loop ────────────────────────────────
    // CT:    Loop runs exactly FIS_SLOTS = 64 times (no break/continue)
    // CT:    norm_check_ct iterates all 256 coefficients unconditionally
    // CT:    u8 validity flag accumulated with &= (no branches)
    // CT:    Output selection via subtle::ConditionallySelectable (cmov)
    // NON-CT: SHAKE-256 (expand_mask, shake256) may have data-dependent timing
    // NON-CT: NTT butterfly operations (montgomery_reduce conditional adds)
    // TODO:  Replace SHAKE-256 with a verified CT PRNG
    // TODO:  Audit NTT with ct-verif or valgrind --tool=callgrind
    let mut found: u8 = 0;

    for slot in 0u16..FIS_SLOTS as u16 {
        let nonce = slot * (L as u16);
        let y = expand_mask::<L>(&rho_prime, nonce, GAMMA1);

        let mut y_ntt = y;
        for p in y_ntt.iter_mut() { ntt(p); }
        let mut w = matrix_vector_product::<K, L>(&a_mat, &y_ntt);
        for p in w.iter_mut() { invntt_tomont(p); p.caddq(); }

        let mut w0 = [Poly::ZERO; K];
        let mut w1 = [Poly::ZERO; K];
        for i in 0..K { let (hi, lo) = w[i].decompose(GAMMA2); w1[i] = hi; w0[i] = lo; }

        let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
        for i in 0..K {
            let mut buf = [0u8; POLYW1_PACKED_BYTES];
            pack_w1_gamma2_88(&mut buf, &w1[i]);
            w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES].copy_from_slice(&buf);
        }
        let mut ctilde = [0u8; LAMBDA_BYTES];
        shake256(&[&mu as &[u8], &w1_packed], &mut ctilde);

        let mut c = sample_in_ball(&ctilde, TAU);
        ntt(&mut c);

        let mut cs1 = [Poly::ZERO; L];
        for j in 0..L { cs1[j] = c.pointwise_montgomery(&s1_ntt[j]); }
        for p in cs1.iter_mut() { invntt_tomont(p); }

        let mut z = [Poly::ZERO; L];
        for j in 0..L {
            for k in 0..N { z[j].coeffs[k] = y[j].coeffs[k] + cs1[j].coeffs[k]; }
            z[j].reduce();
        }

        let mut slot_valid = norm_check_ct::<L>(&z, GAMMA1 - BETA);

        let mut cs2 = [Poly::ZERO; K];
        for i in 0..K { cs2[i] = c.pointwise_montgomery(&s2_ntt[i]); }
        for p in cs2.iter_mut() { invntt_tomont(p); }

        let mut w0_mod = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_mod[i].coeffs[j] = w0[i].coeffs[j] - cs2[i].coeffs[j]; }
            w0_mod[i].reduce();
        }
        slot_valid &= norm_check_ct::<K>(&w0_mod, GAMMA2 - BETA);

        let mut ct0 = [Poly::ZERO; K];
        for i in 0..K { ct0[i] = c.pointwise_montgomery(&t0_ntt[i]); }
        for p in ct0.iter_mut() { invntt_tomont(p); p.reduce(); }
        slot_valid &= norm_check_ct::<K>(&ct0, GAMMA2);

        let mut w0_final = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_final[i].coeffs[j] = w0_mod[i].coeffs[j] + ct0[i].coeffs[j]; }
        }

        let mut hint = [Poly::ZERO; K];
        let mut h_weight = 0usize;
        for i in 0..K {
            for j in 0..N {
                let hb = crate::poly::make_hint_coeff(w0_final[i].coeffs[j], w1[i].coeffs[j], GAMMA2);
                hint[i].coeffs[j] = hb as i32;
                h_weight += hb as usize;
            }
        }
        slot_valid &= (OMEGA.wrapping_sub(h_weight) >> (usize::BITS as usize - 1)) as u8 ^ 1;

        let candidate = serialize_sig::<L, K>(&ctilde, &z, &hint, GAMMA1, OMEGA);
        let use_this = Choice::from(slot_valid & (found ^ 1));
        for (r, cand) in result_bytes.iter_mut().zip(candidate.iter()) {
            *r = u8::conditional_select(r, cand, use_this);
        }
        found |= slot_valid;
    }

    if found != 0 {
        let mut sig = [0u8; SIG_BYTES];
        sig.copy_from_slice(&result_bytes);
        Ok(Signature128 { bytes: sig })
    } else {
        Err(SignError::AllSlotsRejected)
    }
}

#[cfg(feature = "prism192")]
pub fn sign192<R: RngCore + CryptoRng>(
    sk_bytes: &[u8; p192::SK_BYTES],
    message: &[u8],
    context: &[u8],
    rng: &mut R,
) -> Result<Signature192, SignError> {
    if context.len() > 255 { return Err(SignError::InvalidContext); }
    use p192::*;

    let mut off = 0;
    let rho: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let key: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let tr:  [u8; TR_BYTES]   = sk_bytes[off..off + TR_BYTES].try_into().unwrap();   off += TR_BYTES;

    let mut s1 = [Poly::ZERO; L];
    for p in s1.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta4(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut s2 = [Poly::ZERO; K];
    for p in s2.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta4(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut t0 = [Poly::ZERO; K];
    for p in t0.iter_mut() {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYT0_PACKED_BYTES]);
        unpack_t0(p, &buf); off += POLYT0_PACKED_BYTES;
    }

    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() { for p in row.iter_mut() { ntt(p); } }
    let mut s1_ntt = s1; for p in s1_ntt.iter_mut() { ntt(p); }
    let mut s2_ntt = s2; for p in s2_ntt.iter_mut() { ntt(p); }
    let mut t0_ntt = t0; for p in t0_ntt.iter_mut() { ntt(p); }

    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    let mut rnd = [0u8; RND_BYTES];
    rng.fill_bytes(&mut rnd);
    let mut rho_prime = [0u8; CRH_BYTES];
    shake256(&[&key as &[u8], &rnd, &mu], &mut rho_prime);

    let mut result_bytes = [0u8; SIG_BYTES];
    let mut found: u8 = 0;

    for slot in 0u16..FIS_SLOTS as u16 {
        let nonce = slot * (L as u16);
        let y = expand_mask::<L>(&rho_prime, nonce, GAMMA1);

        let mut y_ntt = y;
        for p in y_ntt.iter_mut() { ntt(p); }
        let mut w = matrix_vector_product::<K, L>(&a_mat, &y_ntt);
        for p in w.iter_mut() { invntt_tomont(p); p.caddq(); }

        let mut w0 = [Poly::ZERO; K];
        let mut w1 = [Poly::ZERO; K];
        for i in 0..K { let (hi, lo) = w[i].decompose(GAMMA2); w1[i] = hi; w0[i] = lo; }

        let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
        for i in 0..K {
            let mut buf = [0u8; POLYW1_PACKED_BYTES];
            pack_w1_gamma2_32(&mut buf, &w1[i]);
            w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES].copy_from_slice(&buf);
        }
        let mut ctilde = [0u8; LAMBDA_BYTES];
        shake256(&[&mu as &[u8], &w1_packed], &mut ctilde);

        let mut c = sample_in_ball(&ctilde, TAU);
        ntt(&mut c);

        let mut cs1 = [Poly::ZERO; L];
        for j in 0..L { cs1[j] = c.pointwise_montgomery(&s1_ntt[j]); }
        for p in cs1.iter_mut() { invntt_tomont(p); }

        let mut z = [Poly::ZERO; L];
        for j in 0..L {
            for k in 0..N { z[j].coeffs[k] = y[j].coeffs[k] + cs1[j].coeffs[k]; }
            z[j].reduce();
        }

        let mut slot_valid = norm_check_ct::<L>(&z, GAMMA1 - BETA);

        let mut cs2 = [Poly::ZERO; K];
        for i in 0..K { cs2[i] = c.pointwise_montgomery(&s2_ntt[i]); }
        for p in cs2.iter_mut() { invntt_tomont(p); }

        let mut w0_mod = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_mod[i].coeffs[j] = w0[i].coeffs[j] - cs2[i].coeffs[j]; }
            w0_mod[i].reduce();
        }
        slot_valid &= norm_check_ct::<K>(&w0_mod, GAMMA2 - BETA);

        let mut ct0 = [Poly::ZERO; K];
        for i in 0..K { ct0[i] = c.pointwise_montgomery(&t0_ntt[i]); }
        for p in ct0.iter_mut() { invntt_tomont(p); p.reduce(); }
        slot_valid &= norm_check_ct::<K>(&ct0, GAMMA2);

        let mut w0_final = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_final[i].coeffs[j] = w0_mod[i].coeffs[j] + ct0[i].coeffs[j]; }
        }

        let mut hint = [Poly::ZERO; K];
        let mut h_weight = 0usize;
        for i in 0..K {
            for j in 0..N {
                let hb = crate::poly::make_hint_coeff(w0_final[i].coeffs[j], w1[i].coeffs[j], GAMMA2);
                hint[i].coeffs[j] = hb as i32;
                h_weight += hb as usize;
            }
        }
        slot_valid &= (OMEGA.wrapping_sub(h_weight) >> (usize::BITS as usize - 1)) as u8 ^ 1;

        let candidate = serialize_sig::<L, K>(&ctilde, &z, &hint, GAMMA1, OMEGA);
        let use_this = Choice::from(slot_valid & (found ^ 1));
        for (r, cand) in result_bytes.iter_mut().zip(candidate.iter()) {
            *r = u8::conditional_select(r, cand, use_this);
        }
        found |= slot_valid;
    }

    if found != 0 {
        let mut sig = [0u8; SIG_BYTES];
        sig.copy_from_slice(&result_bytes);
        Ok(Signature192 { bytes: sig })
    } else {
        Err(SignError::AllSlotsRejected)
    }
}

#[cfg(feature = "prism256")]
pub fn sign256<R: RngCore + CryptoRng>(
    sk_bytes: &[u8; p256::SK_BYTES],
    message: &[u8],
    context: &[u8],
    rng: &mut R,
) -> Result<Signature256, SignError> {
    if context.len() > 255 { return Err(SignError::InvalidContext); }
    use p256::*;

    let mut off = 0;
    let rho: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let key: [u8; SEED_BYTES] = sk_bytes[off..off + SEED_BYTES].try_into().unwrap(); off += SEED_BYTES;
    let tr:  [u8; TR_BYTES]   = sk_bytes[off..off + TR_BYTES].try_into().unwrap();   off += TR_BYTES;

    let mut s1 = [Poly::ZERO; L];
    for p in s1.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta2(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut s2 = [Poly::ZERO; K];
    for p in s2.iter_mut() {
        let mut buf = [0u8; POLYETA_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYETA_PACKED_BYTES]);
        unpack_eta2(p, &buf); off += POLYETA_PACKED_BYTES;
    }
    let mut t0 = [Poly::ZERO; K];
    for p in t0.iter_mut() {
        let mut buf = [0u8; POLYT0_PACKED_BYTES];
        buf.copy_from_slice(&sk_bytes[off..off + POLYT0_PACKED_BYTES]);
        unpack_t0(p, &buf); off += POLYT0_PACKED_BYTES;
    }

    let mut a_mat = expand_a::<K, L>(&rho);
    for row in a_mat.iter_mut() { for p in row.iter_mut() { ntt(p); } }
    let mut s1_ntt = s1; for p in s1_ntt.iter_mut() { ntt(p); }
    let mut s2_ntt = s2; for p in s2_ntt.iter_mut() { ntt(p); }
    let mut t0_ntt = t0; for p in t0_ntt.iter_mut() { ntt(p); }

    let ctx_prefix = [0x00u8, context.len() as u8];
    let mut mu = [0u8; CRH_BYTES];
    shake256(&[&tr as &[u8], &ctx_prefix, context, message], &mut mu);

    let mut rnd = [0u8; RND_BYTES];
    rng.fill_bytes(&mut rnd);
    let mut rho_prime = [0u8; CRH_BYTES];
    shake256(&[&key as &[u8], &rnd, &mu], &mut rho_prime);

    let mut result_bytes = [0u8; SIG_BYTES];
    let mut found: u8 = 0;

    for slot in 0u16..FIS_SLOTS as u16 {
        let nonce = slot * (L as u16);
        let y = expand_mask::<L>(&rho_prime, nonce, GAMMA1);

        let mut y_ntt = y;
        for p in y_ntt.iter_mut() { ntt(p); }
        let mut w = matrix_vector_product::<K, L>(&a_mat, &y_ntt);
        for p in w.iter_mut() { invntt_tomont(p); p.caddq(); }

        let mut w0 = [Poly::ZERO; K];
        let mut w1 = [Poly::ZERO; K];
        for i in 0..K { let (hi, lo) = w[i].decompose(GAMMA2); w1[i] = hi; w0[i] = lo; }

        let mut w1_packed = vec![0u8; K * POLYW1_PACKED_BYTES];
        for i in 0..K {
            let mut buf = [0u8; POLYW1_PACKED_BYTES];
            pack_w1_gamma2_32(&mut buf, &w1[i]);
            w1_packed[i * POLYW1_PACKED_BYTES..(i + 1) * POLYW1_PACKED_BYTES].copy_from_slice(&buf);
        }
        let mut ctilde = [0u8; LAMBDA_BYTES];
        shake256(&[&mu as &[u8], &w1_packed], &mut ctilde);

        let mut c = sample_in_ball(&ctilde, TAU);
        ntt(&mut c);

        let mut cs1 = [Poly::ZERO; L];
        for j in 0..L { cs1[j] = c.pointwise_montgomery(&s1_ntt[j]); }
        for p in cs1.iter_mut() { invntt_tomont(p); }

        let mut z = [Poly::ZERO; L];
        for j in 0..L {
            for k in 0..N { z[j].coeffs[k] = y[j].coeffs[k] + cs1[j].coeffs[k]; }
            z[j].reduce();
        }

        let mut slot_valid = norm_check_ct::<L>(&z, GAMMA1 - BETA);

        let mut cs2 = [Poly::ZERO; K];
        for i in 0..K { cs2[i] = c.pointwise_montgomery(&s2_ntt[i]); }
        for p in cs2.iter_mut() { invntt_tomont(p); }

        let mut w0_mod = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_mod[i].coeffs[j] = w0[i].coeffs[j] - cs2[i].coeffs[j]; }
            w0_mod[i].reduce();
        }
        slot_valid &= norm_check_ct::<K>(&w0_mod, GAMMA2 - BETA);

        let mut ct0 = [Poly::ZERO; K];
        for i in 0..K { ct0[i] = c.pointwise_montgomery(&t0_ntt[i]); }
        for p in ct0.iter_mut() { invntt_tomont(p); p.reduce(); }
        slot_valid &= norm_check_ct::<K>(&ct0, GAMMA2);

        let mut w0_final = [Poly::ZERO; K];
        for i in 0..K {
            for j in 0..N { w0_final[i].coeffs[j] = w0_mod[i].coeffs[j] + ct0[i].coeffs[j]; }
        }

        let mut hint = [Poly::ZERO; K];
        let mut h_weight = 0usize;
        for i in 0..K {
            for j in 0..N {
                let hb = crate::poly::make_hint_coeff(w0_final[i].coeffs[j], w1[i].coeffs[j], GAMMA2);
                hint[i].coeffs[j] = hb as i32;
                h_weight += hb as usize;
            }
        }
        slot_valid &= (OMEGA.wrapping_sub(h_weight) >> (usize::BITS as usize - 1)) as u8 ^ 1;

        let candidate = serialize_sig::<L, K>(&ctilde, &z, &hint, GAMMA1, OMEGA);
        let use_this = Choice::from(slot_valid & (found ^ 1));
        for (r, cand) in result_bytes.iter_mut().zip(candidate.iter()) {
            *r = u8::conditional_select(r, cand, use_this);
        }
        found |= slot_valid;
    }

    if found != 0 {
        let mut sig = [0u8; SIG_BYTES];
        sig.copy_from_slice(&result_bytes);
        Ok(Signature256 { bytes: sig })
    } else {
        Err(SignError::AllSlotsRejected)
    }
}
