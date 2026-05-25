//! Sampling operations for PRISM-DSA
//!
//! ExpandA: deterministically expand ρ into the module matrix A ∈ R_q^{K×L}
//! ExpandS: expand ρ' into short secret vectors s1, s2
//! ExpandMask: expand ρ' into uniform nonce y ∈ [-γ1, γ1]^L
//! SampleInBall: sample sparse ternary challenge c ∈ R_q with ||c||_0 = τ, coefficients in {-1,0,1}

use sha3::digest::XofReader;
use crate::hash::xof_reader;
use crate::params::{N, Q, SEED_BYTES};
use crate::poly::Poly;

#[allow(dead_code)]

/// Rejection-sample a uniform coefficient in [0, q-1] from a byte stream.
/// Uses the rej_uniform technique: draw 3 bytes, interpret as 23-bit value,
/// accept if < q. Expected draws: 1/(q/2^23) ≈ 1.11 bytes per accepted coeff.
fn rej_uniform<R: XofReader>(reader: &mut R) -> i32 {
    loop {
        let mut buf = [0u8; 3];
        reader.read(&mut buf);
        let val = (buf[0] as u32) | ((buf[1] as u32) << 8) | (((buf[2] & 0x7f) as u32) << 16);
        if val < Q as u32 {
            return val as i32;
        }
    }
}

/// ExpandA: generate the K×L module matrix from public seed ρ
///
/// A[i][j] = SHAKE256(ρ || j || i) → Poly uniform in R_q
/// Both i (row) and j (col) are encoded as single bytes.
pub fn expand_a<const K: usize, const L: usize>(rho: &[u8; SEED_BYTES]) -> [[Poly; L]; K] {
    let zero = Poly::ZERO;
    let mut mat = [[zero; L]; K];

    for i in 0..K {
        for j in 0..L {
            let mut xof = xof_reader(&[rho as &[u8], &[j as u8, i as u8]]);
            let mut p = Poly::ZERO;
            let mut count = 0;
            while count < N {
                let mut buf = [0u8; 3];
                xof.read(&mut buf);
                let val = (buf[0] as u32)
                    | ((buf[1] as u32) << 8)
                    | (((buf[2] & 0x7f) as u32) << 16);
                if val < Q as u32 {
                    p.coeffs[count] = val as i32;
                    count += 1;
                }
            }
            mat[i][j] = p;
        }
    }
    mat
}

/// ExpandS: generate short secret vectors from extended seed ρ'
/// Coefficients drawn from {-η, ..., η} using CBD-like rejection.
/// eta=2: rej_eta2; eta=4: rej_eta4
pub fn expand_s<const SIZE: usize>(rho_prime: &[u8], nonce: u16, eta: i32) -> [Poly; SIZE] {
    let zero = Poly::ZERO;
    let mut out = [zero; SIZE];

    for i in 0..SIZE {
        let nonce_bytes = (nonce + i as u16).to_le_bytes();
        let mut xof = xof_reader(&[rho_prime, &nonce_bytes]);
        let p = &mut out[i];
        let mut count = 0;
        while count < N {
            let mut buf = [0u8; 1];
            xof.read(&mut buf);
            let b = buf[0];
            let b0 = (b & 0x0F) as i32;
            let b1 = (b >> 4) as i32;

            if eta == 2 {
                if b0 < 15 {
                    let t = b0 - ((205 * b0) >> 10) * 5; // b0 mod 5
                    p.coeffs[count] = 2 - t;
                    count += 1;
                }
                if count < N && b1 < 15 {
                    let t = b1 - ((205 * b1) >> 10) * 5;
                    p.coeffs[count] = 2 - t;
                    count += 1;
                }
            } else {
                // eta=4: coefficients in {0,...,8} mapped to {-4,...,4}
                if b0 < 9 {
                    p.coeffs[count] = 4 - b0;
                    count += 1;
                }
                if count < N && b1 < 9 {
                    p.coeffs[count] = 4 - b1;
                    count += 1;
                }
            }
        }
    }
    out
}

/// ExpandMask: generate nonce y from ρ' and nonce counter
/// Coefficients in (-γ1, γ1], encoded as (γ1 - c) in the bitstream.
///
/// gamma1 = 2^17: uses 18 bits per coefficient (9/4 bytes average, 576 bytes packed)
/// gamma1 = 2^19: uses 20 bits per coefficient (640 bytes packed)
pub fn expand_mask<const L: usize>(
    rho_prime: &[u8],
    nonce: u16,
    gamma1: i32,
) -> [Poly; L] {
    let zero = Poly::ZERO;
    let mut y = [zero; L];

    let bits = if gamma1 == (1 << 17) { 18usize } else { 20usize };

    for i in 0..L {
        let kappa = (nonce + i as u16).to_le_bytes();
        let mut xof = xof_reader(&[rho_prime, &kappa]);

        let bytes_needed = N * bits / 8;
        let mut buf = vec![0u8; bytes_needed];
        xof.read(&mut buf);

        let p = &mut y[i];
        if bits == 18 {
            // 18 bits per coefficient, packed in groups of 9 bytes → 4 coefficients
            for j in 0..N / 4 {
                let b = &buf[j * 9..(j + 1) * 9];
                let z0 = gamma1 - (((b[0] as i32) | ((b[1] as i32) << 8) | ((b[2] as i32) << 16)) & 0x3FFFF);
                let z1 = gamma1 - ((((b[2] as i32) >> 2) | ((b[3] as i32) << 6) | ((b[4] as i32) << 14)) & 0x3FFFF);
                let z2 = gamma1 - ((((b[4] as i32) >> 4) | ((b[5] as i32) << 4) | ((b[6] as i32) << 12)) & 0x3FFFF);
                let z3 = gamma1 - ((((b[6] as i32) >> 6) | ((b[7] as i32) << 2) | ((b[8] as i32) << 10)) & 0x3FFFF);
                p.coeffs[4 * j] = z0;
                p.coeffs[4 * j + 1] = z1;
                p.coeffs[4 * j + 2] = z2;
                p.coeffs[4 * j + 3] = z3;
            }
        } else {
            // 20 bits per coefficient, packed in groups of 5 bytes → 2 coefficients
            for j in 0..N / 2 {
                let b = &buf[j * 5..(j + 1) * 5];
                let z0 = gamma1 - (((b[0] as i32) | ((b[1] as i32) << 8) | ((b[2] as i32) << 16)) & 0xFFFFF);
                let z1 = gamma1 - ((((b[2] as i32) >> 4) | ((b[3] as i32) << 4) | ((b[4] as i32) << 12)) & 0xFFFFF);
                p.coeffs[2 * j] = z0;
                p.coeffs[2 * j + 1] = z1;
            }
        }
    }
    y
}

/// SampleInBall: sample sparse ternary challenge polynomial
/// Output: c ∈ R_q with exactly τ non-zero coefficients, each ±1
///
/// Algorithm (Dilithium-compatible):
///   Use SHAKE256(c̃) as XOF
///   First 8 bytes → bitmask for sign bits
///   Then sample positions i in {0..255} with Fisher-Yates-like rejection
///
/// Constant time: no (position sampling has variable-time rejection)
/// This is acceptable since c̃ is public.
pub fn sample_in_ball(ctilde: &[u8], tau: usize) -> Poly {
    let mut c = Poly::ZERO;
    let mut xof = xof_reader(&[ctilde]);

    let mut sign_bytes = [0u8; 8];
    xof.read(&mut sign_bytes);
    let mut signs: u64 = u64::from_le_bytes(sign_bytes);

    // Fisher-Yates: for i in (256-τ)..256, swap c[i] with c[uniform(0..=i)]
    // c starts as [0; 256]; placing a ±1 at position j means c[j] = ±1 (displaced from i)
    for i in (N - tau)..N {
        let j = loop {
            let mut b = [0u8; 1];
            xof.read(&mut b);
            if (b[0] as usize) <= i {
                break b[0] as usize;
            }
        };
        c.coeffs[i] = c.coeffs[j];
        c.coeffs[j] = 1 - 2 * ((signs & 1) as i32);
        signs >>= 1;
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_in_ball_weight() {
        let ctilde = [42u8; 32];
        let c = sample_in_ball(&ctilde, 39);
        let weight: usize = c.coeffs.iter().filter(|&&x| x != 0).count();
        assert_eq!(weight, 39);
        for &coeff in &c.coeffs {
            assert!(coeff == -1 || coeff == 0 || coeff == 1);
        }
    }

    #[test]
    fn expand_a_deterministic() {
        let rho = [0u8; 32];
        let a1 = expand_a::<4, 4>(&rho);
        let a2 = expand_a::<4, 4>(&rho);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(a1[i][j].coeffs, a2[i][j].coeffs);
            }
        }
    }
}
