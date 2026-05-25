//! Bit-packing for PRISM-DSA keys and signatures
//!
//! All pack/unpack functions are inverse of each other:
//! unpack(pack(x)) == x for all valid x.

use crate::params::{N, D};
use crate::poly::Poly;

// ─────────────────────────────────────────────────────────
// Public key component t1: 10 bits per coefficient
// t1 coefficients ∈ [0, (q-1)/2^D] = [0, 1023]
// 256 coefficients × 10 bits = 320 bytes
// ─────────────────────────────────────────────────────────

pub const POLYT1_PACKED: usize = 320;

pub fn pack_t1(out: &mut [u8; POLYT1_PACKED], p: &Poly) {
    for i in 0..N / 4 {
        let a = [
            p.coeffs[4 * i] as u32,
            p.coeffs[4 * i + 1] as u32,
            p.coeffs[4 * i + 2] as u32,
            p.coeffs[4 * i + 3] as u32,
        ];
        out[5 * i] = (a[0]) as u8;
        out[5 * i + 1] = ((a[0] >> 8) | (a[1] << 2)) as u8;
        out[5 * i + 2] = ((a[1] >> 6) | (a[2] << 4)) as u8;
        out[5 * i + 3] = ((a[2] >> 4) | (a[3] << 6)) as u8;
        out[5 * i + 4] = (a[3] >> 2) as u8;
    }
}

pub fn unpack_t1(p: &mut Poly, buf: &[u8; POLYT1_PACKED]) {
    for i in 0..N / 4 {
        p.coeffs[4 * i] = ((buf[5 * i] as u32 | ((buf[5 * i + 1] as u32) << 8)) & 0x3FF) as i32;
        p.coeffs[4 * i + 1] = (((buf[5 * i + 1] as u32) >> 2 | ((buf[5 * i + 2] as u32) << 6)) & 0x3FF) as i32;
        p.coeffs[4 * i + 2] = (((buf[5 * i + 2] as u32) >> 4 | ((buf[5 * i + 3] as u32) << 4)) & 0x3FF) as i32;
        p.coeffs[4 * i + 3] = (((buf[5 * i + 3] as u32) >> 6 | ((buf[5 * i + 4] as u32) << 2)) & 0x3FF) as i32;
    }
}

// ─────────────────────────────────────────────────────────
// Secret key component t0: 13 bits per coefficient (signed)
// t0 coefficients ∈ (-2^{D-1}, 2^{D-1}] = (-4096, 4096]
// Stored as (2^{D-1} - t0), unsigned 13-bit value
// 256 × 13 bits = 416 bytes
// ─────────────────────────────────────────────────────────

pub const POLYT0_PACKED: usize = 416;

pub fn pack_t0(out: &mut [u8; POLYT0_PACKED], p: &Poly) {
    for i in 0..N / 8 {
        let a: [u32; 8] = core::array::from_fn(|j| {
            ((1 << (D - 1)) - p.coeffs[8 * i + j]) as u32
        });
        out[13 * i] = a[0] as u8;
        out[13 * i + 1] = ((a[0] >> 8) | (a[1] << 5)) as u8;
        out[13 * i + 2] = (a[1] >> 3) as u8;
        out[13 * i + 3] = ((a[1] >> 11) | (a[2] << 2)) as u8;
        out[13 * i + 4] = ((a[2] >> 6) | (a[3] << 7)) as u8;
        out[13 * i + 5] = (a[3] >> 1) as u8;
        out[13 * i + 6] = ((a[3] >> 9) | (a[4] << 4)) as u8;
        out[13 * i + 7] = (a[4] >> 4) as u8;
        out[13 * i + 8] = ((a[4] >> 12) | (a[5] << 1)) as u8;
        out[13 * i + 9] = ((a[5] >> 7) | (a[6] << 6)) as u8;
        out[13 * i + 10] = (a[6] >> 2) as u8;
        out[13 * i + 11] = ((a[6] >> 10) | (a[7] << 3)) as u8;
        out[13 * i + 12] = (a[7] >> 5) as u8;
    }
}

pub fn unpack_t0(p: &mut Poly, buf: &[u8; POLYT0_PACKED]) {
    for i in 0..N / 8 {
        let b = &buf[13 * i..13 * i + 13];
        let a: [u32; 8] = [
            (b[0] as u32) | (((b[1] as u32) & 0x1F) << 8),
            ((b[1] as u32) >> 5) | ((b[2] as u32) << 3) | (((b[3] as u32) & 0x03) << 11),
            ((b[3] as u32) >> 2) | (((b[4] as u32) & 0x7F) << 6),
            ((b[4] as u32) >> 7) | ((b[5] as u32) << 1) | (((b[6] as u32) & 0x0F) << 9),
            ((b[6] as u32) >> 4) | ((b[7] as u32) << 4) | (((b[8] as u32) & 0x01) << 12),
            ((b[8] as u32) >> 1) | (((b[9] as u32) & 0x3F) << 7),
            ((b[9] as u32) >> 6) | ((b[10] as u32) << 2) | (((b[11] as u32) & 0x07) << 10),
            ((b[11] as u32) >> 3) | ((b[12] as u32) << 5),
        ];
        for j in 0..8 {
            p.coeffs[8 * i + j] = (1 << (D - 1)) as i32 - (a[j] & 0x1FFF) as i32;
        }
    }
}

// ─────────────────────────────────────────────────────────
// Secret key coefficients eta=2: 3 bits per coefficient
// Values in {-2,-1,0,1,2} stored as (2-coeff) ∈ {0,1,2,3,4}
// 256 × 3 bits = 96 bytes
// ─────────────────────────────────────────────────────────

pub const POLYETA2_PACKED: usize = 96;

pub fn pack_eta2(out: &mut [u8; POLYETA2_PACKED], p: &Poly) {
    for i in 0..N / 8 {
        let a: [u8; 8] = core::array::from_fn(|j| (2 - p.coeffs[8 * i + j]) as u8);
        out[3 * i] = a[0] | (a[1] << 3) | (a[2] << 6);
        out[3 * i + 1] = (a[2] >> 2) | (a[3] << 1) | (a[4] << 4) | (a[5] << 7);
        out[3 * i + 2] = (a[5] >> 1) | (a[6] << 2) | (a[7] << 5);
    }
}

pub fn unpack_eta2(p: &mut Poly, buf: &[u8; POLYETA2_PACKED]) {
    for i in 0..N / 8 {
        let b = &buf[3 * i..3 * i + 3];
        let a = [
            (b[0] & 0x07) as i32,
            ((b[0] >> 3) & 0x07) as i32,
            ((b[0] >> 6) | ((b[1] & 0x01) << 2)) as i32,
            ((b[1] >> 1) & 0x07) as i32,
            ((b[1] >> 4) & 0x07) as i32,
            ((b[1] >> 7) | ((b[2] & 0x03) << 1)) as i32,
            ((b[2] >> 2) & 0x07) as i32,
            ((b[2] >> 5) & 0x07) as i32,
        ];
        for j in 0..8 {
            p.coeffs[8 * i + j] = 2 - a[j];
        }
    }
}

// ─────────────────────────────────────────────────────────
// Signature z component: 18 bits per coefficient (gamma1=2^17)
// z coefficients ∈ (-γ1, γ1], stored as (γ1 - z) ∈ [0, 2γ1)
// 256 × 18 bits = 576 bytes
// ─────────────────────────────────────────────────────────

pub const POLYZ_PACKED_18: usize = 576;

pub fn pack_z_gamma1_17(out: &mut [u8; POLYZ_PACKED_18], p: &Poly, gamma1: i32) {
    for i in 0..N / 4 {
        let a: [u32; 4] = core::array::from_fn(|j| (gamma1 - p.coeffs[4 * i + j]) as u32);
        out[9 * i] = a[0] as u8;
        out[9 * i + 1] = (a[0] >> 8) as u8;
        out[9 * i + 2] = ((a[0] >> 16) | (a[1] << 2)) as u8;
        out[9 * i + 3] = (a[1] >> 6) as u8;
        out[9 * i + 4] = ((a[1] >> 14) | (a[2] << 4)) as u8;
        out[9 * i + 5] = (a[2] >> 4) as u8;
        out[9 * i + 6] = ((a[2] >> 12) | (a[3] << 6)) as u8;
        out[9 * i + 7] = (a[3] >> 2) as u8;
        out[9 * i + 8] = (a[3] >> 10) as u8;
    }
}

pub fn unpack_z_gamma1_17(p: &mut Poly, buf: &[u8; POLYZ_PACKED_18], gamma1: i32) {
    for i in 0..N / 4 {
        let b = &buf[9 * i..9 * i + 9];
        let a = [
            (b[0] as u32) | ((b[1] as u32) << 8) | (((b[2] as u32) & 0x03) << 16),
            ((b[2] as u32) >> 2) | ((b[3] as u32) << 6) | (((b[4] as u32) & 0x0F) << 14),
            ((b[4] as u32) >> 4) | ((b[5] as u32) << 4) | (((b[6] as u32) & 0x3F) << 12),
            ((b[6] as u32) >> 6) | ((b[7] as u32) << 2) | ((b[8] as u32) << 10),
        ];
        for j in 0..4 {
            p.coeffs[4 * i + j] = gamma1 - (a[j] as i32);
        }
    }
}

// ─────────────────────────────────────────────────────────
// w1 packing (gamma2 = (q-1)/88): 6 bits per coefficient
// w1 ∈ [0, (q-1)/(2*gamma2) - 1] = [0, 43]
// 256 × 6 bits = 192 bytes
// ─────────────────────────────────────────────────────────

pub const POLYW1_PACKED_88: usize = 192;

pub fn pack_w1_gamma2_88(out: &mut [u8; POLYW1_PACKED_88], p: &Poly) {
    for i in 0..N / 4 {
        out[3 * i] = (p.coeffs[4 * i] | (p.coeffs[4 * i + 1] << 6)) as u8;
        out[3 * i + 1] = ((p.coeffs[4 * i + 1] >> 2) | (p.coeffs[4 * i + 2] << 4)) as u8;
        out[3 * i + 2] = ((p.coeffs[4 * i + 2] >> 4) | (p.coeffs[4 * i + 3] << 2)) as u8;
    }
}

// ─────────────────────────────────────────────────────────
// Hint vector packing: OMEGA + K bytes
// First K bytes: running offset (how many hints in poly[0..i])
// Next OMEGA bytes: positions of set bits (sorted within each polynomial)
// ─────────────────────────────────────────────────────────

pub fn pack_hint<const K: usize>(
    out: &mut Vec<u8>,
    h: &[Poly; K],
    omega: usize,
) {
    out.resize(omega + K, 0u8);
    let mut idx = 0usize; // index into the positions array

    for i in 0..K {
        for j in 0..N {
            if h[i].coeffs[j] != 0 && idx < omega {
                out[idx] = j as u8;
                idx += 1;
            }
        }
        out[omega + i] = idx as u8;
    }
}

pub fn unpack_hint<const K: usize>(
    h: &mut [Poly; K],
    buf: &[u8],
    omega: usize,
) -> bool {
    if buf.len() != omega + K {
        return false;
    }

    for i in 0..N * K {
        let poly_idx = i / N;
        let coeff_idx = i % N;
        h[poly_idx].coeffs[coeff_idx] = 0;
    }

    let mut k = 0usize;
    for i in 0..K {
        let end = buf[omega + i] as usize;
        if end < k || end > omega {
            return false;
        }
        for j in k..end {
            if j > k && buf[j] <= buf[j - 1] {
                return false; // must be strictly increasing within polynomial
            }
            h[i].coeffs[buf[j] as usize] = 1;
        }
        k = end;
    }
    true
}

// ─────────────────────────────────────────────────────────
// Secret key coefficients eta=4: 4 bits per coefficient
// Values in {-4,...,4} stored as (4-coeff) ∈ {0,...,8}
// 256 × 4 bits = 128 bytes
// ─────────────────────────────────────────────────────────

pub const POLYETA4_PACKED: usize = 128;

pub fn pack_eta4(out: &mut [u8; POLYETA4_PACKED], p: &Poly) {
    for i in 0..N / 2 {
        let a0 = (4 - p.coeffs[2 * i]) as u8;
        let a1 = (4 - p.coeffs[2 * i + 1]) as u8;
        out[i] = a0 | (a1 << 4);
    }
}

pub fn unpack_eta4(p: &mut Poly, buf: &[u8; POLYETA4_PACKED]) {
    for i in 0..N / 2 {
        p.coeffs[2 * i]     = 4 - (buf[i] & 0x0F) as i32;
        p.coeffs[2 * i + 1] = 4 - ((buf[i] >> 4) & 0x0F) as i32;
    }
}

// ─────────────────────────────────────────────────────────
// Signature z component: 20 bits per coefficient (gamma1=2^19)
// z coefficients ∈ (-γ1, γ1], stored as (γ1 - z) ∈ [0, 2γ1)
// 256 × 20 bits = 640 bytes
// ─────────────────────────────────────────────────────────

pub const POLYZ_PACKED_20: usize = 640;

pub fn pack_z_gamma1_19(out: &mut [u8; POLYZ_PACKED_20], p: &Poly, gamma1: i32) {
    for i in 0..N / 4 {
        let a: [u32; 4] = core::array::from_fn(|j| (gamma1 - p.coeffs[4 * i + j]) as u32);
        out[10 * i]     =  a[0] as u8;
        out[10 * i + 1] = (a[0] >> 8) as u8;
        out[10 * i + 2] = ((a[0] >> 16) | (a[1] << 4)) as u8;
        out[10 * i + 3] = (a[1] >> 4) as u8;
        out[10 * i + 4] = (a[1] >> 12) as u8;
        out[10 * i + 5] =  a[2] as u8;
        out[10 * i + 6] = (a[2] >> 8) as u8;
        out[10 * i + 7] = ((a[2] >> 16) | (a[3] << 4)) as u8;
        out[10 * i + 8] = (a[3] >> 4) as u8;
        out[10 * i + 9] = (a[3] >> 12) as u8;
    }
}

pub fn unpack_z_gamma1_19(p: &mut Poly, buf: &[u8; POLYZ_PACKED_20], gamma1: i32) {
    for i in 0..N / 4 {
        let b = &buf[10 * i..10 * i + 10];
        let a = [
            ((b[0] as u32) | ((b[1] as u32) << 8) | (((b[2] as u32) & 0x0F) << 16)) & 0xFFFFF,
            (((b[2] as u32) >> 4) | ((b[3] as u32) << 4) | ((b[4] as u32) << 12)) & 0xFFFFF,
            ((b[5] as u32) | ((b[6] as u32) << 8) | (((b[7] as u32) & 0x0F) << 16)) & 0xFFFFF,
            (((b[7] as u32) >> 4) | ((b[8] as u32) << 4) | ((b[9] as u32) << 12)) & 0xFFFFF,
        ];
        for j in 0..4 {
            p.coeffs[4 * i + j] = gamma1 - (a[j] as i32);
        }
    }
}

// ─────────────────────────────────────────────────────────
// w1 packing (gamma2 = (q-1)/32): 4 bits per coefficient
// w1 ∈ [0, (q-1)/(2*gamma2) - 1] = [0, 15]
// 256 × 4 bits = 128 bytes
// ─────────────────────────────────────────────────────────

pub const POLYW1_PACKED_32: usize = 128;

pub fn pack_w1_gamma2_32(out: &mut [u8; POLYW1_PACKED_32], p: &Poly) {
    for i in 0..N / 2 {
        out[i] = (p.coeffs[2 * i] as u8) | ((p.coeffs[2 * i + 1] as u8) << 4);
    }
}
