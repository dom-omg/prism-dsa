//! Polynomial operations over R_q = Z_q[X]/(X^N + 1)

use crate::params::{N, Q, D};
use crate::reduce::{reduce32, caddq, montgomery_reduce, MONT};

#[derive(Clone, Copy)]
pub struct Poly {
    pub coeffs: [i32; N],
}

impl Poly {
    pub const ZERO: Self = Poly { coeffs: [0i32; N] };

    /// Reduce all coefficients to [-6283008, 6283008]
    pub fn reduce(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = reduce32(*c);
        }
    }

    /// Map all coefficients to [0, q-1]
    pub fn freeze(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = reduce32(*c);
            *c = caddq(*c);
        }
    }

    /// Add q to any negative coefficient (CT, maps negative to positive representative)
    pub fn caddq(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = caddq(*c);
        }
    }

    /// a + b mod q
    pub fn add(&self, other: &Poly) -> Poly {
        let mut r = Poly::ZERO;
        for i in 0..N {
            r.coeffs[i] = self.coeffs[i] + other.coeffs[i];
        }
        r
    }

    /// a - b mod q
    pub fn sub(&self, other: &Poly) -> Poly {
        let mut r = Poly::ZERO;
        for i in 0..N {
            r.coeffs[i] = self.coeffs[i] - other.coeffs[i];
        }
        r
    }

    /// Shift left by D bits: a[i] <<= D
    pub fn shiftl(&self) -> Poly {
        let mut r = *self;
        for c in r.coeffs.iter_mut() {
            *c <<= D;
        }
        r
    }

    /// Power2Round: decompose into (a1, a0) where a = a1·2^D + a0
    /// a0 ∈ (-2^{D-1}, 2^{D-1}], a1 = (a - a0) / 2^D
    pub fn power2round(&self) -> (Poly, Poly) {
        let mut a1 = Poly::ZERO;
        let mut a0 = Poly::ZERO;
        for i in 0..N {
            let (r1, r0) = power2round_coeff(self.coeffs[i]);
            a1.coeffs[i] = r1;
            a0.coeffs[i] = r0;
        }
        (a1, a0)
    }

    /// Decompose: split into (a1, a0) where a = a1·α + a0,
    /// α = 2·γ2, with special case at the top
    pub fn decompose(&self, gamma2: i32) -> (Poly, Poly) {
        let mut a1 = Poly::ZERO;
        let mut a0 = Poly::ZERO;
        for i in 0..N {
            let (r1, r0) = decompose_coeff(self.coeffs[i], gamma2);
            a1.coeffs[i] = r1;
            a0.coeffs[i] = r0;
        }
        (a1, a0)
    }

    /// HighBits via Decompose
    pub fn high_bits(&self, gamma2: i32) -> Poly {
        self.decompose(gamma2).0
    }

    /// LowBits via Decompose
    pub fn low_bits(&self, gamma2: i32) -> Poly {
        self.decompose(gamma2).1
    }

    /// ||a||_∞: maximum absolute coefficient value
    pub fn norm_inf(&self) -> i32 {
        self.coeffs.iter().map(|&c| c.abs()).max().unwrap_or(0)
    }

    /// Check ||a||_∞ < bound (constant-time: no early exit)
    ///
    /// Returns true iff all coefficients have |c| < bound.
    /// Runs in constant time — iterates all N coefficients.
    pub fn check_norm(&self, bound: i32) -> bool {
        let mut ok = true;
        for &c in &self.coeffs {
            let r = reduce32(c);
            // reduce32 maps to approximately [-q/2, q/2]; take absolute value
            ok &= r.abs() < bound;
        }
        ok
    }

    /// Make hint: h[i] = 1 iff HighBits(r+z)[i] ≠ HighBits(r)[i]
    pub fn make_hint(z: &Poly, r: &Poly, gamma2: i32) -> Poly {
        let mut h = Poly::ZERO;
        for i in 0..N {
            h.coeffs[i] = make_hint_coeff(z.coeffs[i], r.coeffs[i], gamma2) as i32;
        }
        h
    }

    /// Apply hints to recover HighBits(r)
    pub fn use_hint(&self, h: &Poly, gamma2: i32) -> Poly {
        let mut out = Poly::ZERO;
        for i in 0..N {
            out.coeffs[i] = use_hint_coeff(h.coeffs[i], self.coeffs[i], gamma2);
        }
        out
    }

    /// Number of non-zero entries (for hint weight check)
    pub fn popcount(&self) -> usize {
        self.coeffs.iter().filter(|&&c| c != 0).count()
    }

    /// Point-wise multiplication in Montgomery domain.
    /// Both self and other must be in NTT (frequency) domain.
    /// Each coefficient must be < 2^{31}/q to avoid overflow.
    pub fn pointwise_montgomery(&self, other: &Poly) -> Poly {
        let mut r = Poly::ZERO;
        for i in 0..N {
            r.coeffs[i] = montgomery_reduce((self.coeffs[i] as i64) * (other.coeffs[i] as i64));
        }
        r
    }

    /// Convert to Montgomery domain: a[i] ← a[i] * MONT mod q
    pub fn to_mont(&mut self) {
        for c in self.coeffs.iter_mut() {
            *c = montgomery_reduce((*c as i64) * (MONT as i64));
        }
    }
}

impl Default for Poly {
    fn default() -> Self {
        Poly::ZERO
    }
}

// ─────────────────────────────────────────────────────────
// Scalar coefficient operations
// ─────────────────────────────────────────────────────────

pub fn power2round_coeff(a: i32) -> (i32, i32) {
    // Exact match of CRYSTALS-Dilithium/FIPS-204 rounding:
    // a1 = ceil(a / 2^D) using round-up for tie at 2^{D-1}
    let a1 = (a + (1 << (D - 1)) - 1) >> D;
    let a0 = a - (a1 << D);
    (a1, a0)
}

pub fn decompose_coeff(a: i32, gamma2: i32) -> (i32, i32) {
    // Normalize to [0, q-1]
    let a = caddq(reduce32(caddq(a)));

    let alpha = 2 * gamma2;
    let mut a0 = a % alpha;

    // Center a0 in (-gamma2, gamma2] = (-alpha/2, alpha/2]:
    // subtract alpha when a0 > gamma2 (i.e., a0 > alpha/2)
    a0 -= ((gamma2 - a0) >> 31) & alpha;

    let a1 = (a - a0) / alpha;
    let top = (Q - 1) / alpha;

    // Constant-time equality check: a1 == top
    // ((a1-top)|(top-a1)) has bit 31 set iff a1 != top (one side is negative)
    let ne = ((a1 - top) | (top - a1)).wrapping_shr(31) & 1; // 1 if a1 != top
    let is_top = 1 - ne;                                       // 1 if a1 == top

    // FIPS 204: if a1 == top → a1 = 0, a0 = a0 - 1
    let a1_final = a1 * ne;      // a1 if not top, 0 if top
    let a0_final = a0 - is_top;  // a0 if not top, a0-1 if top

    (a1_final, a0_final)
}

/// MakeHint in Dilithium reference format:
/// a0 = LowBits of the expression (w0 - cs2 + ct0), a1 = HighBits(w) = w1
/// Returns 1 iff the low bits overflow into the high bits
pub fn make_hint_coeff(a0: i32, a1: i32, gamma2: i32) -> bool {
    a0 > gamma2 || a0 < -gamma2 || (a0 == -gamma2 && a1 != 0)
}

/// UseHint in Dilithium reference format:
/// a = element in [0, q-1], hint = 0 or 1
/// Returns corrected HighBits
pub fn use_hint_coeff(hint: i32, a: i32, gamma2: i32) -> i32 {
    let (a1, a0) = decompose_coeff(a, gamma2);
    if hint == 0 {
        return a1;
    }
    let m = (Q - 1) / (2 * gamma2); // = 44 for gamma2=(q-1)/88
    if a0 > 0 {
        if a1 == m - 1 { 0 } else { a1 + 1 }
    } else {
        if a1 == 0 { m - 1 } else { a1 - 1 }
    }
}

// ─────────────────────────────────────────────────────────
// Module vectors: PolyVecL (length L) and PolyVecK (length K)
// ─────────────────────────────────────────────────────────

/// Vector of polynomials, generic length
#[derive(Clone)]
pub struct PolyVec<const SIZE: usize> {
    pub polys: [Poly; SIZE],
}

impl<const SIZE: usize> PolyVec<SIZE> {
    pub fn zero() -> Self {
        PolyVec {
            polys: [Poly::ZERO; SIZE],
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..SIZE {
            r.polys[i] = self.polys[i].add(&other.polys[i]);
        }
        r
    }

    pub fn sub(&self, other: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..SIZE {
            r.polys[i] = self.polys[i].sub(&other.polys[i]);
        }
        r
    }

    pub fn reduce(&mut self) {
        for p in self.polys.iter_mut() {
            p.reduce();
        }
    }

    pub fn caddq(&mut self) {
        for p in self.polys.iter_mut() {
            p.caddq();
        }
    }

    pub fn freeze(&mut self) {
        for p in self.polys.iter_mut() {
            p.freeze();
        }
    }

    pub fn norm_inf(&self) -> i32 {
        self.polys.iter().map(|p| p.norm_inf()).max().unwrap_or(0)
    }

    /// Constant-time norm check: true iff all ||p[i]||_∞ < bound
    pub fn check_norm(&self, bound: i32) -> bool {
        let mut ok = true;
        for p in &self.polys {
            ok &= p.check_norm(bound);
        }
        ok
    }

    pub fn shiftl(&self) -> Self {
        let mut r = Self::zero();
        for i in 0..SIZE {
            r.polys[i] = self.polys[i].shiftl();
        }
        r
    }

    pub fn to_mont(&mut self) {
        for p in self.polys.iter_mut() {
            p.to_mont();
        }
    }

    pub fn power2round(&self) -> (Self, Self) {
        let mut r1 = Self::zero();
        let mut r0 = Self::zero();
        for i in 0..SIZE {
            let (a1, a0) = self.polys[i].power2round();
            r1.polys[i] = a1;
            r0.polys[i] = a0;
        }
        (r1, r0)
    }

    pub fn decompose(&self, gamma2: i32) -> (Self, Self) {
        let mut r1 = Self::zero();
        let mut r0 = Self::zero();
        for i in 0..SIZE {
            let (a1, a0) = self.polys[i].decompose(gamma2);
            r1.polys[i] = a1;
            r0.polys[i] = a0;
        }
        (r1, r0)
    }

    pub fn high_bits(&self, gamma2: i32) -> Self {
        self.decompose(gamma2).0
    }

    pub fn low_bits(&self, gamma2: i32) -> Self {
        self.decompose(gamma2).1
    }

    pub fn make_hint(z: &Self, r: &Self, gamma2: i32) -> Self {
        let mut h = Self::zero();
        for i in 0..SIZE {
            h.polys[i] = Poly::make_hint(&z.polys[i], &r.polys[i], gamma2);
        }
        h
    }

    pub fn use_hint(&self, h: &Self, gamma2: i32) -> Self {
        let mut out = Self::zero();
        for i in 0..SIZE {
            out.polys[i] = self.polys[i].use_hint(&h.polys[i], gamma2);
        }
        out
    }

    /// Total hint weight: sum of popcount across all polynomials
    pub fn hint_weight(&self) -> usize {
        self.polys.iter().map(|p| p.popcount()).sum()
    }

    /// Scale each polynomial by a scalar polynomial c (in NTT domain)
    pub fn scale_by(&self, c: &Poly) -> Self {
        let mut r = Self::zero();
        for i in 0..SIZE {
            r.polys[i] = c.pointwise_montgomery(&self.polys[i]);
        }
        r
    }
}
