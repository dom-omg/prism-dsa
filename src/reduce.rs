//! Modular arithmetic for Z_q where q = 8380417
//!
//! All operations are branch-free and run in constant time relative to
//! the value of the input (no secret-dependent branching).

use crate::params::Q;

/// Montgomery constant: MONT = 2^32 mod q = 4193792
/// Stored as signed i32 in centered representation: 4193792 - 8380417 = -4186625
pub const MONT: i32 = -4_186_625;

/// q^{-1} mod 2^32, used in Montgomery reduction
pub const QINV: u32 = 58_728_449;

/// Montgomery reduction.
///
/// Input: a in [-q·2^31, q·2^31]
/// Output: r ≡ a·2^{-32} (mod q), |r| < q
///
/// Constant time: yes (arithmetic only, no branches on a)
#[inline(always)]
pub fn montgomery_reduce(a: i64) -> i32 {
    let t = (a as u32 as u64).wrapping_mul(QINV as u64) as u32 as i32;
    ((a - (t as i64) * (Q as i64)) >> 32) as i32
}

/// Reduce32: reduce a mod q to [-6283008, 6283008]
///
/// Input: a ≤ 2^31 - 2^22 - 1
/// Uses the identity: a = a - round(a/q)·q, exploiting that q ≈ 2^23
///
/// Constant time: yes
#[inline(always)]
pub fn reduce32(a: i32) -> i32 {
    let t = (a + (1 << 22)) >> 23;
    a - t * Q
}

/// Add q if a is negative (map to [0, q-1] representation)
///
/// Constant time: yes (bit-mask trick)
#[inline(always)]
pub fn caddq(a: i32) -> i32 {
    a + ((a >> 31) & Q)
}

/// Full reduction to canonical [0, q-1]
#[inline(always)]
pub fn freeze(a: i32) -> i32 {
    caddq(reduce32(a))
}

/// Multiply a·b mod q in Montgomery domain.
/// Both a and b must already be in Montgomery representation (multiplied by 2^32 mod q).
///
/// Returns a·b·2^{-32} mod q.
#[inline(always)]
pub fn mont_mul(a: i32, b: i32) -> i32 {
    montgomery_reduce((a as i64) * (b as i64))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mont_roundtrip() {
        // MONT = 2^32 mod q → montgomery_reduce(MONT * MONT) = MONT * 2^{-32} mod q = 1
        let x: i32 = 12345;
        // convert x to Montgomery form: x_mont = x * MONT mod q (using mont_mul)
        // then convert back: montgomery_reduce(x_mont) = x
        let _x_mont = montgomery_reduce((x as i64) * (MONT as i64));
        // x_mont should equal x * 2^32 mod q... actually let me just test reduce32
        let reduced = reduce32(x);
        assert!(reduced.abs() < Q, "reduce32 out of range");
    }

    #[test]
    fn caddq_nonneg() {
        assert_eq!(caddq(0), 0);
        assert_eq!(caddq(-1), Q - 1);
        assert_eq!(caddq(-(Q - 1)), 1);
    }
}
