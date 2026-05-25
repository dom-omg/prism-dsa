//! SHAKE-256 based hash functions for PRISM-DSA
//!
//! PRISM-DSA uses SHAKE256 (XOF) for all hash operations:
//!   H_λ: SHAKE256 with λ output bytes (challenge hash, tr, μ)
//!   ExpandA: SHAKE256(ρ || j || i) → polynomial in Z_q
//!   ExpandMask: SHAKE256(ρ' || nonce) → masked polynomial
//!   H1: SHAKE256(μ || pack(w1)) → challenge hash c̃

use sha3::{Shake256, digest::{Update, ExtendableOutput, XofReader}};

/// One-shot SHAKE256: absorb all inputs, squeeze into output
pub fn shake256(inputs: &[&[u8]], output: &mut [u8]) {
    let mut h = Shake256::default();
    for chunk in inputs {
        h.update(chunk);
    }
    h.finalize_xof().read(output);
}

/// Build a SHAKE256 XOF reader from multiple input slices.
/// The caller squeezes bytes as needed.
pub fn xof_reader(inputs: &[&[u8]]) -> impl XofReader {
    let mut h = Shake256::default();
    for chunk in inputs {
        h.update(chunk);
    }
    h.finalize_xof()
}
