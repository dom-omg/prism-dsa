// Fuzz verify128: should never panic, only return Ok or Err
#![no_main]
use libfuzzer_sys::fuzz_target;
use prism_dsa::verify128;
use prism_dsa::params::p128::{PK_BYTES, SIG_BYTES};

fuzz_target!(|data: &[u8]| {
    if data.len() < PK_BYTES + SIG_BYTES + 2 {
        return;
    }
    let pk: [u8; PK_BYTES] = data[..PK_BYTES].try_into().unwrap();
    let sig: [u8; SIG_BYTES] = data[PK_BYTES..PK_BYTES + SIG_BYTES].try_into().unwrap();
    let msg_len = data[PK_BYTES + SIG_BYTES] as usize;
    let ctx_len = data[PK_BYTES + SIG_BYTES + 1] as usize;
    let rest = &data[PK_BYTES + SIG_BYTES + 2..];
    if rest.len() < msg_len + ctx_len {
        return;
    }
    let msg = &rest[..msg_len];
    let ctx = &rest[msg_len..msg_len + ctx_len];
    // Must not panic — any return value is acceptable
    let _ = verify128(&pk, msg, ctx, &sig);
});
