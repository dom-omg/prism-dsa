// Fuzz sign+verify roundtrip: any valid sk must produce a verifiable sig
#![no_main]
use libfuzzer_sys::fuzz_target;
use prism_dsa::{keygen128_from_seed, sign128, verify128};
use rand_core::{RngCore, CryptoRng};

struct DeterministicRng<'a> {
    data: &'a [u8],
    pos: usize,
}
impl<'a> RngCore for DeterministicRng<'a> {
    fn next_u32(&mut self) -> u32 { 0 }
    fn next_u64(&mut self) -> u64 { 0 }
    fn fill_bytes(&mut self, dest: &mut [u8]) {
        for b in dest.iter_mut() {
            *b = if self.pos < self.data.len() {
                let v = self.data[self.pos]; self.pos += 1; v
            } else { 0 };
        }
    }
    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest); Ok(())
    }
}
impl<'a> CryptoRng for DeterministicRng<'a> {}

fuzz_target!(|data: &[u8]| {
    if data.len() < 34 { return; }
    let seed: [u8; 32] = data[..32].try_into().unwrap();
    let msg_len = data[32] as usize;
    let ctx_len = data[33] as usize;
    if data.len() < 34 + msg_len + ctx_len { return; }
    let msg = &data[34..34 + msg_len];
    let ctx = &data[34 + msg_len..34 + msg_len + ctx_len];
    if ctx_len > 255 { return; }

    let (pk, sk) = match keygen128_from_seed(&seed) {
        Ok(kp) => kp,
        Err(_) => return,
    };
    let mut rng = DeterministicRng { data: &data[34 + msg_len + ctx_len..], pos: 0 };
    if let Ok(sig) = sign128(&sk.bytes, msg, ctx, &mut rng) {
        // Every produced signature MUST verify
        assert!(verify128(&pk.bytes, msg, ctx, &sig.bytes).is_ok(),
            "sign produced unverifiable signature!");
    }
});
