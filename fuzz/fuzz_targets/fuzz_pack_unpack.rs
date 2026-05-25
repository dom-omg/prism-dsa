// Fuzz packing roundtrips: pack(unpack(x)) == x and must not panic
#![no_main]
use libfuzzer_sys::fuzz_target;
use prism_dsa::packing::*;
use prism_dsa::poly::Poly;

fuzz_target!(|data: &[u8]| {
    // Fuzz unpack_t1 + pack_t1 roundtrip
    if data.len() >= POLYT1_PACKED {
        let buf: [u8; POLYT1_PACKED] = data[..POLYT1_PACKED].try_into().unwrap();
        let mut p = Poly::ZERO;
        unpack_t1(&mut p, &buf);
        let mut out = [0u8; POLYT1_PACKED];
        pack_t1(&mut out, &p);
        // unpack then pack must be idempotent
        assert_eq!(out, buf, "t1 pack/unpack not idempotent");
    }

    // Fuzz hint unpack (should not panic on arbitrary bytes, just return false for invalid)
    if data.len() >= 84 { // OMEGA + K = 80 + 4 = 84
        let mut h = [Poly::ZERO; 4];
        let _ = unpack_hint::<4>(&mut h, &data[..84], 80);
    }
});
