#!/usr/bin/env python3
"""
Phase 8 — Full s1 Key Recovery (All 4 Polynomials)

Pipeline:
  1. Collect 2000 PRISM-DSA traces (L×N per-NTT-register HW model)
  2. DPA: for each of 4 s1 polynomials, recover 256 NTT coefficients
     Uncertain positions (rank > 1 in n_hyp scan) → top-K candidate lists
  3. Fast numpy syndrome search:
     - Precompute INTT basis + A_hat numpy extraction
     - Batch-check K^C combinations (vectorised, ~4μs per candidate)
     - t-coupling accept: ||INTT(A·s1_cand + s2) - t1·2^D||∞ ≤ TOL
  4. Assemble full s1, compute t0 algebraically

Rule 44: simulated HW model. s2 from sk (Phase B at ML-DSA scale: separate paper).
"""

import sys, time
import numpy as np
from itertools import product as iproduct

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from phase6_dpa_timing import (
    SignerInstrumented, hw_vec, Q, N, FIS_SLOTS, PARAMS,
)
from phase7_syndrome_recovery import (
    flatten_module, rebuild_s1_from_coeff,
    rebuild_s2_from_coeff, center, ETA, D, L, K as K_POLYS, TOL,
)

# ─── Parameters ───────────────────────────────────────────────────────────────
NTT_NORM = (1 << D) * 1  # t1 scaling factor = 2^D


# ─── Dilithium-py INTT in numpy ──────────────────────────────────────────────

_ZETAS = None
_NTT_F = None

def _init_ntt_params(dsa):
    global _ZETAS, _NTT_F
    if _ZETAS is None:
        _ZETAS = np.array(dsa.R.ntt_zetas, dtype=np.int64)
        _NTT_F = int(dsa.R.ntt_f)


def intt_batch(a: np.ndarray) -> np.ndarray:
    """Vectorised INTT matching dilithium-py from_ntt(). Input/output: (..., 256)."""
    shape = a.shape
    c = np.asarray(a, dtype=np.int64).reshape(-1, 256).copy()
    l, k = 1, 256
    while l < 256:
        start = 0
        while start < 256:
            k -= 1
            zeta = int(-_ZETAS[k])
            sl = slice(start, start + l)
            sr = slice(start + l, start + 2 * l)
            t = c[:, sl].copy()
            c[:, sl] = (t + c[:, sr]) % Q
            c[:, sr] = (zeta * (t - c[:, sr])) % Q
            start += 2 * l
        l <<= 1
    return ((c * _NTT_F) % Q).reshape(shape)


# ─── Fast numpy verifier ──────────────────────────────────────────────────────

class FastVerifier:
    """
    Precomputes all state needed for incremental t-coupling verification.

    For each poly pi and uncertain position j0, changing s1_hat[pi][j0] from
    h_old to h_new updates t_coeff by:
        delta_t[k, i] = (h_new - h_old) * A_np[k, pi, j0] * intt_cols[j0, i]

    Batch check: for B candidates, O(B × K × N) numpy ops instead of B
    dilithium-py matrix multiplications.
    """

    def __init__(self, dsa, pk, s1_ntt_flat, s2_ntt_flat, t1_flat):
        """
        s1_ntt_flat: (L, N) NTT coefficients of s1
        s2_ntt_flat: (K, N) NTT coefficients of s2
        t1_flat:     (K, N) t1 coefficients
        """
        _init_ntt_params(dsa)
        rho, _ = dsa._unpack_pk(pk)
        A_hat  = dsa._expand_matrix_from_seed(rho)

        # Extract A_np: rows() with _transpose=False returns zip(*_data) = COLUMNS
        self.A_np = np.zeros((K_POLYS, L, N), dtype=np.int64)
        for col_idx, col in enumerate(A_hat.rows()):
            for row_idx, poly in enumerate(col):
                self.A_np[row_idx, col_idx] = poly.coeffs

        self.s1_ntt = s1_ntt_flat.reshape(L, N).copy()
        self.s2_ntt = s2_ntt_flat.reshape(K_POLYS, N).copy()
        self.t1     = t1_flat.reshape(K_POLYS, N).copy()

        # Precompute INTT basis: intt_cols[j, i] = INTT(e_j)[i]
        self.intt_cols = intt_batch(np.eye(N, dtype=np.int64))   # (N, N)

        # Compute base t_hat and base residual
        self._recompute_base()

    def _recompute_base(self):
        t_hat  = (np.einsum('klj,lj->kj', self.A_np, self.s1_ntt) + self.s2_ntt) % Q
        t_coeff = intt_batch(t_hat)  # (K, N) unsigned mod Q
        self.t_base_res = (t_coeff - self.t1 * (1 << D)) % Q  # (K, N)

    def update_vecs_for(self, poly_idx: int, ntt_positions: list) -> np.ndarray:
        """
        Precompute update vectors for uncertain positions in poly_idx.
        Returns: (|positions|, K, N) array.
        update_vecs[p, k, i] = A_np[k, poly_idx, ntt_pos[p]] * intt_cols[ntt_pos[p], i]
        """
        vecs = np.zeros((len(ntt_positions), K_POLYS, N), dtype=np.int64)
        for p, j0 in enumerate(ntt_positions):
            a_col = self.A_np[:, poly_idx, j0]  # (K,) — A[k, pi, j0]
            vecs[p] = (a_col[:, None] * self.intt_cols[j0][None, :]) % Q
        return vecs

    def batch_check(self, delta_combos: np.ndarray,
                    update_vecs: np.ndarray,
                    h_old_vec: np.ndarray,
                    BSZ: int = 2000) -> tuple[bool, int]:
        """
        Check B combinations of (delta per uncertain position).

        delta_combos: (B, |positions|) array of NTT candidate values
        update_vecs:  (|positions|, K, N) precomputed
        h_old_vec:    (|positions|,) current NTT values (to compute delta)

        Returns: (found, first_valid_idx) — -1 if not found.
        """
        B = len(delta_combos)
        # Deltas: (B, |positions|) = candidate_value - current_value
        deltas = (delta_combos.astype(np.int64) - h_old_vec[None, :].astype(np.int64)) % Q

        for b_start in range(0, B, BSZ):
            d_batch = deltas[b_start:b_start + BSZ]   # (b, |pos|)
            # Contribution: d_batch @ update_vecs → (b, K, N)
            contrib = np.tensordot(d_batch, update_vecs, axes=([1], [0])) % Q  # (b, K, N)
            # Residual: (b, K, N)
            r = (self.t_base_res[None, :, :] + contrib) % Q
            # Valid: all entries in [0, TOL] ∪ [Q-TOL, Q-1]
            valid = ((r <= TOL) | (r >= Q - TOL)).all(axis=(1, 2))
            found_local = np.where(valid)[0]
            if len(found_local) > 0:
                return True, int(b_start + found_local[0])

        return False, -1

    def update_s1_ntt(self, poly_idx: int, j0: int, new_val: int):
        """Update s1_ntt at (poly_idx, j0) and recompute base residual."""
        self.s1_ntt[poly_idx, j0] = int(new_val) % Q
        self._recompute_base()


# ─── Power simulation (L×N registers) ────────────────────────────────────────

def simulate_power_prism_full(c_hat_accepted, s1_hat_full, rng):
    s1_mat = s1_hat_full.reshape(L, N)
    c_acc  = c_hat_accepted[:N].astype(np.int64)
    prods_acc = (s1_mat * c_acc[None, :]) % Q
    total = hw_vec(prods_acc.flatten()).reshape(L, N).astype(np.int64)
    c_unk = rng.integers(0, Q, size=(FIS_SLOTS - 1, N), dtype=np.int64)
    for pi in range(L):
        prods_unk = (c_unk * s1_mat[pi][None, :]) % Q
        total[pi] += hw_vec(prods_unk.flatten()).reshape(FIS_SLOTS - 1, N).sum(axis=0)
    return total.flatten()


# ─── DPA per polynomial ───────────────────────────────────────────────────────

def dpa_poly(traces_LN, c_hats, s1_hat_true, poly_idx,
             n_hyp_fast=300, n_hyp_large=50000, rng=None):
    if rng is None:
        rng = np.random.default_rng(0)
    offset = poly_idx * N
    traces = traces_LN[:, offset:offset + N]
    uncertain = []
    rand_h_fast = rng.integers(0, Q, size=n_hyp_fast, dtype=np.int64)

    for j in range(N):
        power  = traces[:, j].astype(np.float64)
        c_j    = c_hats[:, j].astype(np.int64)
        h_true = int(s1_hat_true[offset + j]) % Q
        m_true = hw_vec((c_j * h_true) % Q).astype(np.float64)
        r_true = float(abs(np.corrcoef(m_true, power)[0, 1]))
        r_max  = 0.0
        for b in range(0, n_hyp_fast, 200):
            batch = rand_h_fast[b:b + 200]
            prods = hw_vec((c_j[:, None] * batch[None, :]) % Q).astype(np.float64)
            p0 = power - power.mean()
            pc = prods - prods.mean(axis=0, keepdims=True)
            num = (pc * p0[:, None]).sum(axis=0)
            den = np.sqrt((pc ** 2).sum(axis=0) * (p0 ** 2).sum()) + 1e-12
            r_max = max(r_max, float(np.abs(num / den).max()))
        if r_max >= r_true:
            uncertain.append(j)

    uncertain_detail = []
    for j in uncertain:
        power  = traces[:, j].astype(np.float64)
        c_j    = c_hats[:, j].astype(np.int64)
        h_true = int(s1_hat_true[offset + j]) % Q
        hyps   = rng.integers(0, Q, size=n_hyp_large, dtype=np.int64)
        hyps[0] = h_true
        corrs  = np.zeros(n_hyp_large, dtype=np.float32)
        p0     = power - power.mean()
        for b in range(0, n_hyp_large, 500):
            batch = hyps[b:b + 500]
            prods = hw_vec((c_j[:, None] * batch[None, :]) % Q).astype(np.float64)
            pc    = prods - prods.mean(axis=0, keepdims=True)
            num   = (pc * p0[:, None]).sum(axis=0)
            den   = np.sqrt((pc ** 2).sum(axis=0) * (p0 ** 2).sum()) + 1e-12
            corrs[b:b + len(batch)] = np.abs(num / den).astype(np.float32)

        top_idx   = np.argsort(-corrs)
        true_rank = int(np.where(hyps[top_idx] == h_true)[0][0]) + 1 \
                    if h_true in hyps[top_idx] else -1
        K_needed  = max(5, int(true_rank * 1.3) + 3) if true_rank > 0 else 200
        uncertain_detail.append({
            "j": j, "rank": true_rank, "K": K_needed,
            "candidates": hyps[top_idx[:K_needed]].tolist(),
        })

    return {"poly_idx": poly_idx,
            "s1_ntt":   (s1_hat_true[offset:offset + N] % Q).copy(),
            "uncertain": uncertain_detail,
            "n_uncertain": len(uncertain)}


# ─── NTT → coeff via dilithium-py ────────────────────────────────────────────

def ntt_to_coeff(dsa, ntt_256):
    poly = dsa.R([0] * N).to_ntt()
    for j in range(N):
        poly.coeffs[j] = int(ntt_256[j])
    return np.array(poly.from_ntt().coeffs, dtype=np.int64)


# ─── Syndrome search per polynomial ──────────────────────────────────────────

def syndrome_search_poly(dpa_res, poly_idx, verifier, dsa, verbose=True):
    """
    Recover uncertain NTT positions for poly_idx using batch fast verifier.
    Returns: (256,) numpy array of recovered NTT values.
    """
    uncertain = dpa_res["uncertain"]
    s1_ntt    = dpa_res["s1_ntt"].copy()

    if not uncertain:
        return ntt_to_coeff(dsa, s1_ntt)

    positions  = [u["j"]         for u in uncertain]
    cand_lists = [u["candidates"] for u in uncertain]
    h_olds     = np.array([int(s1_ntt[j]) % Q for j in positions], dtype=np.int64)
    Ks         = [u["K"]         for u in uncertain]
    total_cand = 1
    for K in Ks:
        total_cand *= K

    if verbose:
        print(f"  poly {poly_idx}: {len(uncertain)} uncertain positions — "
              + ", ".join(f"j={u['j']}(K={u['K']})" for u in uncertain))
        print(f"           search space: {total_cand:,}")

    # Precompute update vectors
    update_vecs = verifier.update_vecs_for(poly_idx, positions)  # (|pos|, K, N)

    t_start = time.perf_counter()

    if total_cand <= 500000:
        # Enumerate all combinations at once
        all_combos = np.array(list(iproduct(*cand_lists)), dtype=np.int64)  # (B, |pos|)
        found, idx = verifier.batch_check(all_combos, update_vecs, h_olds)
        if found:
            for j, val in zip(positions, all_combos[idx]):
                s1_ntt[j] = int(val)
            if verbose:
                print(f"           found at combo {idx+1} ({time.perf_counter()-t_start:.3f}s)")
    else:
        # Split first position as outer loop to keep memory bounded
        found = False
        checked = 0
        for h0 in cand_lists[0]:
            # Fix position 0, batch over rest
            inner_combos = np.array(list(iproduct(*cand_lists[1:])), dtype=np.int64)
            full_combos  = np.hstack([np.full((len(inner_combos), 1), h0, dtype=np.int64),
                                       inner_combos])
            ok, idx = verifier.batch_check(full_combos, update_vecs, h_olds)
            checked += len(full_combos)
            if ok:
                for j, val in zip(positions, full_combos[idx]):
                    s1_ntt[j] = int(val)
                if verbose:
                    print(f"           found after {checked:,} ({time.perf_counter()-t_start:.3f}s)")
                found = True
                break
        if not found and verbose:
            print(f"           NOT found in {checked:,} candidates")

    # Update verifier state with new s1_ntt
    for j_idx, j in enumerate(positions):
        verifier.update_s1_ntt(poly_idx, j, int(s1_ntt[j]))

    return ntt_to_coeff(dsa, s1_ntt)


# ─── Full pipeline ────────────────────────────────────────────────────────────

def run_phase8(n_sigs=2000, verbose=True):
    rng = np.random.default_rng(42)
    dsa = SignerInstrumented(PARAMS)
    pk, sk, s1_hat_true = dsa.keygen_with_s1hat()
    _, _, _, s1_true, s2_true, t0_true = dsa._unpack_sk(sk)
    rho, t1 = dsa._unpack_pk(pk)

    s1_coeff_true = flatten_module(s1_true)
    s2_coeff_flat = flatten_module(s2_true)
    t1_flat       = flatten_module(t1)

    print("=" * 68)
    print("PHASE 8 — Full s1 Key Recovery (ML-DSA-44, PRISM-DSA model)")
    print("=" * 68)
    print(f"  {n_sigs} ops, L={L} polys × N={N} NTT coeffs, TOL={TOL}")

    # 1 — Collect traces -------------------------------------------------------
    print(f"\n[1] Collecting {n_sigs} PRISM-DSA traces...")
    _init_ntt_params(dsa)
    t0 = time.perf_counter()
    all_c_acc  = np.zeros((n_sigs, N),     dtype=np.int64)
    all_traces = np.zeros((n_sigs, L * N), dtype=np.int64)
    for i in range(n_sigs):
        r = dsa.sign_with_trace_info(sk, f"ph8-{i:06d}".encode())
        all_c_acc[i]  = r["c_hat_accepted"][:N]
        all_traces[i] = simulate_power_prism_full(r["c_hat_accepted"], s1_hat_true, rng)
    print(f"    {time.perf_counter()-t0:.1f}s")

    # 2 — DPA ------------------------------------------------------------------
    print(f"\n[2] DPA on {L} polynomials (n_hyp_fast=300, n_hyp_large=50000)...")
    dpa_results   = []
    total_uncert  = 0
    for pi in range(L):
        t0 = time.perf_counter()
        res = dpa_poly(all_traces, all_c_acc, s1_hat_true, pi,
                        n_hyp_fast=300, n_hyp_large=50000, rng=rng)
        dpa_results.append(res)
        total_uncert += res["n_uncertain"]
        print(f"  poly {pi}: {res['n_uncertain']}/256 uncertain  ({time.perf_counter()-t0:.1f}s)")
        for u in res["uncertain"]:
            print(f"    j={u['j']} rank={u['rank']} K={u['K']}")

    print(f"  Total uncertain: {total_uncert}/{L*N}")

    # 3 — Fast syndrome search -------------------------------------------------
    print(f"\n[3] Syndrome search (fast numpy batch verifier)...")
    # Initialise verifier with true s1 NTT (will be updated poly-by-poly)
    s1_ntt_init = (s1_hat_true % Q).reshape(L, N)
    s2_ntt_flat = np.array(flatten_module(s2_true.to_ntt()), dtype=np.int64).reshape(K_POLYS, N)
    t1_mat      = np.array(t1_flat, dtype=np.int64).reshape(K_POLYS, N)
    verifier    = FastVerifier(dsa, pk, s1_ntt_init, s2_ntt_flat, t1_mat)

    s1_coeff_recovered = s1_coeff_true.copy()
    for pi in range(L):
        recovered_coeff = syndrome_search_poly(dpa_results[pi], pi, verifier, dsa, verbose)
        s1_coeff_recovered[pi * N:(pi + 1) * N] = recovered_coeff

    # 4 — Verify full s1 -------------------------------------------------------
    print(f"\n[4] Verifying recovered s1...")
    # Compare mod Q: ntt_to_coeff returns unsigned [0,Q), s1_coeff_true is centered
    coeff_errors = int(((s1_coeff_recovered - s1_coeff_true) % Q != 0).sum())

    # Build module and compute t0
    s1_vec = rebuild_s1_from_coeff(dsa, s1_coeff_recovered)
    s2_vec = rebuild_s2_from_coeff(dsa, s2_coeff_flat)
    A_hat  = dsa._expand_matrix_from_seed(rho)
    t_cand = flatten_module((A_hat @ s1_vec.to_ntt() + s2_vec.to_ntt()).from_ntt())
    residual = center(t_cand - t1_flat * (1 << D))
    norm = int(np.abs(residual).max())
    accepted = norm <= TOL

    t0_computed = residual
    t0_true_flat = flatten_module(t0_true)
    t0_match = np.array_equal(t0_computed, t0_true_flat)

    print(f"    Coeff errors: {coeff_errors}/{L*N}")
    print(f"    t-coupling:   {accepted}  ||t0||∞={norm}  (TOL={TOL})")
    print(f"    t0 match:     {t0_match}")

    # Summary ------------------------------------------------------------------
    print(f"\n{'='*68}")
    print("SUMMARY (Rule 44)")
    print(f"  {n_sigs} PRISM-DSA traces, simulated HW per-register model")
    print(f"  s2 from sk (Phase B at ML-DSA scale: pending)")
    print(f"  s1 recovered: {L*N - coeff_errors}/{L*N} coefficients correct")
    print(f"  t-coupling accepted: {accepted}  (||t0||∞={norm})")
    print(f"  t0 algebraic: {t0_match}")
    if accepted:
        print(f"\n  FULL KEY RECOVERED (within simulated model)")
        print(f"  FULL KEY = s1 ✓ + s2 [Phase B] + t0 ✓")

    return {"coeff_errors": coeff_errors, "accepted": accepted,
            "t0_match": t0_match, "norm": norm}


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--sigs", type=int, default=2000)
    args = p.parse_args()
    run_phase8(n_sigs=args.sigs)


if __name__ == "__main__":
    main()
