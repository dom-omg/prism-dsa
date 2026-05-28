#!/usr/bin/env python3
"""
Phase B — s2 Recovery via EP-Boundary Constraints (ML-DSA-44)

MECHANISM:
    EP[k][j] = (A·z)[k][j] - (c·C)[k][j]     where C = A·s1 - t1·2^D (known from s1)
    EP = w + c·t1·2^D  (NOT w itself — EP is unconstrained by signing rejection)

    Tight constraint fires when |LowBits(EP)[j]| ∈ (γ₂ - β, γ₂].
    Fires for both HB_Az == HB_EP and HB_Az != HB_EP cases (no boundary filter needed).
    Oracle accuracy: ~95% for both cases (h-direct correlation).

    dist = γ₂ - |EP_mod|; tight when dist < β (= τ·η = 78)

    EP_mod > 0 (near upper boundary):  h=1 → cs2 < -dist;  h=0 → cs2 ≥ -dist
    EP_mod < 0 (near lower boundary):  h=1 → cs2 > dist;   h=0 → cs2 ≤ dist

ORACLE ACCURACY vs ALTERNATIVE BASES:
    Az-based with HB_Az==HB_EP filter: 90.9% accuracy, 47% constraints discarded
    EP-based (this code):              95.0% accuracy, 0% constraints discarded
    W-based (w = V_ver - cC):          never fires — signing rejection |LB(w-cs2)| < γ₂-β
                                       ensures |LB(w)| ≤ γ₂-β-1 < tight threshold

MLWE HARDNESS NOTE (ML-DSA-44, η=2):
    The score for zero estimate (cs2=0) exceeds score for true s2 because h=0
    dominates (93.4% of tight constraints) and zero trivially satisfies cs2≥-dist.
    score(zero) ≈ 93.4% > score(true s2) ≈ 90.6% — coordinate descent cannot converge.
    This corresponds to SNR ≈ 0.27 from the interval-oracle analysis.
    For larger η (ML-DSA-65/87 with η=4), the SNR improves; recovery may be feasible.

    Coordinate descent is kept for completeness and as Phase B framework.
"""

import sys, time
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dilithium_py.ml_dsa.ml_dsa import ML_DSA
from dilithium_py.ml_dsa.default_parameters import DEFAULT_PARAMETERS
from phase7_syndrome_recovery import (
    flatten_module, rebuild_s1_from_coeff, center, ETA, D, L, K as K_POLYS,
)

Q = 8380417
N = 256


def decompose_scalar(r, alpha, q):
    r = int(r) % q
    r0 = r % alpha
    if r0 > alpha // 2:
        r0 -= alpha
    if (r - r0) % q == q - 1:
        r0 -= 1
        r1 = 0
    else:
        r1 = (r - r0) // alpha
    return r1, r0


def highbits_scalar(r, alpha, q):
    return decompose_scalar(r, alpha, q)[0]


def poly_mul_neg(a: np.ndarray, b: np.ndarray, q: int) -> np.ndarray:
    """Negacyclic polynomial multiply: a × b in Z_q[X]/(X^N+1)."""
    n = len(a)
    idx  = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    sign = np.where(np.arange(n)[:, None] >= np.arange(n)[None, :], 1, -1)
    mat  = (sign * a[idx]).astype(np.int64)
    return (mat @ b.astype(np.int64)) % q


def negacyclic_matrix(c_coeffs: np.ndarray) -> np.ndarray:
    """Build (N, N) negacyclic convolution matrix for polynomial c."""
    n = len(c_coeffs)
    idx  = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
    sign = np.where(np.arange(n)[:, None] >= np.arange(n)[None, :], 1, -1)
    return (sign * c_coeffs[idx]).astype(np.int64)


def collect_ep_constraints(dsa, sk, A_hat, C_mat, n_sigs, poly_idx, beta_s2):
    """
    Collect EP-based tight constraints for one polynomial.

    C_mat[k][j] = (A·s1)[k][j] - t1[k][j]·2^D  (coefficient domain, centered)

    Returns:
        ep_mod_arr: (M,) EP_mod values (LowBits of EP centred)
        dist_arr:   (M,) distance from |EP_mod| to gamma2 — tight when < beta_s2
        h_arr:      (M,) h bits (0 or 1)
        A_mat:      (M, N) negacyclic rows for cs2 = (c·s2[k])[j]
    """
    gamma2 = dsa.gamma_2
    alpha  = 2 * gamma2
    offset = poly_idx * N

    max_c = max(50000, n_sigs * N * 3 // 1000)
    ep_mod_arr = np.empty(max_c, dtype=np.int64)
    dist_arr   = np.empty(max_c, dtype=np.int64)
    h_arr      = np.empty(max_c, dtype=np.int8)
    A_mat      = np.empty((max_c, N), dtype=np.int64)
    ptr = 0

    t0 = time.perf_counter()
    for i in range(n_sigs):
        msg = f"phB-{i:06d}".encode()
        sig = dsa.sign(sk, msg)
        c_tilde, z_mod, h_mod = dsa._unpack_sig(sig)
        c_poly   = dsa.R.sample_in_ball(c_tilde, dsa.tau)
        c_coeffs = np.array(c_poly.coeffs, dtype=np.int64)

        # A·z for this polynomial row k
        Az_hat = A_hat @ z_mod.to_ntt()
        Az_flat = flatten_module(Az_hat.from_ntt())   # K*N coefficient domain
        Az_row  = center(Az_flat)[offset:offset + N]  # centered (K,)[j]

        h_flat = flatten_module(h_mod).reshape(K_POLYS, N)[poly_idx]  # (N,)

        # c·C[k] = negacyclic mult of c with C[poly_idx]
        C_k    = C_mat[poly_idx]  # (N,) centered
        cC_k   = center(poly_mul_neg(c_coeffs, np.asarray(C_k, dtype=np.int64), Q))  # (N,)

        # EP[k][j] = Az[k][j] - cC[k][j]
        EP_row = center(np.array(Az_row, dtype=np.int64) - np.array(cC_k, dtype=np.int64))  # centered

        # Negacyclic matrix for this challenge
        neg_mat = negacyclic_matrix(c_coeffs)  # (N, N)

        for j in range(N):
            ep = int(EP_row[j])
            # EP_mod = LowBits(EP, alpha) centred
            ep_mod = ep % alpha
            if ep_mod > alpha // 2:
                ep_mod -= alpha
            dist = int(gamma2) - abs(ep_mod)  # > 0 when near boundary

            if dist >= 0 and dist < beta_s2:
                # Tight constraint
                if ptr >= max_c:
                    extra = max_c
                    ep_mod_arr = np.concatenate([ep_mod_arr, np.empty(extra, dtype=np.int64)])
                    dist_arr   = np.concatenate([dist_arr,   np.empty(extra, dtype=np.int64)])
                    h_arr      = np.concatenate([h_arr,      np.empty(extra, dtype=np.int8)])
                    A_mat      = np.vstack([A_mat, np.empty((extra, N), dtype=np.int64)])
                    max_c += extra

                ep_mod_arr[ptr] = ep_mod
                dist_arr[ptr]   = dist
                h_arr[ptr]      = int(h_flat[j])
                A_mat[ptr]      = neg_mat[j]
                ptr += 1

        if (i + 1) % 500 == 0:
            dt = time.perf_counter() - t0
            print(f"    [{i+1}/{n_sigs}]  tight={ptr}  ({dt:.1f}s)")

    dt = time.perf_counter() - t0
    print(f"  poly {poly_idx}: {ptr} tight constraints from {n_sigs} sigs ({dt:.1f}s)")
    return (ep_mod_arr[:ptr].copy(), dist_arr[:ptr].copy(),
            h_arr[:ptr].copy(),      A_mat[:ptr].copy())


def constraint_valid(cs2_val, ep_mod, dist, h_bit, beta_s2):
    """
    Check if cs2_val ∈ [-beta_s2, beta_s2] satisfies the EP-boundary constraint.

    EP_mod > 0 (near upper boundary +gamma2):
        h=1: cs2 < -dist   (crossing → cs2 < -dist)
        h=0: cs2 >= -dist
    EP_mod < 0 (near lower boundary -gamma2):
        h=1: cs2 > dist
        h=0: cs2 <= dist
    """
    if ep_mod > 0:
        if h_bit == 1:
            return cs2_val < -dist
        else:
            return cs2_val >= -dist
    else:
        if h_bit == 1:
            return cs2_val > dist
        else:
            return cs2_val <= dist


def batch_constraint_score(cs2_arr, ep_mod_arr, dist_arr, h_arr, beta_s2):
    """
    Vectorised constraint scoring.
    cs2_arr: (M,) current cs2 values (one per constraint)
    Returns: fraction satisfied
    """
    pos_mask = ep_mod_arr > 0

    # h=1, EP_mod > 0: cs2 < -dist
    valid_pos_h1 = pos_mask & (h_arr == 1) & (cs2_arr < -dist_arr)
    # h=0, EP_mod > 0: cs2 >= -dist
    valid_pos_h0 = pos_mask & (h_arr == 0) & (cs2_arr >= -dist_arr)
    # h=1, EP_mod < 0: cs2 > dist
    valid_neg_h1 = ~pos_mask & (h_arr == 1) & (cs2_arr > dist_arr)
    # h=0, EP_mod < 0: cs2 <= dist
    valid_neg_h0 = ~pos_mask & (h_arr == 0) & (cs2_arr <= dist_arr)

    return int((valid_pos_h1 | valid_pos_h0 | valid_neg_h1 | valid_neg_h0).sum())


def coordinate_descent(ep_mod_arr, dist_arr, h_arr, A_mat, beta_s2,
                        n_iters=8, true_s2=None):
    """
    Recover s2[k] via coordinate descent on EP-boundary constraints.
    """
    M = len(ep_mod_arr)
    s2_est = np.zeros(N, dtype=np.int64)
    cs2_cur = np.zeros(M, dtype=np.int64)  # A_mat @ s2_est, centered

    t0 = time.perf_counter()
    for it in range(n_iters):
        n_changed = 0
        for l in range(N):
            a_col   = A_mat[:, l]          # (M,)
            base_cs2 = cs2_cur             # current, before subtracting contribution of l

            # Remove contribution of s2_est[l] from cs2_cur
            base_minus_l = center((cs2_cur - int(s2_est[l]) * a_col) % Q)

            best_v     = int(s2_est[l])
            best_score = -1
            for v in range(-ETA, ETA + 1):
                cs2_test  = center((base_minus_l + v * a_col) % Q)
                score     = batch_constraint_score(cs2_test, ep_mod_arr, dist_arr, h_arr, beta_s2)
                if score > best_score:
                    best_score = score
                    best_v = v

            if best_v != s2_est[l]:
                cs2_cur   = center((base_minus_l + best_v * a_col) % Q)
                s2_est[l] = best_v
                n_changed += 1

        score_total = batch_constraint_score(cs2_cur, ep_mod_arr, dist_arr, h_arr, beta_s2)
        errors = int(((s2_est - true_s2) % Q != 0).sum()) if true_s2 is not None else -1
        print(f"  iter {it+1}: changed={n_changed:3d}  score={score_total}/{M}  "
              f"errors={'?' if errors < 0 else errors}  ({time.perf_counter()-t0:.1f}s)")
        if n_changed == 0:
            break

    return s2_est


def run_phase_b(n_sigs=2000):
    dsa = ML_DSA(DEFAULT_PARAMETERS["ML_DSA_44"])
    gamma2   = dsa.gamma_2   # 95232
    alpha    = 2 * gamma2
    beta_s2  = dsa.tau * dsa.eta   # τ·η = 39·2 = 78

    pk, sk = dsa.keygen()
    rho, t1 = dsa._unpack_pk(pk)
    _, _, _, s1, s2, _ = dsa._unpack_sk(sk)
    A_hat = dsa._expand_matrix_from_seed(rho)

    s2_true = center(flatten_module(s2)).reshape(K_POLYS, N)

    # Compute C = A·s1 - t1·2^D  (coefficient domain, K×N)
    As1_flat  = center(flatten_module((A_hat @ s1.to_ntt()).from_ntt()))
    t1_flat   = np.array(flatten_module(t1), dtype=np.int64)
    C_flat    = center(As1_flat - t1_flat * (1 << D))
    C_mat     = C_flat.reshape(K_POLYS, N)  # (K, N) centered

    print("=" * 60)
    print("PHASE B — s2 Recovery via EP-Boundary Oracle (ML-DSA-44)")
    print("=" * 60)
    print(f"  n_sigs={n_sigs}, β_s2={beta_s2}, γ₂={gamma2}")
    print(f"  Tight threshold: |EP_mod| ∈ ({gamma2-beta_s2}, {gamma2}]")

    s2_recovered = np.zeros((K_POLYS, N), dtype=np.int64)

    for k in range(K_POLYS):
        print(f"\n[poly {k}] Collecting constraints...")
        ep_mod, dist, h, A_mat = collect_ep_constraints(
            dsa, sk, A_hat, C_mat, n_sigs, poly_idx=k, beta_s2=beta_s2)

        print(f"[poly {k}] Coordinate descent ({len(ep_mod)} constraints)...")
        s2_k = coordinate_descent(ep_mod, dist, h, A_mat, beta_s2,
                                   n_iters=10, true_s2=s2_true[k])
        errors = int(((s2_k - s2_true[k]) % Q != 0).sum())
        print(f"  → poly {k}: {N - errors}/256 correct  ({errors} errors)")
        s2_recovered[k] = s2_k

    total_errors = int(((s2_recovered - s2_true) % Q != 0).sum())
    print(f"\n{'='*60}")
    print(f"s2 TOTAL: {K_POLYS*N - total_errors}/{K_POLYS*N}  ({total_errors} errors)")
    return s2_recovered, total_errors


def verify_oracle_accuracy(n_sigs: int = 300) -> None:
    """Verify EP-based oracle accuracy for both HB_Az==HB_EP and HB_Az!=HB_EP."""
    dsa = ML_DSA(DEFAULT_PARAMETERS["ML_DSA_44"])
    gamma2 = dsa.gamma_2
    alpha  = 2 * gamma2
    beta_s2 = dsa.tau * dsa.eta

    pk, sk = dsa.keygen()
    rho, t1 = dsa._unpack_pk(pk)
    _, _, _, s1, s2, _ = dsa._unpack_sk(sk)
    A_hat = dsa._expand_matrix_from_seed(rho)
    As1_flat  = center(flatten_module((A_hat @ s1.to_ntt()).from_ntt()))
    t1_flat   = np.array(flatten_module(t1), dtype=np.int64)
    C_flat    = center(As1_flat - t1_flat * (1 << D))
    s2_true   = center(flatten_module(s2)).reshape(K_POLYS, N)

    print(f"=== Oracle verification: {n_sigs} sigs, k=0 ===")
    stats = {"ep_HB_same": [0, 0], "ep_HB_diff": [0, 0], "zero_score": [0, 0]}
    k = 0; offset = k * N

    for i in range(n_sigs):
        msg = f"verify-{i:06d}".encode()
        sig = dsa.sign(sk, msg)
        c_tilde, z_mod, h_mod = dsa._unpack_sig(sig)
        c_poly   = dsa.R.sample_in_ball(c_tilde, dsa.tau)
        c_coeffs = np.array(c_poly.coeffs, dtype=np.int64)
        Az_flat  = center(flatten_module((A_hat @ z_mod.to_ntt()).from_ntt()))
        Az_row   = Az_flat[offset:offset + N]
        h_flat   = flatten_module(h_mod).reshape(K_POLYS, N)[k]
        C_k      = C_flat[offset:offset + N]
        cC_k     = center(poly_mul_neg(c_coeffs, np.asarray(C_k, np.int64), Q))
        EP_row   = center(Az_row.astype(np.int64) - cC_k)
        cs2_true_row = center(poly_mul_neg(c_coeffs, s2_true[k], Q))

        az_mod = Az_row % alpha; az_mod = np.where(az_mod > alpha//2, az_mod - alpha, az_mod)
        ep_mod = EP_row % alpha; ep_mod = np.where(ep_mod > alpha//2, ep_mod - alpha, ep_mod)
        hb_Az, _ = decompose_scalar_arr(Az_row % Q, alpha, Q)
        hb_EP, _ = decompose_scalar_arr(EP_row % Q, alpha, Q)
        dist_ep   = gamma2 - np.abs(ep_mod)
        tight     = (dist_ep >= 0) & (dist_ep < beta_s2)
        same_hb   = (hb_Az == hb_EP)

        for j in np.where(tight)[0]:
            ep = int(ep_mod[j]); d = int(dist_ep[j])
            hb = int(h_flat[j]); cs = int(cs2_true_row[j])
            is_same = bool(same_hb[j])
            if ep > 0:
                sat_true = (hb == 1 and cs < -d) or (hb == 0 and cs >= -d)
                sat_zero = (hb == 0)  # cs2=0 always ≥ -d
            else:
                sat_true = (hb == 1 and cs > d)  or (hb == 0 and cs <= d)
                sat_zero = (hb == 0)  # cs2=0 always ≤ d
            key = "ep_HB_same" if is_same else "ep_HB_diff"
            stats[key][0] += 1; stats[key][1] += int(sat_true)
            stats["zero_score"][0] += 1; stats["zero_score"][1] += int(sat_zero)

    for name, (tot, sat) in stats.items():
        if name != "zero_score":
            print(f"  {name}: {sat}/{tot} = {sat/max(tot,1)*100:.1f}% accuracy")
    tot, sat = stats["zero_score"]
    print(f"  zero_score: {sat}/{tot} = {sat/max(tot,1)*100:.1f}% (vs true s2 above)")
    print(f"  [HB_Az!=HB_EP is NOT filtered — EP handles both cases correctly]")


def decompose_scalar_arr(r: np.ndarray, alpha: int, q: int):
    r = np.asarray(r, np.int64) % q
    r0 = r % alpha
    r0 = np.where(r0 > alpha // 2, r0 - alpha, r0)
    r1 = (r - r0) // alpha
    boundary = (r - np.where(r0 < 0, r0 + alpha, r0)) % q == q - 1
    r0 = np.where(boundary, r0 - 1, r0)
    r1 = np.where(boundary, 0, r1)
    return r1, r0


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--sigs", type=int, default=2000)
    p.add_argument("--verify", action="store_true", help="Run oracle accuracy verification")
    args = p.parse_args()
    if args.verify:
        verify_oracle_accuracy(n_sigs=args.sigs)
    else:
        run_phase_b(n_sigs=args.sigs)
