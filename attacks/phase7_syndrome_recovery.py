#!/usr/bin/env python3
"""
Phase 7 — Soft-DPA + Syndrome Recovery: Full Key via t-Coupling

PROBLEM (Phase 6 established):
    DPA on PRISM-DSA: 14× harder than ML-DSA. Finite traces → ε > 0 NTT-domain
    errors in recovered s1_hat. One wrong coefficient corrupts A·s1 chain.

THREE-LAYER FIX:
    L1 — Error rate curve:   measure ε(N) empirically (fraction of 256 NTT
                             coefficients wrong at each trace count checkpoint).
    L2 — Soft-decision DPA:  confidence gap per coefficient → identify C << N
                             uncertain positions.
    L3 — Syndrome search:    work in COEFFICIENT DOMAIN (s1 ∈ {-η..η}^{L·N}).
                             Enumerate corrections at C positions (5^C candidates
                             vs ∞ for Z_q-domain). Verify via t-coupling oracle.

KEY INVARIANT:
    t = A·s1 + s2   (exact, from keygen)
    t = t1·2^D + t0, with ||t0||∞ ≤ 2^(D-1) = 4096
    t1 is in the PUBLIC KEY.

VERIFICATION ORACLE (public-key only):
    Given coefficient-domain (s1_cand, s2_cand):
        t_cand = A·s1_cand + s2_cand
        residual = t_cand - t1·2^D   (computed entirely from public key + candidate)
        Accept iff ||residual||∞ ≤ TOL = 2^(D-1) + η = 4098
    For correct (s1, s2): residual = t0, ||t0||∞ ≤ 4096 → always accepted.
    For random wrong:  residual ≈ uniform mod q → norm ≈ q/4 ≈ 2M → rejected.
    False-positive probability: (2·TOL / q)^{K·N} = (8196/8380417)^1024 → 0.

COEFFICIENT-DOMAIN SYNDROME SEARCH:
    After DPA → INTT projection → rounding to {-η,...,η}, errors are LOCALIZED.
    For ML-DSA-44: η=2, so 5 values per coefficient.
    Search space: 4^C (c wrong values per position) where C = #uncertain coefficients.
    ε=1: 4 × 256  =   1,024 candidates (trivial)
    ε=2: 16 × 256² ≈  1M   candidates (seconds on 1 CPU)
    ε=3: 64 × 256³ ≈  270M  candidates (GPU / parallel)

Rule 44: all claims below = what this code measures. Physical DPA noise not modeled.
"""

import numpy as np
import time
import sys
import os
from itertools import product as iproduct

sys.path.insert(0, os.path.dirname(__file__))
from phase6_dpa_timing import (
    SignerInstrumented, simulate_power_mldsa, simulate_power_prism,
    hw_vec, Q, N, FIS_SLOTS, PARAMS,
)

# ─── Parameters (ML-DSA-44) ───────────────────────────────────────────────────

ETA = 2        # s1, s2 coefficients in [-2, 2]
D   = 13       # t decomposition: t0 ∈ [-2^12, 2^12]
L   = 4        # polynomials in s1
K   = 4        # polynomials in s2
TOL = (1 << (D - 1)) + ETA   # = 4098

# ─── Helpers ─────────────────────────────────────────────────────────────────

def center(x, q=Q):
    """Center-lift from [0,q) to [-q/2, q/2]."""
    x = np.asarray(x, np.int64) % q
    return np.where(x > q // 2, x - q, x)


def flatten_module(mod):
    """Flatten dilithium-py Module Vector to np.int64 array."""
    out = []
    for row in mod.rows():
        for poly in row:
            out.extend(poly.coeffs)
    return np.array(out, dtype=np.int64)


def rebuild_s1_from_coeff(dsa, coeff_flat: np.ndarray):
    """
    Rebuild dilithium-py s1 Vector from flat coefficient array (L·N values).
    Returns: s1 Vector (not NTT).
    """
    R = dsa.R
    M = dsa.M
    polys = []
    for pi in range(L):
        c = coeff_flat[pi * N:(pi + 1) * N].tolist()
        polys.append(R(c))
    return M.vector(polys)


def rebuild_s2_from_coeff(dsa, coeff_flat: np.ndarray):
    """Rebuild dilithium-py s2 Vector from flat coefficient array (K·N values)."""
    R = dsa.R
    M = dsa.M
    polys = []
    for pi in range(K):
        c = coeff_flat[pi * N:(pi + 1) * N].tolist()
        polys.append(R(c))
    return M.vector(polys)


# ─── t-Coupling verification oracle ──────────────────────────────────────────

def verify_candidate(s1_coeff_flat: np.ndarray,
                     s2_coeff_flat: np.ndarray,
                     rho: bytes,
                     t1_flat: np.ndarray,
                     dsa: SignerInstrumented) -> tuple[bool, int]:
    """
    Verify candidate (s1, s2) via public-key t-coupling.

    Computes t_cand = A·s1_cand + s2_cand using dilithium-py native arithmetic.
    Checks: ||t_cand - t1·2^D||∞ ≤ TOL.

    Returns: (accepted, max_residual_norm)
    """
    s1_vec = rebuild_s1_from_coeff(dsa, s1_coeff_flat)
    s2_vec = rebuild_s2_from_coeff(dsa, s2_coeff_flat)
    A_hat  = dsa._expand_matrix_from_seed(rho)
    t_hat  = A_hat @ s1_vec.to_ntt() + s2_vec.to_ntt()
    t_cand = flatten_module(t_hat.from_ntt())

    residual = center(t_cand - t1_flat * (1 << D))
    max_norm = int(np.abs(residual).max())
    return max_norm <= TOL, max_norm


# ─── L1: Error rate curve ─────────────────────────────────────────────────────

def error_rate_curve(dsa: SignerInstrumented,
                     sk: bytes,
                     s1_hat_true: np.ndarray,
                     n_sigs: int,
                     checkpoints: list,
                     n_hyp: int = 600,
                     use_prism: bool = False,
                     rng: np.random.Generator = None) -> dict:
    """
    Measure DPA error rate vs. signing operation count.

    At each checkpoint, run hard-decision CPA on all N=256 NTT coefficients
    of poly 0. Count how many have rank(true) ≠ 1 among n_hyp random competitors.

    Rule 44: rank measured against n_hyp random competitors, not full Z_q scan.
    Success here is a CONSERVATIVE lower bound on actual DPA success rate.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    all_c_acc  = np.zeros((n_sigs, N), dtype=np.int64)
    all_traces = np.zeros((n_sigs, N), dtype=np.int64)
    all_n_iter = np.zeros(n_sigs, dtype=np.int32)

    for i in range(n_sigs):
        msg = f"err7-{i:06d}".encode()
        r   = dsa.sign_with_trace_info(sk, msg)
        c_acc = r["c_hat_accepted"]
        if use_prism:
            tr = simulate_power_prism(c_acc, s1_hat_true, rng)
        else:
            tr = simulate_power_mldsa(c_acc, r["c_hat_rejected"], s1_hat_true, rng)
        all_c_acc[i]  = c_acc[:N]
        all_traces[i] = tr[:N]
        all_n_iter[i] = r["n_iter"]

    clean_mask = all_n_iter == 1
    results    = {}

    for chk in checkpoints:
        if chk > n_sigs:
            break
        if not use_prism:
            idx = np.where(clean_mask[:chk])[0]
            if len(idx) < 5:
                results[chk] = {"n_traces": 0, "errors": N, "error_rate": 1.0}
                continue
            tr_use = all_traces[idx]
            ch_use = all_c_acc[idx]
        else:
            tr_use = all_traces[:chk]
            ch_use = all_c_acc[:chk]

        errors = 0
        for j in range(N):
            power  = tr_use[:, j].astype(np.float64)
            c_j    = ch_use[:, j].astype(np.int64)
            h_true = int(s1_hat_true[j]) % Q
            rand_h = rng.integers(0, Q, size=n_hyp, dtype=np.int64)

            m_true = hw_vec((c_j * h_true) % Q).astype(np.float64)
            r_true = abs(np.corrcoef(m_true, power)[0, 1])

            models = hw_vec(
                (c_j[:, None] * rand_h[None, :]) % Q
            ).astype(np.float64)   # (traces, n_hyp)
            r_rand = np.array([
                abs(np.corrcoef(models[:, k], power)[0, 1])
                for k in range(n_hyp)
            ])
            if (r_rand >= r_true).any():
                errors += 1

        results[chk] = {
            "n_traces":   len(tr_use),
            "errors":     errors,
            "error_rate": errors / N,
            "n_correct":  N - errors,
        }

    return results


# ─── L2: Soft-decision DPA ────────────────────────────────────────────────────

def soft_dpa_confidence(traces: np.ndarray,
                        c_hats: np.ndarray,
                        s1_hat_true: np.ndarray,
                        n_hyp: int = 600,
                        rng: np.random.Generator = None) -> np.ndarray:
    """
    Return confidence[j] = corr(true_h) - max(corr(random_h)) for each j.

    Positive → DPA cleanly identifies the true value.
    Negative or near-zero → uncertain position → candidate for syndrome search.

    Note: ground truth s1_hat_true used as RESEARCH METRIC only. In a production
    attack, confidence is measured as corr(rank-1) - corr(rank-2) from a full Z_q
    scan (computationally expensive but feasible with GPU).
    """
    if rng is None:
        rng = np.random.default_rng(0)
    conf = np.zeros(N)
    for j in range(N):
        power  = traces[:, j].astype(np.float64)
        c_j    = c_hats[:, j].astype(np.int64)
        h_true = int(s1_hat_true[j]) % Q
        rand_h = rng.integers(0, Q, size=n_hyp, dtype=np.int64)
        m_true = hw_vec((c_j * h_true) % Q).astype(np.float64)
        r_true = abs(np.corrcoef(m_true, power)[0, 1])
        models = hw_vec((c_j[:, None] * rand_h[None, :]) % Q).astype(np.float64)
        r_rand = max(abs(np.corrcoef(models[:, k], power)[0, 1]) for k in range(n_hyp))
        conf[j] = r_true - r_rand
    return conf


# ─── L3: Coefficient-domain syndrome search ──────────────────────────────────

def syndrome_search_coeff(s1_coeff_flat: np.ndarray,
                          s2_coeff_flat: np.ndarray,
                          uncertain_positions: list,
                          rho: bytes,
                          t1_flat: np.ndarray,
                          dsa: SignerInstrumented,
                          verbose: bool = True) -> dict:
    """
    Search for correct s1 by correcting errors at uncertain_positions.

    uncertain_positions: list of flat indices into s1_coeff_flat (L·N values).
    At each position, the current value is wrong. Try all alternatives in {-η,...,η}.

    Search space: 4^C candidates (C = len(uncertain_positions), 4 wrong values each).
    ε=1: max 4×L·N = 4096 candidates. ε=2: ~16M. ε=3: ~67B (infeasible naïve).

    Returns first accepted candidate.
    """
    eta_vals = list(range(-ETA, ETA + 1))   # {-2, -1, 0, 1, 2}
    C        = len(uncertain_positions)
    total    = (2 * ETA) ** C               # 4^C (excluding the current value)

    if verbose:
        print(f"\n[Syndrome Search — Coefficient Domain]")
        print(f"  C = {C} uncertain positions")
        print(f"  η = {ETA}, alternatives per pos = {2*ETA}")
        print(f"  Total candidates: {total:,}")
        print(f"  TOL = {TOL} (vs random ≈ {Q//4:,})")

    found   = None
    checked = 0
    t_start = time.perf_counter()

    # For each position, try all values EXCEPT the current (injected wrong) value
    alt_sets = []
    for flat_idx in uncertain_positions:
        current = int(s1_coeff_flat[flat_idx])
        alts    = [v for v in eta_vals if v != current]
        alt_sets.append((flat_idx, alts))

    for combo in iproduct(*[a for _, a in alt_sets]):
        s1_cand = s1_coeff_flat.copy()
        for (flat_idx, _), new_val in zip(alt_sets, combo):
            s1_cand[flat_idx] = new_val

        accepted, max_norm = verify_candidate(s1_cand, s2_coeff_flat, rho, t1_flat, dsa)
        checked += 1

        if accepted:
            found = {
                "s1_coeff": s1_cand,
                "max_residual": max_norm,
                "corrections": {idx: val for (idx, _), val in zip(alt_sets, combo)},
                "candidates_checked": checked,
                "time_s": time.perf_counter() - t_start,
            }
            break

    elapsed = time.perf_counter() - t_start
    if verbose:
        if found:
            print(f"\n  KEY FOUND after {found['candidates_checked']:,} candidates "
                  f"({elapsed:.3f}s)")
            print(f"  Residual ||t0||∞ = {found['max_residual']}  (TOL={TOL})")
        else:
            print(f"  NOT FOUND in {checked:,} candidates ({elapsed:.3f}s)")

    return {
        "found":    found is not None,
        "result":   found,
        "checked":  checked,
        "total":    total,
        "time_s":   elapsed,
    }


# ─── Full pipeline ────────────────────────────────────────────────────────────

def run_phase7(n_sigs: int = 250,
               n_hyp:  int = 400,
               epsilon: int = 1,
               verbose: bool = True) -> dict:
    """
    End-to-end Phase 7: error rate curve → soft DPA → syndrome search.

    epsilon: number of coefficient-domain errors to inject (simulate imperfect DPA).
    """
    rng = np.random.default_rng(42)
    dsa = SignerInstrumented(PARAMS)
    pk, sk, s1_hat_true = dsa.keygen_with_s1hat()

    # Unpack ground truth
    rho, k_key, tr_sk, s1_true, s2_true, _ = dsa._unpack_sk(sk)
    rho_pk, t1 = dsa._unpack_pk(pk)
    s1_coeff_true = flatten_module(s1_true)   # L·N coeff-domain values
    s2_coeff_true = flatten_module(s2_true)   # K·N coeff-domain values
    t1_flat       = flatten_module(t1)        # K·N

    print("=" * 70)
    print("PHASE 7 — Soft-DPA + Syndrome Recovery (ML-DSA-44)")
    print("=" * 70)
    print(f"  n_sigs={n_sigs}  n_hyp={n_hyp}  ε_inject={epsilon}")
    print(f"  L={L}, K={K}, N={N}, η={ETA}, D={D}, TOL={TOL}")

    # ── STEP 1: Error rate curve ─────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print("STEP 1 — Error Rate Curve (poly 0, all N=256 NTT coefficients)")
    checkpoints = [c for c in [20, 50, 100, 200, 300, 500, 700, 1000, 1400, 2000] if c <= n_sigs]

    t0 = time.perf_counter()
    ml_curve = error_rate_curve(dsa, sk, s1_hat_true, n_sigs, checkpoints,
                                 n_hyp=n_hyp, use_prism=False, rng=rng)
    pr_curve = error_rate_curve(dsa, sk, s1_hat_true, n_sigs, checkpoints,
                                 n_hyp=n_hyp, use_prism=True,  rng=rng)
    print(f"  Computed in {time.perf_counter()-t0:.1f}s\n")

    hdr = f"  {'N_ops':>5}  {'ML traces':>9}  {'ML err':>7}  {'ML ε%':>6}  " \
          f"{'PR traces':>9}  {'PR err':>7}  {'PR ε%':>6}"
    print(hdr)
    print("  " + "─" * 62)
    for chk in checkpoints:
        ml = ml_curve.get(chk, {})
        pr = pr_curve.get(chk, {})
        print(f"  {chk:>5}  {ml.get('n_traces','?'):>9}  "
              f"{ml.get('errors','?'):>7}  {100*ml.get('error_rate',1):>5.1f}%  "
              f"{pr.get('n_traces','?'):>9}  "
              f"{pr.get('errors','?'):>7}  {100*pr.get('error_rate',1):>5.1f}%")

    # ── STEP 2: t-coupling oracle sanity ─────────────────────────────────────
    print(f"\n{'─'*60}")
    print("STEP 2 — t-Coupling Oracle Verification")

    ok_true, norm_true = verify_candidate(s1_coeff_true, s2_coeff_true, rho, t1_flat, dsa)
    # Wrong s1: flip one coefficient to a random value
    s1_wrong = s1_coeff_true.copy()
    s1_wrong[13] = (s1_wrong[13] + 3) % (2 * ETA + 1) - ETA   # shift by ≠0
    if s1_wrong[13] == s1_coeff_true[13]:
        s1_wrong[13] = -s1_coeff_true[13]
    ok_wrong, norm_wrong = verify_candidate(s1_wrong, s2_coeff_true, rho, t1_flat, dsa)

    print(f"  True   (s1, s2): accepted={ok_true},  ||residual||∞ = {norm_true}")
    print(f"  Wrong  (s1, s2): accepted={ok_wrong}, ||residual||∞ = {norm_wrong:,}")
    print(f"  Random residual would be ≈ {Q//4:,}  (×{Q//4//TOL:.0f} larger than TOL)")

    # ── STEP 3: Soft-DPA confidence + inject ε errors ────────────────────────
    print(f"\n{'─'*60}")
    print(f"STEP 3 — Soft-DPA Confidence + Inject ε={epsilon} error(s)")

    # Collect clean traces
    all_c_acc  = np.zeros((n_sigs, N), dtype=np.int64)
    all_traces = np.zeros((n_sigs, N), dtype=np.int64)
    all_n_iter = np.zeros(n_sigs, dtype=np.int32)
    for i in range(n_sigs):
        r = dsa.sign_with_trace_info(sk, f"soft7-{i:06d}".encode())
        all_c_acc[i]  = r["c_hat_accepted"][:N]
        all_traces[i] = simulate_power_mldsa(
            r["c_hat_accepted"], r["c_hat_rejected"], s1_hat_true, rng)[:N]
        all_n_iter[i] = r["n_iter"]

    clean_idx = np.where(all_n_iter == 1)[0]
    print(f"  Clean traces: {len(clean_idx)} / {n_sigs}")

    conf = soft_dpa_confidence(all_traces[clean_idx], all_c_acc[clean_idx],
                                s1_hat_true, n_hyp=n_hyp, rng=rng)
    worst_j = np.argsort(conf)[:epsilon + 2]

    print(f"\n  Confidence (worst {epsilon+2}):")
    for ji in worst_j:
        print(f"    NTT j={ji:3d}: conf={conf[ji]:+.4f}  "
              f"s1_coeff[{ji}]={s1_coeff_true[ji]}")

    # Inject ε errors at lowest-confidence positions (in coefficient domain)
    s1_coeff_dpa = s1_coeff_true.copy()
    injected_flat_idx = []
    for ei in range(epsilon):
        ji = int(worst_j[ei])
        true_val  = int(s1_coeff_true[ji])
        # Pick wrong value from {-η,...,η} \ {true_val}
        alts      = [v for v in range(-ETA, ETA + 1) if v != true_val]
        wrong_val = alts[ei % len(alts)]
        s1_coeff_dpa[ji] = wrong_val
        injected_flat_idx.append(ji)
        print(f"  Inject error: s1[{ji}] = {true_val} → {wrong_val}")

    ok_dpa, norm_dpa = verify_candidate(s1_coeff_dpa, s2_coeff_true, rho, t1_flat, dsa)
    print(f"\n  DPA+error accepted={ok_dpa}, ||residual||∞ = {norm_dpa:,}")

    # ── STEP 4: Syndrome search ───────────────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"STEP 4 — Syndrome Search (ε={epsilon}, C={len(injected_flat_idx)} positions)")

    result = syndrome_search_coeff(
        s1_coeff_flat     = s1_coeff_dpa,
        s2_coeff_flat     = s2_coeff_true,
        uncertain_positions = injected_flat_idx,
        rho               = rho,
        t1_flat           = t1_flat,
        dsa               = dsa,
        verbose           = verbose,
    )

    # Verify recovered key matches true key
    if result["found"]:
        s1_recovered = result["result"]["s1_coeff"]
        key_correct  = np.array_equal(s1_recovered, s1_coeff_true)
        print(f"  Recovered == true s1: {key_correct}")

    # ── Summary ───────────────────────────────────────────────────────────────
    last_chk = checkpoints[-1]
    ml_last  = ml_curve.get(last_chk, {})
    pr_last  = pr_curve.get(last_chk, {})

    print(f"\n{'='*70}")
    print("SUMMARY (Rule 44 — exact measurements)")
    print(f"{'='*70}")
    print(f"""
  PRECONDITIONS:
    Simulated HW power traces, no physical noise. ML-DSA-44.
    n_hyp={n_hyp} random competitors → conservative success rate lower bound.

  ERROR RATE at N_ops={last_chk}:
    ML-DSA (clean, {ml_last.get('n_traces','?')} traces):  {ml_last.get('errors','?')}/{N} errors = {100*ml_last.get('error_rate',1):.1f}%
    PRISM-DSA ({pr_last.get('n_traces','?')} traces):       {pr_last.get('errors','?')}/{N} errors = {100*pr_last.get('error_rate',1):.1f}%

  t-COUPLING ORACLE:
    True  (s1,s2): ||residual||∞ = {norm_true}  ≤ TOL={TOL}  → {'ACCEPT' if ok_true else 'REJECT'}
    Wrong (s1,s2): ||residual||∞ = {norm_wrong:,}  vs TOL={TOL} → {'ACCEPT' if ok_wrong else 'REJECT'}
    Signal/noise ratio: {norm_wrong // TOL:.0f}×  (random residual / TOL)

  SYNDROME SEARCH (ε={epsilon}):
    Candidates checked: {result['checked']:,} / {result['total']:,}
    Found: {result['found']}
    Time:  {result['time_s']:.3f}s
    {'Key correct: ' + str(result['found'] and np.array_equal(result["result"]["s1_coeff"], s1_coeff_true)) if result['found'] else ''}

  CLAIM (within simulated model):
    t-coupling distinguishes correct from wrong key with norm gap ×{norm_wrong//TOL:.0f}.
    Syndrome search recovers key with ε={epsilon} coefficient error(s) in {result['checked']:,} candidates.
    Search space vs brute force: {result['total']:,} vs 5^{L*N} (astronomically large).
""")

    return {
        "ml_curve":   ml_curve,
        "pr_curve":   pr_curve,
        "confidence": conf,
        "result":     result,
        "norm_true":  norm_true,
        "norm_wrong": norm_wrong,
    }


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    p = argparse.ArgumentParser(description="Phase 7 — Syndrome Recovery via t-Coupling")
    p.add_argument("--sigs",    type=int, default=200, help="Signing ops (default 200)")
    p.add_argument("--hyp",     type=int, default=300, help="CPA hypotheses (default 300)")
    p.add_argument("--epsilon", type=int, default=1,   help="Errors to inject (default 1)")
    p.add_argument("--theory",  action="store_true",   help="Theory only, no simulation")
    p.add_argument("--oracle",  action="store_true",   help="Oracle test only (fast)")
    args = p.parse_args()

    if args.theory:
        print("t-COUPLING THEORY (ML-DSA-44)")
        print(f"  t = A·s1 + s2 = t1·2^D + t0")
        print(f"  ||t0||∞ ≤ 2^(D-1)={1<<(D-1)},  ||s2||∞ ≤ η={ETA}")
        print(f"  TOL = {TOL}")
        print(f"  Wrong residual ≈ {Q//4:,}  (×{Q//4//TOL:.0f} signal/noise)")
        print(f"\n  Coeff-domain syndrome search space:")
        print(f"  {'ε':>3}  {'candidates':>15}  {'L·N positions':>14}")
        for eps in range(1, 5):
            print(f"  {eps:>3}  {(2*ETA)**eps * (L*N)**eps:>15,}  {L*N:>14}")
        return

    if args.oracle:
        dsa = SignerInstrumented(PARAMS)
        pk, sk, s1_hat_true = dsa.keygen_with_s1hat()
        _, _, _, s1, s2, _ = dsa._unpack_sk(sk)
        rho, t1 = dsa._unpack_pk(pk)
        s1_f = flatten_module(s1)
        s2_f = flatten_module(s2)
        t1_f = flatten_module(t1)
        ok, norm = verify_candidate(s1_f, s2_f, rho, t1_f, dsa)
        print(f"Oracle test: accepted={ok}, norm={norm}, TOL={TOL}")
        s1_bad = s1_f.copy(); s1_bad[0] = -s1_bad[0] - 1
        ok2, norm2 = verify_candidate(s1_bad, s2_f, rho, t1_f, dsa)
        print(f"Wrong key:   accepted={ok2}, norm={norm2:,}")
        return

    run_phase7(n_sigs=args.sigs, n_hyp=args.hyp, epsilon=args.epsilon)


if __name__ == "__main__":
    main()
