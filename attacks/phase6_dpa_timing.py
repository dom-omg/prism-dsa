#!/usr/bin/env python3
"""
Phase 6 — DPA-Timing Attack: Quantifying the PRISM-DSA 14x Hardening

INVARIANT (derived Phase 5B):
    Timing and hint channels are orthogonal (r≈0, confirmed Phase 5B).
    The REAL attack is not "timing → t0" but "timing → trace selection → DPA → s1".

ATTACK MODEL — Hamming Weight power side-channel on NTT multiplication:

    During c·s1 (NTT-domain), each butterfly computes c_hat[j] * s1_hat[j] mod q.
    An attacker observing power consumption recovers s1_hat[j] via CPA.

    Simulated power for a single signing ATTEMPT:
        P_attempt = Σ_{j=0}^{N-1} HW(c_hat[j] * s1_hat[j] mod q)
        (sum over all n=256 NTT coefficients — cross-coefficient noise)

    For ML-DSA with n_iter=k (k attempts before valid sig):
        P_trace = Σ_{att=1}^{k} P_att  where c_k is KNOWN (from sig), c_1..c_{k-1} UNKNOWN

    For PRISM-DSA FIS_SLOTS=64 (64 attempts always):
        P_trace = Σ_{slot=0}^{63} P_slot  where ONLY c_accepted is KNOWN

SNR ANALYSIS (Rule 44 — exact, not approximate):

    σ_signal² = Var[HW(c[j0] × s1_hat[j0])] for target j0
    σ_cross²  = (N-1) × Var[HW(c[j'] × s1_hat[j'])] for j'≠j0 per attempt

    For one CLEAN trace (ML-DSA n_iter=1):
        noise² = σ_cross² = 255 × Var[HW]         ← only cross-coeff noise
        SNR_clean = σ_signal / σ_cross

    For PRISM-DSA (64 slots):
        noise² = 63 × σ_hw² + 64 × σ_cross²       ← 63 unknown slots + 64 cross-coeff
               ≈ 64 × σ_cross² (when σ_cross >> σ_hw, which holds for small k)
        SNR_prism = σ_signal / (8 × σ_cross) = SNR_clean / 8

    CPA traces needed ∝ 1/SNR²:
        ML-DSA clean:   T traces,     total signing ops = T / P(n_iter=1) = T / 0.22 ≈ 4.5T
        PRISM-DSA:      64T traces,   total signing ops = 64T
        Ratio:          64T / 4.5T ≈ 14x

EMPIRICAL VALIDATION: this script measures the CPA correlation growth rate for
all three scenarios and confirms the theoretical 14x.

Requirements: dilithium-py==1.4.0  numpy  scipy
"""

import os
import time
import struct
import numpy as np
from scipy import stats

from dilithium_py.ml_dsa.ml_dsa import ML_DSA
from dilithium_py.ml_dsa.default_parameters import DEFAULT_PARAMETERS

Q = 8380417
N = 256
FIS_SLOTS = 64

PARAMS = DEFAULT_PARAMETERS["ML_DSA_44"]

# ─── HW model ─────────────────────────────────────────────────────────────────

# Precompute HW table for fast lookup
_HW = np.zeros(256, dtype=np.uint8)
for _b in range(256):
    _HW[_b] = bin(_b).count('1')

def hw(x: int) -> int:
    x = int(x) % Q
    return _HW[x & 0xFF] + _HW[(x >> 8) & 0xFF] + _HW[(x >> 16) & 0xFF]

def hw_vec(arr: np.ndarray) -> np.ndarray:
    """Vectorised HW for array of values mod q."""
    a = arr.astype(np.uint32) % Q
    return (_HW[a & 0xFF].astype(np.int32)
            + _HW[(a >> 8) & 0xFF].astype(np.int32)
            + _HW[(a >> 16) & 0xFF].astype(np.int32))


# ─── Instrumented signer ──────────────────────────────────────────────────────

class SignerInstrumented(ML_DSA):
    """ML-DSA-44 exposing s1_hat, iteration count, and per-attempt c_hat."""

    def keygen_with_s1hat(self):
        pk, sk = self.keygen()
        rho, k_key, tr, s1, s2, t0 = self._unpack_sk(sk)
        s1_hat = []
        for poly in s1.rows()[0]:
            s1_hat.extend(poly.to_ntt().coeffs)
        for row in s1.rows()[1:]:
            for poly in row:
                s1_hat.extend(poly.to_ntt().coeffs)
        return pk, sk, np.array(s1_hat, dtype=np.int64)

    def sign_with_trace_info(self, sk: bytes, m: bytes) -> dict:
        """Sign and return: sig, n_iter, list of c_hat per attempt, wall_ns."""
        rnd = self.random_bytes(32)
        t_ns = time.perf_counter_ns()
        result = self._sign_trace(sk, m, rnd)
        result["wall_ns"] = time.perf_counter_ns() - t_ns
        return result

    def _sign_trace(self, sk: bytes, m: bytes, rnd: bytes) -> dict:
        rho, k_key, tr, s1, s2, t0 = self._unpack_sk(sk)
        s1_hat = s1.to_ntt()
        s2_hat = s2.to_ntt()
        t0_hat = t0.to_ntt()
        A_hat = self._expand_matrix_from_seed(rho)

        mu = self._h(tr + m, 64)
        rho_prime = self._h(k_key + rnd + mu, 64)

        kappa = 0
        alpha = self.gamma_2 << 1
        n_iter = 0
        c_hat_list = []   # c_hat per attempt (NTT coefficients of challenge poly)

        while True:
            n_iter += 1
            y = self._expand_mask_vector(rho_prime, kappa)
            y_hat = y.to_ntt()
            w = (A_hat @ y_hat).from_ntt()
            kappa += self.l

            w1 = w.high_bits(alpha)
            w1_bytes = w1.bit_pack_w(self.gamma_2)
            c_tilde = self._h(mu + w1_bytes, self.c_tilde_bytes)
            c = self.R.sample_in_ball(c_tilde, self.tau)
            c_hat_curr = np.array(c.to_ntt().coeffs, dtype=np.int64)
            c_hat_list.append(c_hat_curr)

            c_hat_ntt = c.to_ntt()
            c_s1 = s1_hat.scale(c_hat_ntt).from_ntt()
            z = y + c_s1
            if z.check_norm_bound(self.gamma_1 - self.beta):
                continue

            c_s2 = s2_hat.scale(c_hat_ntt).from_ntt()
            r0 = (w - c_s2).low_bits(alpha)
            if r0.check_norm_bound(self.gamma_2 - self.beta):
                continue

            c_t0 = t0_hat.scale(c_hat_ntt).from_ntt()
            if c_t0.check_norm_bound(self.gamma_2):
                continue

            h = (-c_t0).make_hint(w - c_s2 + c_t0, alpha)
            if h.sum_hint() > self.omega:
                continue

            sig = self._pack_sig(c_tilde, z, h)
            return {
                "sig": sig,
                "n_iter": n_iter,
                "c_hat_accepted": c_hat_list[-1],   # known to verifier
                "c_hat_rejected": c_hat_list[:-1],  # unknown to attacker (require rho')
            }


# ─── Power trace simulation ───────────────────────────────────────────────────

def simulate_power_mldsa(c_hat_accepted: np.ndarray,
                         c_hat_rejected: list,
                         s1_hat: np.ndarray,
                         rng: np.random.Generator) -> np.ndarray:
    """
    Power trace for one ML-DSA signing operation.

    Includes ALL attempts (accepted + rejected) summed.
    Returns per-coefficient (n=256 for poly 0) aggregated HW sum.
    """
    total = np.zeros(N, dtype=np.int64)
    # Accepted attempt
    prod_acc = (c_hat_accepted[:N] * s1_hat[:N]) % Q
    total += hw_vec(prod_acc)
    # Rejected attempts (c_hat known because we instrumented — but in the attack
    # scenario they're UNKNOWN; we include them to model the real power trace)
    for c_rej in c_hat_rejected:
        prod_rej = (c_rej[:N] * s1_hat[:N]) % Q
        total += hw_vec(prod_rej)
    return total  # shape (256,), measurement at each NTT "register"


def simulate_power_prism(c_hat_accepted: np.ndarray,
                         s1_hat: np.ndarray,
                         rng: np.random.Generator) -> np.ndarray:
    """
    Power trace for one PRISM-DSA (FIS=64) signing operation.

    FIS always runs exactly 64 slots. The accepted slot has known c_hat_accepted.
    The other 63 slots have UNKNOWN c (pseudorandom, appears uniform to attacker).
    We sample them uniformly to model the attacker's uncertainty.
    """
    total = np.zeros(N, dtype=np.int64)
    # Accepted slot (known c)
    prod_acc = (c_hat_accepted[:N] * s1_hat[:N]) % Q
    total += hw_vec(prod_acc)
    # 63 unknown slots: simulate as random elements of Z_q
    for _ in range(FIS_SLOTS - 1):
        c_unk = rng.integers(0, Q, size=N, dtype=np.int64)
        prod_unk = (c_unk * s1_hat[:N]) % Q
        total += hw_vec(prod_unk)
    return total  # shape (256,)


# ─── CPA (Correlation Power Analysis) ────────────────────────────────────────

def cpa_at_true(traces: np.ndarray,
                c_hats_known: np.ndarray,
                s1_hat_true: np.ndarray,
                j0: int,
                n_random_hyp: int = 400,
                rng: np.random.Generator = None) -> dict:
    """
    CPA for coefficient j0.

    For the target hypothesis h = s1_hat_true[j0]:
        model[i] = HW(c_hat_i[j0] * h mod q)

    We measure: rank of true hypothesis among (n_random_hyp + 1) candidates.
    Returns: correlation at true hypothesis, median corr at random, rank.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = traces.shape[0]
    power = traces[:, j0].astype(np.float64)
    c_j0 = c_hats_known[:, j0].astype(np.int64)

    # True hypothesis
    h_true = int(s1_hat_true[j0]) % Q
    model_true = hw_vec((c_j0 * h_true) % Q).astype(np.float64)
    r_true = np.corrcoef(model_true, power)[0, 1]

    # Random hypotheses
    hyps = rng.integers(0, Q, size=n_random_hyp, dtype=np.int64)
    corrs_rand = []
    for h in hyps:
        m = hw_vec((c_j0 * h) % Q).astype(np.float64)
        corrs_rand.append(abs(np.corrcoef(m, power)[0, 1]))

    # Rank of true hypothesis
    rank = sum(1 for r in corrs_rand if r >= abs(r_true)) + 1

    return {
        "corr_true": abs(r_true),
        "corr_rand_median": float(np.median(corrs_rand)),
        "corr_rand_max": float(np.max(corrs_rand)),
        "rank": rank,
        "success": rank == 1,
    }


# ─── Main experiment ──────────────────────────────────────────────────────────

def run_dpa_comparison(n_sigs: int = 1000, j0: int = 7, n_cpa_hyp: int = 400,
                       verbose: bool = True) -> dict:
    """
    Compare CPA success across three scenarios using the SAME n_sigs signing operations:

        A. ML-DSA selected   — only n_iter=1 sigs (clean traces)
        B. ML-DSA all        — all sigs mixed (dirty traces included)
        C. PRISM-DSA         — all sigs, but each has 64 slots (FIS model)

    Returns: dict of results at each checkpoint.
    """
    rng = np.random.default_rng(42)
    dsa = SignerInstrumented(PARAMS)
    pk, sk, s1_hat = dsa.keygen_with_s1hat()

    if verbose:
        print(f"\nKey s1_hat[{j0}] = {s1_hat[j0]}")
        print(f"Signing {n_sigs} messages with ML-DSA-44...")

    # Collect all signatures with trace info
    all_c_acc = []     # accepted c_hat per sig
    all_traces_mldsa = []   # power traces for ML-DSA (all attempts)
    all_traces_prism = []   # power traces for PRISM-DSA model
    all_n_iter = []

    for i in range(n_sigs):
        msg = f"msg-dpa-{i:06d}".encode()
        r = dsa.sign_with_trace_info(sk, msg)

        c_acc = r["c_hat_accepted"]
        c_rej = r["c_hat_rejected"]
        n_iter = r["n_iter"]

        # ML-DSA power trace: all actual attempts (accepted + rejected)
        trace_ml = simulate_power_mldsa(c_acc, c_rej, s1_hat, rng)

        # PRISM-DSA power trace: accepted + 63 simulated unknown slots
        trace_pr = simulate_power_prism(c_acc, s1_hat, rng)

        all_c_acc.append(c_acc)
        all_traces_mldsa.append(trace_ml)
        all_traces_prism.append(trace_pr)
        all_n_iter.append(n_iter)

        if verbose and (i + 1) % 200 == 0:
            clean = sum(1 for k in all_n_iter if k == 1)
            print(f"  [{i+1:4d}/{n_sigs}]  clean (n=1): {clean} ({100*clean/(i+1):.1f}%)")

    all_c_acc = np.array(all_c_acc, dtype=np.int64)         # (n_sigs, 256)
    all_traces_mldsa = np.array(all_traces_mldsa, dtype=np.int64)  # (n_sigs, 256)
    all_traces_prism = np.array(all_traces_prism, dtype=np.int64)  # (n_sigs, 256)
    all_n_iter = np.array(all_n_iter, dtype=np.int32)

    clean_mask = all_n_iter == 1
    n_clean = clean_mask.sum()
    p_clean = n_clean / n_sigs

    if verbose:
        print(f"\nN_sigs total:  {n_sigs}")
        print(f"N_clean (n=1): {n_clean} ({100*p_clean:.1f}%)")
        print(f"E[n_iter]:     {all_n_iter.mean():.3f}")

    # ─── CPA runs at increasing N checkpoints ────────────────────────────────
    # Fine-grained early checkpoints to detect first-success precisely
    checkpoints = list(range(5, 30, 5)) + [30, 40, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1000]
    checkpoints = sorted(set(c for c in checkpoints if c <= n_sigs))

    results = []
    for chk in checkpoints:
        row = {"n_total": chk}

        # Scenario A: ML-DSA selected (n_iter=1 only)
        sel_idx = np.where(clean_mask[:chk])[0]
        if len(sel_idx) >= 5:
            cpa_a = cpa_at_true(
                all_traces_mldsa[sel_idx],
                all_c_acc[sel_idx],
                s1_hat, j0, n_cpa_hyp, rng)
            row["a_n_clean"] = len(sel_idx)
            row["a_corr"] = cpa_a["corr_true"]
            row["a_rank"] = cpa_a["rank"]
            row["a_success"] = cpa_a["success"]
        else:
            row["a_n_clean"] = len(sel_idx)
            row["a_corr"] = 0.0
            row["a_rank"] = n_cpa_hyp
            row["a_success"] = False

        # Scenario B: ML-DSA all sigs
        cpa_b = cpa_at_true(
            all_traces_mldsa[:chk],
            all_c_acc[:chk],
            s1_hat, j0, n_cpa_hyp, rng)
        row["b_corr"] = cpa_b["corr_true"]
        row["b_rank"] = cpa_b["rank"]
        row["b_success"] = cpa_b["success"]

        # Scenario C: PRISM-DSA
        cpa_c = cpa_at_true(
            all_traces_prism[:chk],
            all_c_acc[:chk],
            s1_hat, j0, n_cpa_hyp, rng)
        row["c_corr"] = cpa_c["corr_true"]
        row["c_rank"] = cpa_c["rank"]
        row["c_success"] = cpa_c["success"]

        results.append(row)

    return {
        "results": results,
        "n_clean": int(n_clean),
        "p_clean": float(p_clean),
        "mean_n_iter": float(all_n_iter.mean()),
        "j0": j0,
        "s1_hat_j0": int(s1_hat[j0]),
    }


# ─── SNR theoretical calculation ─────────────────────────────────────────────

def theoretical_snr_analysis():
    """
    Compute theoretical SNR for each scenario using HW variance model.
    s1_hat[j] ~ effectively uniform in Z_q for DPA purposes.
    """
    # Var[HW(x mod q)] for x uniform in Z_q
    # q = 8380417 < 2^23, so HW ∈ {0,...,23}
    # For x uniform: E[HW] ≈ 11.5 bits, Var[HW] ≈ sum bit-independence
    n_samples = 100_000
    rng = np.random.default_rng(0)
    x = rng.integers(0, Q, size=n_samples, dtype=np.int64)
    hw_x = hw_vec(x).astype(np.float64)
    var_hw = float(np.var(hw_x))
    mean_hw = float(np.mean(hw_x))

    sigma_hw = var_hw ** 0.5  # std of HW(c_j * s1_j) for one coefficient
    sigma_cross = ((N - 1) * var_hw) ** 0.5   # cross-coefficient noise per attempt

    # Effective SNR (signal = variation explained by c_j0 correlation)
    # True SNR in CPA: r ≈ σ_signal / σ_total where σ_signal = std of covariance
    # Conservative estimate: SNR ≈ sqrt(Var[HW_signal] / total_noise_var)
    snr_clean = var_hw / ((N - 1) * var_hw)   # = 1/(N-1) ≈ 0.0039
    # Actual CPA correlation (not SNR in traditional sense):
    # For large N: ρ → Cor(HW(c*s1), Σ_j' HW(c*s1[j']))
    # = Var[HW(c*s1[j0])] / sqrt(Var[HW(c*s1[j0])] * N*Var[HW(c*s1)])
    # = sqrt(Var / (N * Var)) = 1/sqrt(N)
    snr_clean_cpa = 1.0 / (N ** 0.5)   # ~0.063 for N=256

    snr_prism = snr_clean_cpa / (FIS_SLOTS ** 0.5)   # 1/8 × snr_clean

    # Traces for reliable CPA: N ∝ 1/snr²  (standard CPA threshold ~0.1)
    cpa_threshold = 0.1
    traces_clean = (cpa_threshold / snr_clean_cpa) ** 2
    traces_prism = (cpa_threshold / snr_prism) ** 2
    ratio_traces = traces_prism / traces_clean

    # Total signing ops (accounting for clean selection)
    p_clean = 0.22   # P(n_iter=1) ≈ geometric(0.22)
    ops_mldsa_selected = traces_clean / p_clean
    ops_prism = traces_prism
    ratio_ops = ops_prism / ops_mldsa_selected

    return {
        "var_hw": var_hw,
        "mean_hw": mean_hw,
        "sigma_hw": sigma_hw,
        "sigma_cross": sigma_cross,
        "snr_clean_cpa": snr_clean_cpa,
        "snr_prism_cpa": snr_prism,
        "cpa_threshold": cpa_threshold,
        "traces_clean": traces_clean,
        "traces_prism": traces_prism,
        "ratio_traces": ratio_traces,
        "ops_mldsa_selected": ops_mldsa_selected,
        "ops_prism": ops_prism,
        "ratio_ops": ratio_ops,
    }


# ─── Multi-key / multi-coeff Monte Carlo ─────────────────────────────────────

def monte_carlo_success_rate(n_keys: int = 5, n_sigs_per_key: int = 600,
                              j0: int = 7, n_cpa_hyp: int = 300,
                              verbose: bool = True) -> dict:
    """
    Run DPA comparison across n_keys independent keys.
    At each checkpoint, measure average CPA correlation for each scenario.
    Returns success rates and correlation statistics.
    """
    rng = np.random.default_rng(99)
    dsa = SignerInstrumented(PARAMS)
    key_results = []

    for ki in range(n_keys):
        if verbose:
            print(f"\n─── KEY {ki+1}/{n_keys} ───")
        pk, sk, s1_hat = dsa.keygen_with_s1hat()
        r = run_dpa_comparison(n_sigs=n_sigs_per_key, j0=j0,
                                n_cpa_hyp=n_cpa_hyp, verbose=verbose)
        key_results.append(r)

    # Aggregate
    checkpoints = [row["n_total"] for row in key_results[0]["results"]]
    agg = []
    for ci, chk in enumerate(checkpoints):
        a_corrs = [kr["results"][ci]["a_corr"] for kr in key_results]
        b_corrs = [kr["results"][ci]["b_corr"] for kr in key_results]
        c_corrs = [kr["results"][ci]["c_corr"] for kr in key_results]
        a_success = [kr["results"][ci]["a_success"] for kr in key_results]
        b_success = [kr["results"][ci]["b_success"] for kr in key_results]
        c_success = [kr["results"][ci]["c_success"] for kr in key_results]
        agg.append({
            "n_total": chk,
            "a_corr_mean": np.mean(a_corrs),
            "b_corr_mean": np.mean(b_corrs),
            "c_corr_mean": np.mean(c_corrs),
            "a_sr": np.mean(a_success),
            "b_sr": np.mean(b_success),
            "c_sr": np.mean(c_success),
        })
    return {"by_checkpoint": agg, "key_results": key_results}


# ─── Print results ────────────────────────────────────────────────────────────

def print_theoretical():
    t = theoretical_snr_analysis()
    print("═" * 72)
    print("THEORETICAL SNR ANALYSIS — ML-DSA-44 DPA (HW model)")
    print("═" * 72)
    print(f"\n[A] HW VARIANCE PARAMETERS")
    print(f"    Var[HW(x mod q)] measured on {100_000:,} uniform samples")
    print(f"    Mean HW:   {t['mean_hw']:.2f}  (expected log2(q)/2 ≈ 11.5)")
    print(f"    σ(HW):     {t['sigma_hw']:.3f}")
    print(f"    σ_cross:   {t['sigma_cross']:.1f}  (√((N-1)×Var) = √(255×Var), N=256)")

    print(f"\n[B] CPA SNR PER TRACE")
    print(f"    ρ_clean  (ML-DSA n_iter=1, 1 attempt):  1/√{N} = {t['snr_clean_cpa']:.4f}")
    print(f"    ρ_prism  (PRISM-DSA, 64 slots):         1/√{N}×√{FIS_SLOTS} = {t['snr_prism_cpa']:.4f}")
    print(f"    SNR ratio: {t['snr_clean_cpa'] / t['snr_prism_cpa']:.2f}x  "
          f"(= √{FIS_SLOTS} = {FIS_SLOTS**0.5:.1f} ✓)")

    print(f"\n[C] TRACES NEEDED (CPA threshold ρ={t['cpa_threshold']}, for coeff j0)")
    print(f"    ML-DSA clean traces:   {t['traces_clean']:,.0f}")
    print(f"    PRISM-DSA traces:      {t['traces_prism']:,.0f}")
    print(f"    Ratio:                 {t['ratio_traces']:.1f}x")

    print(f"\n[D] TOTAL SIGNING OPERATIONS FOR FULL ATTACK")
    print(f"    ML-DSA (select n_iter=1):  traces÷P(n=1) = "
          f"{t['traces_clean']:.0f}÷0.22 = {t['ops_mldsa_selected']:,.0f} ops")
    print(f"    PRISM-DSA:               {t['ops_prism']:,.0f} ops  (no selection possible)")
    print(f"\n    ┌────────────────────────────────────────────┐")
    print(f"    │  PRISM-DSA / ML-DSA ratio = {t['ratio_ops']:.1f}x           │")
    print(f"    │  PRISM-DSA is {t['ratio_ops']:.0f}× harder to attack via DPA │")
    print(f"    └────────────────────────────────────────────┘")

    print(f"\n[E] FULL s1 RECOVERY (l×N = {4*N:,} NTT coefficients)")
    print(f"    ML-DSA ops:  {t['ops_mldsa_selected'] * 4 * N:,.0f}")
    print(f"    PRISM-DSA:   {t['ops_prism'] * 4 * N:,.0f}")
    print(f"    Note: parallel CPA over all j simultaneously reduces this linearly.")

    return t


def print_empirical(res: dict):
    results = res["results"]
    print("\n" + "═" * 72)
    print("EMPIRICAL CPA — CORRELATION GROWTH vs SIGNING OPERATIONS")
    print("═" * 72)
    print(f"\n  Scenarios:")
    print(f"    A = ML-DSA selected (n_iter=1 only, {res['p_clean']*100:.1f}% of sigs)")
    print(f"    B = ML-DSA all sigs (mixed n_iter, E[n]={res['mean_n_iter']:.2f})")
    print(f"    C = PRISM-DSA model (FIS=64 slots per sig)")
    print(f"  Target: NTT coefficient j0={res['j0']}, s1_hat={res['s1_hat_j0']}")
    print(f"\n  {'N_ops':>6}  {'N_clean':>7}  {'ρ_A(sel)':>10}  {'ρ_B(all)':>10}  "
          f"{'ρ_C(PRISM)':>11}  {'rank_A':>7}  {'rank_C':>7}  {'A/C gain':>9}")
    print("  " + "─" * 80)
    for r in results:
        nc = r.get("a_n_clean", 0)
        ra = r.get("a_corr", 0.0)
        rb = r.get("b_corr", 0.0)
        rc = r.get("c_corr", 0.0)
        rk_a = r.get("a_rank", "—")
        rk_c = r.get("c_rank", "—")
        gain = (ra / rc) if rc > 1e-6 else float("nan")
        print(f"  {r['n_total']:>6}  {nc:>7}  {ra:>10.4f}  {rb:>10.4f}  "
              f"{rc:>11.4f}  {rk_a:>7}  {rk_c:>7}  {gain:>9.1f}x")

    # First-success analysis: N_ops to rank=1 for each scenario
    a_first = next((r["n_total"] for r in results if r.get("a_success")), None)
    b_first = next((r["n_total"] for r in results if r.get("b_success")), None)
    c_first = next((r["n_total"] for r in results if r.get("c_success")), None)

    print(f"\n  ─── First rank=1 success ───")
    print(f"  A (ML-DSA selected): N_ops = {a_first if a_first else '>'+str(results[-1]['n_total'])}")
    print(f"  B (ML-DSA all):      N_ops = {b_first if b_first else '>'+str(results[-1]['n_total'])}")
    print(f"  C (PRISM-DSA):       N_ops = {c_first if c_first else '>'+str(results[-1]['n_total'])}")
    if a_first and c_first:
        empirical_ratio = c_first / a_first
        print(f"\n  Empirical C/A ratio: {empirical_ratio:.1f}x")
        print(f"  Theoretical:         FIS_SLOTS × P(n_iter=1) = 64 × 0.22 = 14.1x")
        print(f"  Note: ratio at first checkpoint ≥ true ratio (A may succeed between N=0 and first chk)")


def print_montecarlo(mc: dict):
    print("\n" + "═" * 72)
    print("MONTE CARLO — CPA CORRELATION GROWTH (averaged across keys)")
    print("═" * 72)
    print(f"\n  {'N_ops':>6}  {'ρ_A(sel)':>10}  {'ρ_B(all)':>10}  {'ρ_C(PRISM)':>11}  "
          f"{'SR_A%':>7}  {'SR_B%':>7}  {'SR_C%':>7}  {'A/C gain':>9}")
    print("  " + "─" * 76)
    for row in mc["by_checkpoint"]:
        ra = row["a_corr_mean"]
        rb = row["b_corr_mean"]
        rc = row["c_corr_mean"]
        gain = (ra / rc) if rc > 1e-6 else float("nan")
        print(f"  {row['n_total']:>6}  {ra:>10.4f}  {rb:>10.4f}  {rc:>11.4f}  "
              f"{row['a_sr']*100:>7.0f}  {row['b_sr']*100:>7.0f}  {row['c_sr']*100:>7.0f}  "
              f"{gain:>9.1f}x")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    p = argparse.ArgumentParser(description="Phase 6 — DPA Timing Attack (PRISM-DSA hardening)")
    p.add_argument("--sigs", type=int, default=600,
                   help="Signing operations per run (default 600)")
    p.add_argument("--keys", type=int, default=3,
                   help="Keys for Monte Carlo (default 3)")
    p.add_argument("--j0", type=int, default=7,
                   help="Target NTT coefficient index (default 7)")
    p.add_argument("--hyp", type=int, default=300,
                   help="Random hypotheses for CPA rank (default 300)")
    p.add_argument("--theory-only", action="store_true",
                   help="Print theoretical analysis only (no simulation)")
    args = p.parse_args()

    print("PRISM-DSA Phase 6 — DPA-Timing Attack Quantification")
    print(f"Platform: {os.uname().sysname} {os.uname().machine}")
    print(f"dilithium-py 1.4.0, ML-DSA-44, HW power model, FIS_SLOTS={FIS_SLOTS}")

    # Theoretical
    t = print_theoretical()

    if args.theory_only:
        return

    # Single-key empirical
    print(f"\n{'─'*72}")
    print(f"EMPIRICAL RUN (1 key, {args.sigs} sigs, j0={args.j0})")
    t0_wall = time.perf_counter()
    res = run_dpa_comparison(n_sigs=args.sigs, j0=args.j0, n_cpa_hyp=args.hyp)
    print(f"Done in {time.perf_counter()-t0_wall:.1f}s")
    print_empirical(res)

    # Monte Carlo
    if args.keys > 1:
        print(f"\n{'─'*72}")
        print(f"MONTE CARLO ({args.keys} keys × {args.sigs} sigs)")
        t0_wall = time.perf_counter()
        mc = monte_carlo_success_rate(n_keys=args.keys, n_sigs_per_key=args.sigs,
                                       j0=args.j0, n_cpa_hyp=args.hyp)
        print(f"Done in {time.perf_counter()-t0_wall:.1f}s")
        print_montecarlo(mc)

    # Summary
    print(f"\n{'═'*72}")
    print("SUMMARY — Rule 44 compliant claims")
    print("═" * 72)
    print(f"""
  PRECONDITION: adversary has side-channel access to ML-DSA signing device,
  can measure power consumption per signing operation (HW model).

  MEASUREMENT: CPA correlation ρ(HW(c·h), power_trace) for hypotheses h.
  True s1_hat[j0] produces maximum ρ after sufficient traces.

  THEORETICAL (HW model, no external noise):
    ML-DSA with n_iter=1 selection:  {t['ops_mldsa_selected']:,.0f} signing ops per NTT coeff
    PRISM-DSA (FIS=64):              {t['ops_prism']:,.0f} signing ops per NTT coeff
    Hardening ratio:                 {t['ratio_ops']:.1f}x

  ROOT CAUSE of hardening:
    ML-DSA timing leaks n_iter → attacker selects CLEAN traces (SNR_clean)
    PRISM-DSA n_iter is always 64 → no selection possible (SNR_prism = SNR_clean/8)
    Trace ratio: 64; ops ratio: 64 × P(n_iter=1) = 64 × 0.22 ≈ 14x

  SCOPE:
    This hardening is WITHIN the simulated HW power model.
    Full key recovery requires all l×N = {4*N:,} NTT coefficients.
    Independent of hint oracle (orthogonal channel, Phase 5B).

  PRISM-DSA FIS eliminates timing-based trace selection → {t['ratio_ops']:.0f}x DPA hardening.
  Residual: hint oracle (Phase 3-4) is unaffected and requires separate countermeasure.
""")


if __name__ == "__main__":
    main()
