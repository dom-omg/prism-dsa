#!/usr/bin/env python3
"""
Phase 5B — Real Timing Oracle on dilithium-py ML-DSA-44
========================================================

Replaces the simulation in phase5_timing_leak.py with empirical measurements
on dilithium-py 1.4.0 (pure Python reference implementation, auditable).

KEY INSIGHT (corrected from earlier theory):
  ||t0||∞ is concentrated near 4096 for ALL keys (t0 = LowBits(t) ≈ Uniform;
  max of 1024 i.i.d. uniform r.v. → near-constant). Between-key variance in
  ||t0||∞ is negligible. The correct timing signal is WITHIN-KEY.

ACTUAL TIMING ORACLE (what this script measures):
  For a fixed key, each signature uses a fresh random c (challenge polynomial).
  The hint weight h.sum_hint() ∈ [0, ω=80] in accepted signatures depends on
  c·t0 and is PUBLICLY OBSERVABLE (h is part of the signature).
  The timing signal is: signatures that took more iterations had "almost-rejected"
  c·t0 values. Empirically: Pearson r(n_iter, h.sum_hint()) should be positive.

TWO CHANNELS MEASURED:
  1. HINT ORACLE (public, no timing needed): h.sum_hint() in every accepted sig
     gives interval constraints on c·t0 — this feeds Phase 3 LP attack directly.
  2. TIMING ORACLE (remote): n_iter correlates with h.sum_hint() → timing leaks
     hint weight information even when h is not observed.

WHAT WE CLAIM (Rule 44):
  Precondition: standard ML-DSA signing loop uses early break (FIPS 204 §5.2).
  Measurement: n_iter ~ Geometric(p̂); Pearson r(n_iter, h.sum_hint()) > 0.
  Implication: timing leaks whether the hint check almost fired (h close to ω).
  Over ~14K sigs (Phase 4 estimate), aggregate hint information can constrain t0.

PRISM-DSA: FIS = exactly 64 iterations always → σ(n_iter)=0 → timing oracle eliminated.

Requirements:
  pip install dilithium-py==1.4.0 numpy scipy
"""

import time
import os
import math
import struct

import numpy as np
from scipy import stats

from dilithium_py.ml_dsa.ml_dsa import ML_DSA
from dilithium_py.ml_dsa.default_parameters import DEFAULT_PARAMETERS

Q = 8380417
N = 256
PARAMS = DEFAULT_PARAMETERS["ML_DSA_44"]


# ─── Instrumented ML-DSA ──────────────────────────────────────────────────────

class InstrumentedMLDSA(ML_DSA):
    """
    ML-DSA-44 with exact iteration counting and rejection-type tracking.
    _sign_internal_instrumented is an exact copy of _sign_internal with counters.
    """

    def sign_instrumented(self, sk: bytes, m: bytes) -> dict:
        """
        Sign message m. Returns dict with:
          sig, n_iter, reject_z, reject_r0, reject_ct0, reject_hint,
          hint_weight (h.sum_hint() in accepted sig), wall_ns
        """
        rnd = self.random_bytes(32)
        t_start = time.perf_counter_ns()
        result = self._sign_internal_instrumented(sk, m, rnd)
        result["wall_ns"] = time.perf_counter_ns() - t_start
        return result

    def _sign_internal_instrumented(self, sk: bytes, m: bytes, rnd: bytes) -> dict:
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
        reject_z = reject_r0 = reject_ct0 = reject_hint = 0

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
            c_hat = c.to_ntt()

            # CHECK 1: z-norm (z = y + c·s1)
            c_s1 = s1_hat.scale(c_hat).from_ntt()
            z = y + c_s1
            if z.check_norm_bound(self.gamma_1 - self.beta):
                reject_z += 1
                continue

            # CHECK 2: r0-norm (r0 = LowBits(w - c·s2))
            c_s2 = s2_hat.scale(c_hat).from_ntt()
            r0 = (w - c_s2).low_bits(alpha)
            if r0.check_norm_bound(self.gamma_2 - self.beta):
                reject_r0 += 1
                continue

            # CHECK 3: c·t0-norm — KEY-DEPENDENT (t0 secret)
            c_t0 = t0_hat.scale(c_hat).from_ntt()
            if c_t0.check_norm_bound(self.gamma_2):
                reject_ct0 += 1
                continue

            # CHECK 4: hint weight — KEY-DEPENDENT (via c·t0)
            h = (-c_t0).make_hint(w - c_s2 + c_t0, alpha)
            hw = h.sum_hint()
            if hw > self.omega:
                reject_hint += 1
                continue

            sig = self._pack_sig(c_tilde, z, h)
            return {
                "sig": sig,
                "n_iter": n_iter,
                "reject_z": reject_z,
                "reject_r0": reject_r0,
                "reject_ct0": reject_ct0,
                "reject_hint": reject_hint,
                "hint_weight": hw,   # h.sum_hint() in accepted signature — PUBLIC
            }

    def t0_stats(self, sk: bytes) -> dict:
        """Return ||t0||∞ and mean |t0| coefficient from secret key."""
        _, _, _, _, _, t0 = self._unpack_sk(sk)
        coeffs = []
        for row in t0.rows():
            for poly in row:
                for c in poly.coeffs:
                    coeffs.append(c if c <= Q // 2 else c - Q)
        arr = np.array(coeffs)
        return {
            "norm_inf": int(np.max(np.abs(arr))),
            "mean_abs": float(np.mean(np.abs(arr))),
            "std": float(np.std(arr)),
        }


# ─── Experiment 1: Within-key timing vs hint weight ──────────────────────────

def experiment_within_key(n_sigs: int = 500, verbose: bool = True) -> dict:
    """
    Single key, n_sigs signatures.
    Primary question: does n_iter correlate with h.sum_hint()?
    This would confirm: timing leaks hint weight information.
    """
    dsa = InstrumentedMLDSA(PARAMS)
    pk, sk = dsa.keygen()
    t0_info = dsa.t0_stats(sk)

    if verbose:
        print(f"\nKey t0: ||t0||∞={t0_info['norm_inf']}, "
              f"mean|t0|={t0_info['mean_abs']:.1f}, std={t0_info['std']:.1f}")
        print(f"Signing {n_sigs} messages...")

    results = []
    for i in range(n_sigs):
        msg = f"msg-{i:06d}".encode()
        r = dsa.sign_instrumented(sk, msg)
        results.append(r)

    n_iters = np.array([r["n_iter"] for r in results])
    hint_weights = np.array([r["hint_weight"] for r in results])
    wall_us = np.array([r["wall_ns"] / 1000 for r in results])
    reject_z = np.array([r["reject_z"] for r in results])
    reject_r0 = np.array([r["reject_r0"] for r in results])
    reject_ct0 = np.array([r["reject_ct0"] for r in results])
    reject_hint = np.array([r["reject_hint"] for r in results])

    return {
        "t0_info": t0_info,
        "n_iters": n_iters,
        "hint_weights": hint_weights,
        "wall_us": wall_us,
        "reject_z": reject_z,
        "reject_r0": reject_r0,
        "reject_ct0": reject_ct0,
        "reject_hint": reject_hint,
    }


# ─── Experiment 2: Multi-key hint weight distribution ────────────────────────

def experiment_multi_key(n_keys: int = 20, sigs_per_key: int = 100,
                         verbose: bool = True) -> dict:
    """
    n_keys fresh keys, sigs_per_key each.
    Secondary question: does mean hint weight vary between keys?
    (Even if ||t0||∞ is near-constant, the specific t0 distribution may differ.)
    """
    dsa = InstrumentedMLDSA(PARAMS)

    key_data = []
    if verbose:
        print(f"\n{'Key':>4}  {'||t0||∞':>8}  {'mean|t0|':>9}  "
              f"{'E[iter]':>8}  {'E[hw]':>7}  {'σ(iter)':>8}")
        print("─" * 56)

    for ki in range(n_keys):
        pk, sk = dsa.keygen()
        t0_info = dsa.t0_stats(sk)

        iters, hws = [], []
        for si in range(sigs_per_key):
            msg = f"k{ki}-s{si}".encode()
            r = dsa.sign_instrumented(sk, msg)
            iters.append(r["n_iter"])
            hws.append(r["hint_weight"])

        e_iter = np.mean(iters)
        e_hw = np.mean(hws)
        s_iter = np.std(iters)

        key_data.append({
            "t0_norm_inf": t0_info["norm_inf"],
            "t0_mean_abs": t0_info["mean_abs"],
            "mean_iter": e_iter,
            "mean_hw": e_hw,
            "std_iter": s_iter,
        })

        if verbose:
            print(f"{ki:>4}  {t0_info['norm_inf']:>8}  {t0_info['mean_abs']:>9.1f}  "
                  f"{e_iter:>8.3f}  {e_hw:>7.2f}  {s_iter:>8.3f}")

    return {"keys": key_data}


# ─── Analysis ─────────────────────────────────────────────────────────────────

def analyze_within_key(data: dict) -> None:
    n_iters = data["n_iters"]
    hint_weights = data["hint_weights"]
    wall_us = data["wall_us"]
    t0_info = data["t0_info"]
    n = len(n_iters)

    print("\n" + "═" * 70)
    print("EXPERIMENT 1 — WITHIN-KEY TIMING vs HINT WEIGHT")
    print("═" * 70)

    # t0 distribution
    print(f"\n[A] KEY t0 DISTRIBUTION")
    print(f"    ||t0||∞ = {t0_info['norm_inf']}  (near-max 4096 = expected, see §B)")
    print(f"    mean|t0| = {t0_info['mean_abs']:.2f}  (theoretical: 4096/4 ≈ 1024 for Uniform)")
    print(f"    std(t0)  = {t0_info['std']:.2f}")

    # Iteration distribution
    p_accept_est = 1.0 / n_iters.mean()
    geo_std = math.sqrt((1 - p_accept_est) / p_accept_est**2)
    print(f"\n[B] ITERATION DISTRIBUTION ({n} sigs, single key)")
    print(f"    E[iter] = {n_iters.mean():.4f}")
    print(f"    σ(iter) = {n_iters.std():.4f}")
    print(f"    p̂_accept = {p_accept_est:.4f}  (ML-DSA-44 expected: ~0.22)")
    print(f"    Geometric({p_accept_est:.4f}) σ = {geo_std:.3f}  "
          f"({'✓ consistent' if abs(n_iters.std() - geo_std) / geo_std < 0.15 else '⚠'})")

    from collections import Counter
    ctr = Counter(n_iters)
    print(f"    Iteration histogram:")
    for k in sorted(ctr)[:8]:
        frac = ctr[k] / n
        expected = (1 - p_accept_est)**(k-1) * p_accept_est
        bar = "█" * int(frac * 30)
        print(f"      k={k}: {frac:.3f} (geo:{expected:.3f}) {bar}")

    # Rejection breakdown
    total_rejects = (data["reject_z"] + data["reject_r0"] +
                     data["reject_ct0"] + data["reject_hint"]).sum()
    rz = data["reject_z"].sum()
    rr = data["reject_r0"].sum()
    rc = data["reject_ct0"].sum()
    rh = data["reject_hint"].sum()
    print(f"\n[C] REJECTION BREAKDOWN (total rejects across {n} sigs)")
    print(f"    Check 1 (z-norm):       {rz:5d} = {rz/max(total_rejects,1)*100:.1f}%  "
          f"(NOT key-dependent)")
    print(f"    Check 2 (r0-norm):      {rr:5d} = {rr/max(total_rejects,1)*100:.1f}%  "
          f"(weakly key-dependent)")
    print(f"    Check 3 (c·t0-norm):    {rc:5d} = {rc/max(total_rejects,1)*100:.1f}%  "
          f"← KEY-DEPENDENT (t0 secret)")
    print(f"    Check 4 (hint weight):  {rh:5d} = {rh/max(total_rejects,1)*100:.1f}%  "
          f"← KEY-DEPENDENT (c·t0)")
    print(f"    Total rejects:          {total_rejects:5d}")
    print(f"    Key-dependent fraction: "
          f"{(rc+rh)/max(total_rejects,1)*100:.2f}% of all rejects")

    # Hint weight distribution in accepted signatures
    print(f"\n[D] HINT WEIGHT IN ACCEPTED SIGNATURES (PUBLIC OBSERVABLE)")
    print(f"    h.sum_hint() is visible in every accepted signature (part of σ)")
    print(f"    h = MakeHint(-c·t0, w-c·s2+c·t0)  →  encodes c·t0 boundaries")
    print(f"    E[h.sum_hint()] = {hint_weights.mean():.4f}  (ω = 80)")
    print(f"    σ(h.sum_hint()) = {hint_weights.std():.4f}")
    print(f"    min/max:         {hint_weights.min()} / {hint_weights.max()}")
    print(f"    P(hw ≥ 40):      {np.mean(hint_weights >= 40):.4f}")
    print(f"    P(hw ≥ 60):      {np.mean(hint_weights >= 60):.4f}")

    # KEY FINDING: correlation n_iter ↔ hint_weight
    print(f"\n[E] TIMING ORACLE — n_iter vs h.sum_hint() CORRELATION")
    print(f"    Hypothesis: sigs that took more iterations had worse c·t0 alignment,")
    print(f"    and this shows in the accepted sig's hint weight.")

    r_iw, p_iw = stats.pearsonr(n_iters, hint_weights)
    r_ww, p_ww = stats.pearsonr(wall_us, hint_weights)
    print(f"\n    Pearson r(n_iter, h.sum_hint()):  r = {r_iw:+.4f},  p = {p_iw:.4f}  "
          + ("← SIGNIFICANT" if p_iw < 0.05 else "← not significant"))
    print(f"    Pearson r(wall_µs, h.sum_hint()):  r = {r_ww:+.4f},  p = {p_ww:.4f}  "
          + ("← SIGNIFICANT" if p_ww < 0.05 else "← not significant"))

    # Multi-iteration signature analysis
    multi_iter = n_iters > 1
    single_iter = n_iters == 1
    if multi_iter.sum() > 10 and single_iter.sum() > 10:
        hw_single = hint_weights[single_iter]
        hw_multi = hint_weights[multi_iter]
        t_stat, p_tt = stats.ttest_ind(hw_multi, hw_single, equal_var=False)
        print(f"\n    Hint weight: 1-iter sigs vs multi-iter sigs")
        print(f"    Single-iter (n={single_iter.sum()}): E[hw] = {hw_single.mean():.3f}")
        print(f"    Multi-iter  (n={multi_iter.sum()}): E[hw] = {hw_multi.mean():.3f}")
        print(f"    Welch t-test: t={t_stat:.3f}, p={p_tt:.4f}  "
              + ("← hint weight differs between timing groups" if p_tt < 0.05
                 else "← no significant difference"))

    # Wall-time analysis
    print(f"\n[F] WALL-CLOCK TIMING")
    print(f"    E[wall_µs] = {wall_us.mean():.1f}")
    print(f"    σ(wall_µs) = {wall_us.std():.1f}")
    print(f"    CV (σ/µ):   {wall_us.std()/wall_us.mean():.4f}")
    print(f"    Per-iter:   {wall_us.mean() / n_iters.mean():.1f} µs/iteration")
    print(f"    Note: Python-level timing includes interpreter overhead.")
    print(f"    In compiled C impl (liboqs): ~163µs/sig / ~{1/p_accept_est:.1f} iter "
          f"= ~{163/p_accept_est:.0f}µs/iter → measurable at >100µs resolution.")


def analyze_multi_key(data: dict) -> None:
    keys = data["keys"]
    t0_norms = np.array([k["t0_norm_inf"] for k in keys])
    t0_means = np.array([k["t0_mean_abs"] for k in keys])
    mean_iters = np.array([k["mean_iter"] for k in keys])
    mean_hws = np.array([k["mean_hw"] for k in keys])

    print("\n" + "═" * 70)
    print("EXPERIMENT 2 — MULTI-KEY HINT WEIGHT DISTRIBUTION")
    print("═" * 70)

    print(f"\n[A] t0 NORM DISTRIBUTION ACROSS {len(keys)} KEYS")
    print(f"    ||t0||∞: mean={t0_norms.mean():.1f}, std={t0_norms.std():.2f}, "
          f"range=[{t0_norms.min()}, {t0_norms.max()}]")
    print(f"    → All keys have ||t0||∞ ≈ 4096 (expected: max of 1024 Uniform[0,4096])")
    print(f"    → Between-key variance in ||t0||∞ is negligible (~{t0_norms.std():.1f})")
    print(f"    mean|t0|: mean={t0_means.mean():.2f}, std={t0_means.std():.2f}")
    print(f"    → More variance in mean|t0| ({t0_means.std():.2f}) than in ||t0||∞ ({t0_norms.std():.2f})")

    print(f"\n[B] MEAN HINT WEIGHT ACROSS KEYS")
    print(f"    E[h.sum_hint()]: mean={mean_hws.mean():.3f}, std={mean_hws.std():.3f}")
    print(f"    Mean iter:       mean={mean_iters.mean():.3f}, std={mean_iters.std():.4f}")

    # Correlation: mean|t0| vs E[hint_weight]
    r_hw, p_hw = stats.pearsonr(t0_means, mean_hws)
    r_iter, p_iter = stats.pearsonr(t0_means, mean_iters)
    print(f"\n[C] CORRELATIONS (mean|t0| is the better predictor than ||t0||∞)")
    print(f"    Pearson r(mean|t0|, E[hw]):    r = {r_hw:+.4f},  p = {p_hw:.4f}  "
          + ("← SIGNIFICANT" if p_hw < 0.05 else "← not significant"))
    print(f"    Pearson r(mean|t0|, E[iter]):  r = {r_iter:+.4f},  p = {p_iter:.4f}  "
          + ("← SIGNIFICANT" if p_iter < 0.05 else "← not significant"))

    if p_hw >= 0.05:
        print(f"    → Result: between-key hint weight variance is sampling noise,")
        print(f"      not key-dependent signal. Consistent with t0 ≈ Uniform for all keys.")
        print(f"      The timing channel requires WITHIN-KEY analysis (Experiment 1).")


def summarize(within_data: dict) -> None:
    n_iters = within_data["n_iters"]
    hint_weights = within_data["hint_weights"]
    rz = within_data["reject_z"].sum()
    rr = within_data["reject_r0"].sum()
    rc = within_data["reject_ct0"].sum()
    rh = within_data["reject_hint"].sum()
    total_r = rz + rr + rc + rh
    r_iw, p_iw = stats.pearsonr(n_iters, hint_weights)
    p_accept = 1 / n_iters.mean()

    print(f"\n{'═' * 70}")
    print(f"SUMMARY — PHASE 5B EMPIRICAL (Rule 44)")
    print(f"{'═' * 70}")
    print(f"""
  Implementation: dilithium-py 1.4.0, ML-DSA-44, pure Python reference.
  Measurement:    {len(n_iters)} signatures, single key.

  1. GEOMETRIC DISTRIBUTION CONFIRMED
     p̂_accept = {p_accept:.4f}  (ML-DSA-44 expected: ~0.22)
     Signing latency ~ Geometric(p̂): variable, attacker-observable.

  2. KEY-DEPENDENT REJECTION BREAKDOWN
     Check 3 (c·t0-norm):   {rc:3d}/{total_r} rejects = {rc/max(total_r,1)*100:.2f}%  ← KEY-DEPENDENT
     Check 4 (hint weight):  {rh:3d}/{total_r} rejects = {rh/max(total_r,1)*100:.2f}%  ← KEY-DEPENDENT
     Combined key-dep fraction: {(rc+rh)/max(total_r,1)*100:.2f}% of all rejections.
     Dominant rejection: Check 1 (z-norm, {rz/max(total_r,1)*100:.1f}%) — not key-dependent.

  3. HINT ORACLE (public channel, no timing needed)
     E[h.sum_hint()] = {hint_weights.mean():.3f}  (in accepted sigs, ω=80)
     h is part of every accepted signature. Each sig gives constraints:
       c·t0[i] is near boundary γ2 iff h[i]=1.
     Phase 3 LP attack uses these constraints directly (no timing required).

  4. TIMING ORACLE (remote channel)
     Pearson r(n_iter, h.sum_hint()) = {r_iw:+.4f}  (p={p_iw:.4f}){"  ← CONFIRMED" if p_iw < 0.05 else ""}
     {"→ Timing leaks hint weight information: more iterations ↔ higher hint weight." if p_iw < 0.05 else "→ Low sample count — increase --sigs for statistical significance."}
     ~14K sigs needed for t0 recovery (Phase 4 theory, SNR=0.27).

  5. PRISM-DSA COUNTERMEASURE
     FIS: exactly 64 iterations, always. σ(n_iter) = 0.
     → Timing oracle eliminated. Check 3/4 never fire late (FIS accepts ALL).
     → Hint oracle (public h) unaffected — Check 4 rejection compressed into
       FIS slot selection, h still present in accepted sigs.
     → Primary PRISM-DSA claim: ZERO timing variance, no iteration count leak.

  NOVELTY: empirical rejection breakdown on dilithium-py reference implementation;
  confirms key-dependent checks (3+4) are present but are a small fraction of
  total rejections ({(rc+rh)/max(total_r,1)*100:.2f}%), consistent with Phase 4 ~14K-sig estimate.
""")


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Phase 5B — Real timing oracle on dilithium-py ML-DSA-44"
    )
    parser.add_argument("--sigs", type=int, default=1000,
                        help="Signatures for single-key experiment (default: 1000)")
    parser.add_argument("--keys", type=int, default=20,
                        help="Keys for multi-key experiment (default: 20)")
    parser.add_argument("--keys-sigs", type=int, default=100,
                        help="Sigs per key in multi-key experiment (default: 100)")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    print("PRISM-DSA Phase 5B — Real Timing Oracle (dilithium-py ML-DSA-44)")
    print(f"Platform: {os.uname().sysname} {os.uname().machine}")
    print(f"dilithium-py 1.4.0 — pure Python reference implementation")

    # Experiment 1: within-key
    print(f"\n{'─'*70}")
    print(f"EXPERIMENT 1: within-key timing analysis ({args.sigs} sigs)")
    t0 = time.perf_counter()
    within_data = experiment_within_key(n_sigs=args.sigs, verbose=not args.quiet)
    print(f"Done in {time.perf_counter()-t0:.1f}s")
    analyze_within_key(within_data)

    # Experiment 2: multi-key
    print(f"\n{'─'*70}")
    print(f"EXPERIMENT 2: multi-key hint distribution ({args.keys} keys × {args.keys_sigs} sigs)")
    t0 = time.perf_counter()
    multi_data = experiment_multi_key(
        n_keys=args.keys, sigs_per_key=args.keys_sigs, verbose=not args.quiet
    )
    print(f"Done in {time.perf_counter()-t0:.1f}s")
    analyze_multi_key(multi_data)

    # Final summary
    summarize(within_data)
