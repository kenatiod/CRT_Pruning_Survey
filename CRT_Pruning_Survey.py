#!/usr/bin/env python3
"""
CRT_Pruning_Survey.py — version 3
================================
Pure CRT partition survey for large omega — no Pell solver required.

For each omega, enumerate all valid support partitions A ⊔ B = P_omega
(bit i of mask = 1 means p_i | m+1, 0 means p_i | m) and compute the CRT
minimum representative n_{A,B} of each class mod M = primorial(omega).

Version 3 changes over version 2
---------------------------------
1. IN-WORKER CENSUS CHECK (the "one candidate per class" conversion).
   Because the pre-flight check enforces B < M, each surviving class
   contains AT MOST ONE integer <= B — its minimal representative.
   Therefore any prime-complete pair with m <= B at ANY level >= omega
   must have m exactly equal to one of the survivor values (a pc pair
   at level omega' >= omega has all primes <= p_omega dividing m(m+1),
   so it lies in a surveyed class, and within [1, B] the class holds
   only its minimum).  V3 checks every survivor, in the worker, for the
   INTRINSIC prime-completeness of m(m+1):
       * trial-divide m(m+1) by the first K primes, where K is minimal
         with primorial(K) > B*(B+1).  If a cofactor remains, then
         GPF(m(m+1)) > p_K, and prime-completeness would force
         m(m+1) >= primorial(K+1) > B*(B+1) >= m(m+1) — impossible.
       * if the cofactor is 1, prime-completeness is decided exactly:
         the set of dividing prime indices must be {1, ..., g} where
         p_g is the GPF.
   This yields the CENSUS CERTIFICATE: "no prime-complete pair with
   m <= B at any level >= omega", valid even when raw survivors exist,
   and INDEPENDENT of all storage caps (every survivor is checked at
   the moment it is found, whether or not it is stored).

2. HEREDITY made explicit.  Level-(omega+1) congruence systems refine
   level-omega systems (forget the new prime's condition), so the floor
   F(omega) = min over classes of the minimal representative is
   nondecreasing in omega, and both certificates are hereditary upward:
   a certificate at omega_0 covers every omega >= omega_0 for the same
   bound B.  The JSON now states this and both certificate statements
   verbatim.

3. BUG FIX: in v2 the CLI flags --survivor_cap_per_task and
   --survivor_cap_total were silently inert (main() assigned them to
   local variables; the bare module-level `global` statement does
   nothing).  V3 passes the caps explicitly as parameters end to end.

4. BUG FIX: in v2 the certificate did not require zero worker errors
   or full task completion; a crashed worker could in principle yield
   a false "all pruned".  V3 certificates require errors == 0 AND
   tasks_done == tasks_total.

5. PERFORMANCE:
   * Worker pool initializer carries the per-omega constants; tasks
     are now slim (high_mask, high_contrib) tuples.
   * Default chunk_size is 1 and H is sized so each task is about
     2^task_log2 Gray-code steps (default 2^22), eliminating the v2
     starvation where 2^H / chunk_size < workers left cores idle.
   * The inner loop replaces `% M` with conditional add/subtract
     (operands are always in [0, M)), avoiding a big-integer division
     per step (~2-3x on the dominant loop).

6. DOCUMENTATION FIX: the exclusion of the two degenerate masks
   (all-A and all-B) is justified by the fact that their minimal
   POSITIVE representatives are M and M-1 respectively, both > B by
   the pre-flight check B < M — not by "both sides must be non-empty".
   (Excluding the all-A mask is in fact necessary: its raw residue is
   0, which would otherwise falsely survive as m = 0.)

Certificates produced (in the result JSON)
-------------------------------------------
all_pruned_certificate:
    No integer m with 1 <= m <= B satisfies p | m(m+1) for every
    p <= p_omega.  (Strongest form; implies the census certificate.)
census_certificate:
    No prime-complete pair (m, m+1) with m <= B exists at any level
    >= omega.  (Holds whenever no survivor passes the intrinsic
    prime-completeness check.)
Both are hereditary: they hold for every omega' >= omega at the same B.

Usage
-----
python3 CRT_Pruning_Survey.py --start_omega 33 --end_omega 40 \\
    --bound_expo 50 --workers 10

python3 CRT_Pruning_Survey.py --omega_list 28,29,30,31,32 \\
    --bound_expo 31 --workers 10

python3 CRT_Pruning_Survey.py --start_omega 30 --end_omega 35 \\
    --bound_expo 0   # bound=0: find absolute min n_CRT only

By Ken Clements, June 2026.
Version 3 changes co-developed with Claude (Anthropic).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from multiprocessing import get_context, cpu_count
from typing import Dict, List, Optional, Tuple

try:
    sys.set_int_max_str_digits(0)
except Exception:
    pass

PROGRAM_NAME = "CRT_Pruning_Survey"
PROGRAM_VERSION = 3

DEFAULT_CAP_PER_TASK = 1_000
DEFAULT_CAP_TOTAL    = 100_000
DEFAULT_PC_CAP       = 10_000   # stored prime-complete survivors (expected 0)


# ─────────────────────────────── utilities ───────────────────────────────────

def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def primes_first_n(n: int) -> List[int]:
    primes: List[int] = []
    c = 2
    while len(primes) < n:
        if all(c % p != 0 for p in primes):
            primes.append(c)
        c += 1 if c == 2 else 2
    return primes

def primes_until_primorial_exceeds(limit: int) -> List[int]:
    """Smallest prefix of the primes whose product exceeds `limit`."""
    primes: List[int] = []
    prod = 1
    c = 2
    while prod <= limit:
        if all(c % p != 0 for p in primes):
            primes.append(c)
            prod *= c
        c += 1 if c == 2 else 2
    return primes

def primorial(primes: List[int]) -> int:
    M = 1
    for p in primes:
        M *= p
    return M

def log10_primorial(primes: List[int]) -> float:
    return sum(math.log10(p) for p in primes)

def format_partition(mask: int, primes: List[int]) -> Dict:
    A = [primes[i] for i in range(len(primes)) if not ((mask >> i) & 1)]
    B = [primes[i] for i in range(len(primes)) if     ((mask >> i) & 1)]
    return {"A_divides_m": A, "B_divides_m_plus_1": B}

# ────────────────────────── CRT basis precomputation ─────────────────────────

def precompute_contrib(primes: List[int], M: int) -> List[int]:
    """
    contrib[i] = (p_i - 1) * e_i mod M, where e_i is the CRT idempotent
    (e_i ≡ 1 mod p_i, e_i ≡ 0 mod p_j for j != i).  The representative
    of partition mask is sum of contrib[i] over set bits, mod M.
    """
    contrib = []
    for p in primes:
        Mi     = M // p
        inv_Mi = pow(Mi, -1, p)
        e_i    = (Mi * inv_Mi) % M
        contrib.append((p - 1) * e_i % M)
    return contrib

# ───────────────────── intrinsic prime-completeness check ────────────────────

def pc_check(m: int, primes_check: List[int]) -> bool:
    """
    Decide whether m*(m+1) is prime-complete (at its own GPF level).

    PRECONDITION: primorial(primes_check) > m*(m+1).  Then:
      * if trial division by primes_check leaves a cofactor > 1, the GPF
        exceeds primes_check[-1]; prime-completeness would require the
        product to be divisible by a primorial larger than itself —
        impossible — so the answer is False.
      * otherwise the factorization over primes_check is complete and
        the census condition is decided exactly.
    """
    N = m * (m + 1)
    top = 0          # highest index (1-based) of a dividing prime
    count = 0        # number of distinct dividing primes
    n = N
    for i, p in enumerate(primes_check, start=1):
        if n % p == 0:
            count += 1
            top = i
            while n % p == 0:
                n //= p
            if n == 1:
                break
    if n != 1:
        return False          # GPF > p_K  ⇒  impossible (see docstring)
    return count == top       # census complete from p_1 up to the GPF

# ──────────────────────────────── worker ─────────────────────────────────────

# Per-process constants installed by the pool initializer.
_W: Dict = {}

def _init_worker(M: int, bound: int, contrib_low: List[int], low_omega: int,
                 omega: int, primes_check: List[int],
                 cap_per_task: int) -> None:
    _W["M"]            = M
    _W["bound"]        = bound
    _W["contrib_low"]  = contrib_low
    _W["low_omega"]    = low_omega
    _W["omega"]        = omega
    _W["primes_check"] = primes_check
    _W["cap"]          = cap_per_task

def _survey_worker(task: Tuple[int, int]):
    """
    Gray-code sweep over 2^L low-group assignments for one fixed
    high-group mask.  Task = (high_mask, high_contrib).

    Degenerate masks: full_mask == 0 has raw residue 0 (its minimal
    POSITIVE member is M) and full_mask == all-ones has residue M-1;
    both minimal positive members exceed B because the pre-flight check
    enforces B < M, so skipping them loses nothing — and skipping
    mask 0 is REQUIRED, since its raw residue 0 would falsely survive.

    Every survivor n_crt <= bound is immediately given the intrinsic
    prime-completeness check (independent of the storage cap).

    Returns a tuple:
      (high_mask, steps, survivor_count, cap_hit, stored_survivors,
       pc_survivors, min_n, min_mask, error)
    """
    high_mask, high_contrib = task
    try:
        t0          = time.time()
        M           = _W["M"]
        bound       = _W["bound"]
        contrib_low = _W["contrib_low"]
        low_omega   = _W["low_omega"]
        omega       = _W["omega"]
        primes_chk  = _W["primes_check"]
        cap         = _W["cap"]

        L_total        = 1 << low_omega
        high_mask_full = high_mask << low_omega
        all_ones       = (1 << omega) - 1

        stored: List[Tuple[int, int]] = []
        pc_found: List[Tuple[int, int]] = []
        survivor_count = 0
        cap_hit  = False
        min_n    = M
        min_mask = -1

        low_val   = 0
        prev_gray = 0

        for step in range(L_total):
            gray = step ^ (step >> 1)

            if step > 0:
                diff = gray ^ prev_gray
                bit  = diff.bit_length() - 1
                if (gray >> bit) & 1:
                    low_val += contrib_low[bit]
                    if low_val >= M:
                        low_val -= M
                else:
                    low_val -= contrib_low[bit]
                    if low_val < 0:
                        low_val += M

            prev_gray = gray

            full_mask = high_mask_full | gray
            if full_mask == 0 or full_mask == all_ones:
                continue

            n_crt = high_contrib + low_val
            if n_crt >= M:
                n_crt -= M

            if n_crt < min_n:
                min_n    = n_crt
                min_mask = full_mask

            if bound > 0 and n_crt <= bound:
                survivor_count += 1
                # Census check happens for EVERY survivor, cap or not.
                if pc_check(n_crt, primes_chk):
                    pc_found.append((n_crt, full_mask))
                if not cap_hit:
                    stored.append((n_crt, full_mask))
                    if len(stored) >= cap:
                        cap_hit = True

        return (high_mask, L_total, survivor_count, cap_hit,
                stored, pc_found, min_n, min_mask,
                None)

    except Exception as e:
        return (high_mask, 0, 0, False, [], [], -1, -1, repr(e))

# ─────────────────────────────── main survey ─────────────────────────────────

def survey_omega(
    omega:        int,
    bound:        int,
    workers:      int,
    logf,
    chunk_size:   int,
    task_log2:    int,
    cap_per_task: int,
    cap_total:    int,
) -> Dict:
    """
    Full partition survey for one omega.

    Task split: H high bits fixed per task, L = omega - H Gray-coded in
    the worker.  H is chosen so each task is about 2^task_log2 steps,
    with a floor that guarantees at least ~4 tasks per worker.

    Pre-flight: aborts if bound >= M.  This check is load-bearing twice
    over: (1) it guarantees the degenerate-mask exclusion is sound, and
    (2) it guarantees each surviving class has at most one member <= B,
    which is what makes the census certificate possible.
    """
    primes = primes_first_n(omega)
    M      = primorial(primes)
    log10M = log10_primorial(primes)

    if bound > 0 and bound >= M:
        msg = (
            f"SAFETY ABORT: bound {bound} >= primorial(omega={omega}) "
            f"(log10(M) = {log10M:.2f}).  Every class would survive, the "
            f"degenerate-mask exclusion would be unsound, and no "
            f"certificate of either kind is possible.  "
            f"Use --bound_expo < {int(log10M)} for omega={omega}."
        )
        logf.write(f"{utc_now_iso()} {msg}\n")
        logf.flush()
        raise ValueError(msg)

    # Census-check prime list: minimal prefix with primorial > B*(B+1).
    primes_check = (primes_until_primorial_exceeds(bound * (bound + 1))
                    if bound > 0 else [])

    # ── Task split ───────────────────────────────────────────────────────────
    H = max(omega - task_log2,
            max(1, math.ceil(math.log2(max(2, workers * 4)))))
    H = min(H, omega - 1)
    L = omega - H
    primes_high = primes[L:]

    contrib_all  = precompute_contrib(primes, M)
    contrib_high = contrib_all[L:]
    contrib_low  = contrib_all[:L]

    num_tasks = 1 << H

    logf.write(
        f"{utc_now_iso()} START omega={omega} p_max={primes[-1]} "
        f"M_digits={len(str(M))} log10(M)={log10M:.2f} H={H} L={L} "
        f"tasks={num_tasks} steps_per_task={1 << L} "
        f"bound={'10^' + str(round(math.log10(bound), 2)) if bound > 0 else 'find_min'} "
        f"census_primes={len(primes_check)} "
        f"caps={cap_per_task}/task,{cap_total} total\n"
    )
    logf.flush()

    def task_generator():
        for hm in range(num_tasks):
            hc = 0
            for i in range(H):
                if (hm >> i) & 1:
                    hc += contrib_high[i]
                    if hc >= M:
                        hc -= M
            yield (hm, hc)

    ctx = get_context("fork") if sys.platform == "darwin" else get_context()

    all_survivors: List[Tuple[int, int]] = []
    pc_survivors:  List[Tuple[int, int]] = []
    global_survivor_count = 0
    pc_survivor_count     = 0
    aggregate_cap_hit = False
    any_task_cap_hit  = False
    pc_cap_hit        = False
    global_min_n      = M
    global_min_mask   = -1
    total_steps       = 0
    errors: List[str] = []
    done              = 0

    t0 = time.time()
    with ctx.Pool(processes=workers,
                  initializer=_init_worker,
                  initargs=(M, bound, contrib_low, L, omega,
                            primes_check, cap_per_task)) as pool:
        for res in pool.imap_unordered(_survey_worker, task_generator(),
                                       chunksize=chunk_size):
            (hm, steps, s_count, c_hit, stored, pc_found,
             min_n, min_mask, err) = res
            done += 1
            total_steps += steps

            if err is not None:
                errors.append(f"task {hm}: {err}")
                continue

            global_survivor_count += s_count
            pc_survivor_count     += len(pc_found)
            if c_hit:
                any_task_cap_hit = True

            if pc_found and not pc_cap_hit:
                room = DEFAULT_PC_CAP - len(pc_survivors)
                pc_survivors.extend(pc_found[:room])
                if len(pc_survivors) >= DEFAULT_PC_CAP:
                    pc_cap_hit = True

            if not aggregate_cap_hit:
                remaining = cap_total - len(all_survivors)
                if remaining > 0:
                    all_survivors.extend(stored[:remaining])
                if len(all_survivors) >= cap_total:
                    aggregate_cap_hit = True

            if 0 < min_n < global_min_n:
                global_min_n    = min_n
                global_min_mask = min_mask

            if done % max(1, num_tasks // 20) == 0:
                logf.write(
                    f"{utc_now_iso()} progress tasks={done}/{num_tasks} "
                    f"elapsed={time.time()-t0:.1f}s "
                    f"survivors={global_survivor_count} "
                    f"pc_survivors={pc_survivor_count} "
                    f"global_min={global_min_n}\n"
                )
                logf.flush()

    elapsed_total    = time.time() - t0
    survivor_cap_hit = any_task_cap_hit or aggregate_cap_hit
    tasks_complete   = (done == num_tasks)
    run_clean        = tasks_complete and not errors

    # ── Certificates ─────────────────────────────────────────────────────────
    all_pruned         = (global_survivor_count == 0)
    all_pruned_cert    = all_pruned and run_clean
    census_cert        = (pc_survivor_count == 0) and run_clean

    all_pruned_statement = (
        f"For every level omega' >= {omega}: no integer m with "
        f"1 <= m <= 10^{round(math.log10(bound), 4) if bound > 0 else 0} "
        f"satisfies p | m(m+1) for all p <= {primes[-1]} "
        f"(hence no prime-complete pair).  Hereditary by refinement of "
        f"congruence systems."
    ) if all_pruned_cert else None

    census_statement = (
        f"For every level omega' >= {omega}: no prime-complete pair "
        f"(m, m+1) with m <= 10^{round(math.log10(bound), 4) if bound > 0 else 0} "
        f"exists.  Each surviving class has exactly one member <= B "
        f"(B < M), every such member was checked for intrinsic "
        f"prime-completeness, and none passed.  Independent of storage "
        f"caps; hereditary by refinement."
    ) if census_cert else None

    # ── Sort/dedupe stored survivors ─────────────────────────────────────────
    all_survivors.sort()
    seen: set = set()
    unique_survivors: List[Tuple[int, int]] = []
    for n, mask in all_survivors:
        if mask not in seen:
            seen.add(mask)
            unique_survivors.append((n, mask))
    unique_survivors.sort()

    valid_total = (1 << omega) - 2

    logf.write(
        f"{utc_now_iso()} DONE omega={omega} total_steps={total_steps} "
        f"valid_partitions={valid_total} "
        f"survivors={global_survivor_count} "
        f"pc_survivors={pc_survivor_count} "
        f"all_pruned={all_pruned} census_cert={census_cert} "
        f"global_min_n_crt={global_min_n} "
        f"elapsed={elapsed_total:.1f}s errors={len(errors)}\n"
    )
    logf.flush()

    survivor_details = []
    for n_crt, mask in unique_survivors[:200]:
        survivor_details.append({
            "n_crt":     n_crt,
            "n_crt_str": str(n_crt),
            "mask":      mask,
            "partition": format_partition(mask, primes),
        })

    pc_details = []
    for n_crt, mask in sorted(pc_survivors)[:200]:
        pc_details.append({
            "m":         n_crt,
            "m_str":     str(n_crt),
            "mask":      mask,
            "partition": format_partition(mask, primes),
        })

    return {
        "omega":                  omega,
        "pmax":                   primes[-1],
        "M_digits":               len(str(M)),
        "log10_M":                float(f"{log10M:.4f}"),
        "bound":                  bound,
        "bound_str":              str(bound),
        "log10_bound":            (float(f"{math.log10(bound):.4f}")
                                   if bound > 0 else None),
        "valid_partitions":       valid_total,
        "survivors":              global_survivor_count,
        "all_partitions_pruned":  all_pruned,
        "prime_complete_survivors": pc_survivor_count,
        "prime_complete_details": pc_details,
        "pc_storage_cap_hit":     pc_cap_hit,
        "census_prime_count":     len(primes_check),
        "census_prime_max":       primes_check[-1] if primes_check else None,
        "tasks_total":            num_tasks,
        "tasks_done":             done,
        "run_clean":              run_clean,
        "certificates": {
            "all_pruned_certificate":  all_pruned_cert,
            "all_pruned_statement":    all_pruned_statement,
            "census_certificate":      census_cert,
            "census_statement":        census_statement,
            "heredity_lemma": (
                "Adding the (omega+1)-th prime's congruence refines every "
                "class, so class minima are nondecreasing in omega and any "
                "certificate at omega_0 holds for all omega >= omega_0 at "
                "the same bound."
            ),
        },
        "cap_statistics": {
            "survivor_cap_hit":          survivor_cap_hit,
            "survivors_capped_per_task": cap_per_task,
            "survivors_capped_total":    cap_total,
            "survivor_count_true":       global_survivor_count,
            "survivor_count_is_exact":   True,   # counting is never capped
        },
        "global_min_n_crt":       global_min_n,
        "global_min_n_crt_str":   str(global_min_n),
        "log10_global_min":       (float(f"{math.log10(global_min_n):.4f}")
                                   if global_min_n > 0 else None),
        "global_min_partition":   (format_partition(global_min_mask, primes)
                                   if global_min_mask >= 0 else {}),
        "survivor_details":       survivor_details,
        "split_H":                H,
        "split_L":                L,
        "chunk_size":             chunk_size,
        "elapsed_sec":            elapsed_total,
        "total_steps":            total_steps,
        "errors":                 errors[:20],
    }

# ──────────────────────────────────── main ───────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=f"{PROGRAM_NAME} v{PROGRAM_VERSION}: "
                    "CRT partition survey with census certificates."
    )
    ap.add_argument("--start_omega",  type=int, default=30)
    ap.add_argument("--end_omega",    type=int, default=36)
    ap.add_argument("--omega_list",   type=str, default="",
                    help="Comma-separated omega values (overrides start/end).")
    ap.add_argument("--bound_expo",   type=int, default=40,
                    help="B = 10^E.  0 = find minimum n_CRT only.  "
                         "Must satisfy 10^E < primorial(omega).")
    ap.add_argument("--workers",      type=int, default=0)
    ap.add_argument("--chunk_size",   type=int, default=1,
                    help="Tasks per pool chunk (default 1; tasks are large).")
    ap.add_argument("--task_log2",    type=int, default=22,
                    help="Target log2(Gray-code steps) per task (default 22).")
    ap.add_argument("--outdir",
                    default=f"{PROGRAM_NAME}_v{PROGRAM_VERSION}_results")
    ap.add_argument("--skip_done",    action="store_true")
    ap.add_argument("--survivor_cap_per_task", type=int,
                    default=DEFAULT_CAP_PER_TASK,
                    help="Max survivors STORED per task (counting and the "
                         "census check are never capped).")
    ap.add_argument("--survivor_cap_total", type=int,
                    default=DEFAULT_CAP_TOTAL,
                    help="Max survivors stored in the aggregate list.")
    args = ap.parse_args()

    workers = args.workers if args.workers > 0 else (cpu_count() or 4)
    bound   = 10 ** args.bound_expo if args.bound_expo > 0 else 0

    if args.omega_list:
        omegas = [int(x.strip()) for x in args.omega_list.split(",")
                  if x.strip()]
    else:
        omegas = list(range(args.start_omega, args.end_omega + 1))

    print(f"[+] {PROGRAM_NAME} version {PROGRAM_VERSION}")
    print(f"[+] omega range: {omegas}")
    if bound > 0:
        print(f"[+] bound B = 10^{args.bound_expo}")
        print(f"[+] survivor storage caps: {args.survivor_cap_per_task}/task, "
              f"{args.survivor_cap_total} aggregate "
              f"(counting and census checks are uncapped)")
    else:
        print("[+] bound: find minimum n_CRT only (no pruning)")
    print(f"[+] workers={workers}  chunk_size={args.chunk_size}  "
          f"task_log2={args.task_log2}")
    print(f"[+] outdir: {args.outdir}")
    print()

    if bound > 0:
        print(f"  {'omega':>5}  {'log10(M)':>10}  {'log10(B)':>9}  "
              f"{'B < M?':>8}  {'est. survivor count':>22}")
        print("  " + "-" * 62)
        for om in omegas:
            ps = primes_first_n(om)
            lm = log10_primorial(ps)
            lb = args.bound_expo
            safe = "YES" if lb < lm else "*** NO — WILL ABORT ***"
            if lb < lm:
                est = (2 ** om) * 10 ** (lb - lm)
                est_str = f"~{est:.2e}"
            else:
                est_str = "ALL"
            print(f"  {om:>5}  {lm:>10.2f}  {lb:>9}  {safe:>8}  {est_str:>22}")
        print()

    print(f"[+] start: {utc_now_iso()}\n")
    ensure_dir(args.outdir)

    for omega in omegas:
        result_path = os.path.join(args.outdir,
                                   f"crt_survey_omega_{omega:02d}.json")
        log_path    = os.path.join(args.outdir,
                                   f"crt_survey_omega_{omega:02d}.log")

        if args.skip_done and os.path.isfile(result_path):
            print(f"[=] omega={omega}: result exists, skipping.")
            continue

        with open(log_path, "a", encoding="utf-8") as logf:
            try:
                result = survey_omega(
                    omega        = omega,
                    bound        = bound,
                    workers      = workers,
                    logf         = logf,
                    chunk_size   = args.chunk_size,
                    task_log2    = args.task_log2,
                    cap_per_task = args.survivor_cap_per_task,
                    cap_total    = args.survivor_cap_total,
                )
            except ValueError as exc:
                print(f"[!] omega={omega}: {exc}")
                continue

        result["sha256_log"] = sha256_file(log_path)
        result["program"]    = PROGRAM_NAME
        result["version"]    = PROGRAM_VERSION
        result["run_utc"]    = utc_now_iso()
        result["platform"]   = platform.platform()

        with open(result_path, "w", encoding="utf-8") as f:
            json.dump(result, f, indent=2, sort_keys=True)

        certs = result["certificates"]
        if certs["all_pruned_certificate"]:
            cert_str = "  *** ALL-PRUNED CERTIFICATE (hereditary) ***"
        elif certs["census_certificate"]:
            cert_str = "  *** CENSUS CERTIFICATE (hereditary) ***"
        else:
            cert_str = ""
        pc_str = (f"  !!! {result['prime_complete_survivors']} "
                  f"PRIME-COMPLETE SURVIVOR(S) !!!"
                  if result["prime_complete_survivors"] > 0 else "")

        print(
            f"[>] omega={omega:2d}  pmax={result['pmax']:4d}  "
            f"M_digits={result['M_digits']:3d}  "
            f"survivors={result['survivors']:>12,}  "
            f"pc={result['prime_complete_survivors']}  "
            f"min_n_crt=10^{result['log10_global_min']:.2f}  "
            f"time={result['elapsed_sec']/60:.2f}min"
            f"{cert_str}{pc_str}"
        )

    print(f"\n[+] Finished. {utc_now_iso()}")


if __name__ == "__main__":
    if sys.platform == "darwin":
        try:
            import multiprocessing as mp
            mp.set_start_method("fork")
        except RuntimeError:
            pass
    main()
