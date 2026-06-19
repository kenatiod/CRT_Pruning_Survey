# CRT_Pruning_Survey

**A pure Chinese-Remainder-Theorem survey that rules out prime-complete products
of consecutive integers at high census levels — no Pell solver required.**

This program supplies one half of the computational evidence behind the
conjecture that

> **633,555 × 633,556 is the last prime-complete product of two consecutive
> integers** (OEIS [A141399](https://oeis.org/A141399)).

For each census level ω it surveys, it produces a *hereditary certificate*: a
machine-checkable statement that no prime-complete pair $m(m+1)$ exists below a
stated bound, at that level **and every level above it**. Run across a contiguous
range of ω, the certificates form a wall that the smooth-pair "ceiling" cannot
climb back over — the structural fact that makes high-ω resumption of solutions
(re-entry) impossible.

---

## Background: prime-complete pairs and the census

Fix ω and let $P_ω = {2, 3, 5, …, p_ω}$ be the first ω primes, with primorial
$M_ω = 2·3·5···p_ω$. A pair of consecutive integers $m, m+1$ is **prime-complete
at level ω** when every census prime divides the product and no larger prime does:

$$
rad( m(m+1) ) = P_ω.
$$

Such pairs are rare and they stop: the last one is $m = 633555$ (at ω = 8). The
question is whether any exist at higher ω. Because $m$ and $m+1$ are coprime, a
prime-complete pair corresponds to a **partition** of the census $P_ω = A ⊔ B$,
with the primes in $A$ dividing $m$ and those in $B$ dividing $m+1$. Each
partition fixes a residue class of $m$ modulo $M_ω$ by the Chinese Remainder
Theorem. So the entire question at level ω is: does any of the $2^ω − 2$
non-degenerate CRT classes contain a prime-complete $m$?

---

## What the program does

For each ω the survey enumerates all valid partitions $A ⊔ B = P_ω$ (bit *i* of a
mask = 1 means $p_i | m+1$, else $p_i | m$) and computes, by a Gray-code sweep
over a precomputed CRT idempotent basis, the **minimal positive representative**
$n_{A,B}$ of each class modulo $M_ω$. Each representative $≤ B$ is then tested,
in the worker, for the intrinsic prime-completeness of $m(m+1)$.

The two degenerate masks (all-A, all-B) are excluded soundly: their minimal
positive representatives are $M_ω$ and $M_ω − 1$, both $> B$ whenever $B < M_ω$.

### The floor that makes it finite

A prime-complete pair of order ω factors as $m(m+1) = M_ω · jk$ with $j, k ≥ 1$,
so $m(m+1) ≥ M_ω$ and therefore

$$
m ≥ √M_ω − 1. \text{           (the "floor")}
$$

This is an exact, unconditional lower bound on where a prime-complete pair could
*first* appear at level ω. It is the reason a bounded survey can certify an
unbounded range, and it is the quantity the safety margin is measured against
(see below).

### Two certificates

Each run records, in its result JSON, whichever of these it can prove:

- **all-pruned certificate** — no integer $m$ with $1 ≤ m ≤ B$ satisfies
  $p | m(m+1)$ for every $p ≤ p_ω$. The strongest form: every CRT class was
  surveyed and none produced a candidate at or below $B$.
- **census certificate** — no prime-complete pair $(m, m+1)$ with $m ≤ B$
  exists at any level $≥ ω$. Holds whenever no surviving representative passes
  the intrinsic prime-completeness check, even if some raw survivors exist.

**Both are hereditary.** Adding the (ω+1)-th prime's congruence only refines the
existing classes, so class minima are non-decreasing in ω: a certificate proved
at $ω₀$ holds for every $ω ≥ ω₀$ at the same bound B. One run buys an infinite
tail of levels.

---

## Soundness conditions

The certificates are theorems, not search summaries, but only when their
preconditions hold. The program enforces them and refuses to emit a misleading
certificate otherwise:

1. **$B < M_ω$ (pre-flight abort).** Guarantees each surviving class holds at most
   one member $≤ B$ (its minimum), which is what makes the census certificate
   valid, and makes the degenerate-mask exclusion sound. A run with $B ≥ M_ω$
   aborts.
2. **Clean completion.** A certificate requires zero worker errors **and** all
   tasks completed; a crashed worker cannot yield a false "all pruned."
3. **Uncapped checking.** Survivor *storage* may be capped for memory, but
   survivor *counting* and the prime-completeness *check* are never capped —
   every survivor is tested at the moment it is found.
4. **Non-vacuous bound — the floor margin (v4).** $B < M_ω$ keeps a certificate
   *sound*, but not necessarily *meaningful*: as ω grows the floor $√M_ω$ rises,
   and a fixed bound eventually drops below it, certifying only a range in which
   the floor identity already guarantees no pair can lie. v4 measures and records
   the **margin above the floor**, $log_{10}(B) − ½·log_{10}(M_ω)$, and can hold it
   uniform across a whole run (see --floor_margin). A certificate is
   proof-relevant only when this margin is $≥ 0$, and carries the intended safety
   factor when it is comfortably positive.

---

## The intrinsic prime-completeness check

For a survivor $m ≤ B$, let $N = m(m+1)$. Trial-divide $N$ by the first $K$
primes, where $K$ is the smallest index with $primorial(K) > B·(B+1) ≥ N$. Then:

- if a cofactor $> 1$ remains, the greatest prime factor of $N$ exceeds $p_K$, so
  prime-completeness would force $N$ to be divisible by a primorial larger than
  itself — impossible — and the answer is **not prime-complete**;
- if the cofactor is $1$, the factorization over the first $K$ primes is complete
  and the census condition is decided **exactly** (the dividing prime indices must
  be ${1, …, g}$ where $p_g$ is the greatest prime factor).

No floating point or probabilistic test enters the decision.

---

## Requirements

- **Python ≥ 3.10** (uses $pow(a, -1, p)$ modular inverse).
- No third-party packages for the core program; only the standard library.

The survey is CPU-bound and embarrassingly parallel across CRT partitions; it
uses $multiprocessing$ and scales with available cores.

---

## Usage

The bound $B$ is the integer below which the survey is exhaustive. There are two
ways to set it.

### Uniform floor margin (recommended for a multi-ω run)

Hold a fixed safety margin (in decades) above the floor $√M_ω$ at every level:

python3 CRT_Pruning_Survey.py --start_omega 34 --end_omega 70 \
    --floor_margin 6 --min_floor_margin 6 --workers 10


This sets $log_{10}(B) = ½·log_{10}(M_ω) + 6$ per level, so every certificate in the
run is six decades clear of the floor. --min_floor_margin stops the run cleanly
the moment that margin cannot be held.

Add --dry_run to print the pre-flight table — floor, bound, $B < M$, and the
margin above the floor for each ω — and exit without surveying. Always worth
running first before a long job.

### Fixed bound

A single exponent $B = 10^E$ for all levels (suitable for one ω, or when you want
an identical bound across a short range):

python3 CRT_Pruning_Survey.py --omega_list 28,29,30,31,32 \
    --bound_expo 31 --workers 10


Note that with a fixed $E$ the floor margin shrinks as ω grows and can go
negative; the pre-flight table shows this, and --min_floor_margin can guard
against it.

### Minimum-only mode

Find the absolute minimum CRT representative at each level (no pruning):

bash
python3 CRT_Pruning_Survey.py --start_omega 30 --end_omega 35 --bound_expo 0


### Key options

| Option | Meaning |
|---|---|
| --start_omega, --end_omega | Inclusive census-level range. |
| --omega_list a,b,c | Explicit levels (overrides start/end). |
| --floor_margin D | Per-ω bound held $D$ decades above the floor $√M_ω$. Overrides --bound_expo. |
| --min_floor_margin D | Abort/stop any level whose floor margin would fall below $D$. |
| --bound_expo E | Fixed $B = 10^E$ for all levels. $0$ = minimum-only. Requires $10^E < M_ω$. |
| --dry_run | Print the pre-flight table and exit. |
| --workers N | Worker processes (default: all cores). |
| --task_log2 K | Target $log₂(steps)$ per task (load balancing; default 22). |
| --survivor_cap_per_task, --survivor_cap_total | Caps on *stored* survivors only; counting and checking are never capped. |
| --skip_done | Skip levels whose result file already exists. |
| --outdir DIR | Output directory (default CRT_Pruning_Survey_results). |

---

## Output

Each level writes crt_survey_omega_NN.json and crt_survey_omega_NN.log to the
output directory. The JSON records:

- the verdict flags $all_pruned_certificate$ and $census_certificate$, each with
  its full statement string;
- floor_log10, floor_margin_decades, floor_margin_nonneg — the
  proof-relevant margin, so each certificate self-documents its own validity;
- survivor and prime-complete-survivor counts (exact — never capped), the global
  minimum representative, partition details for any survivors;
- the heredity-lemma statement, run-clean status (errors and task completion),
  timing, and a SHA-256 of the log for integrity.

Sample certified results are in
[CRT_Pruning_Survey_results/](CRT_Pruning_Survey_results).

A clean, proof-relevant certificate at level ω satisfies, in its JSON:
$all_pruned_certificate$ (or at least census_certificate) = true,
run_clean = true, and floor_margin_nonneg = true.

---

## How this fits the larger proof

This survey is the high-ω instrument of a two-part argument. Cast as a halting
question, the conjecture says a search procedure that looks for a prime-complete
pair at any level beyond ω = 8 **never halts**. The proof of non-halting has two
regimes, and this program covers the first directly:

- **Certified base case (this repository).** For every level up to the highest ω
  surveyed, the all-pruned / census certificate proves — hereditarily and with a
  uniform safety margin above the floor — that the search produces no candidate.
  This is recorded computation, not extrapolation: each level names its margin and
  its clean-run status in an auditable JSON.
- **Asymptotic tail (companion write-up).** Beyond the surveyed range, a single
  classical smooth-pair ceiling bound shows the floor $√M_ω$ has permanently
  overtaken the largest possible smooth pair, so the search guard becomes
  unsatisfiable. The longer the certified base case runs, the deeper that handoff
  sits inside already-exhausted territory.

The companion enumerator LCm_Solver certifies the complementary low-ω range by
the Lehmer–Clements method; together the two cover the finite base case that the
asymptotic argument builds on.

---

## License

MIT — see [LICENSE](LICENSE).

## Author

Ken Clements, 2026. Version 4 algorithmic and safety-margin changes co-developed
with Claude (Anthropic).
