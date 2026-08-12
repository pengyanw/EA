# CLAUDE.md — CDC 2026 Submission 642: Revision Tracker

Paper: *An Evolutionary Algorithm for Actuator-Sensor-Communication Co-Design in Distributed Control*
Authors: Pengyang Wu, Jing Shuang (Lisa) Li
Status: **Accepted as Invited Session Paper** (decision 2026-07-15)
Final submission deadline: **2026-09-18** (PaperPlaza)
Page budget: 6 pages included; up to 8 pages with over-length charge ($140/page beyond 6).

---

## PRIMARY REVISION TARGET

Per the Associate Editor, the two weaknesses that must be addressed above all else:

1. **The thoroughness of the convergence/stability analysis.**
2. **The clarification of the novel contribution brought by the newly developed Algorithm 1.**

The governing reviewer comment (R3-5):

> "Convergence analysis is informative but does not provide a rigorous description of the EA.
> It rather studies a perfected pruning process which the EA can follow or not follow.
> There is a large discrepancy between algorithm and theory."

**Everything in Sections A and B below serves these two points. Treat the rest as secondary.**

> 🛑 **Before any of it: see the BLOCKER section.** There is an unresolved inconsistency
> between the paper's description of the link-count coordinate `ℓ` and the observed simulation
> behavior. It may mean the implemented algorithm differs from the described one, which would
> change what needs writing.
Detailed technical plan: see "Primary Plan" at the end of this file.

> ⚠️ Reviewer 1's comments are referenced by the AE but were NOT shown on the acceptance page.
> Pull them from PaperPlaza before finalizing and fold them in.

---

## Priority A — Theory–Algorithm Consistency (AE point 1 + Reviewer 3)

- [ ] **A1 — Step-size mismatch invalidates Prop. 5.** Prop. 5 assumes the pure Polyak step (37); Algorithm 2 line 11 implements `η ← min(Polyak, 0.5)`. *(R3-6)*
- [ ] **A2 — Repair step vs. convergence.** Sec. VI claims repair "does not affect ... convergence properties of the original algorithm" with no argument. Repair changes the objective the EA minimizes (infinite → finite on unstable genes). *(R3-7, R3-8)*
- [ ] **A3 — Idealized pruning vs. actual EA.** Sec. IV analyzes a monotone 1-D pruning process under `a = 1, s = 1`; Algorithm 1 searches the full `θ = [ℓ, a, s]` space with mask mutation active. *(R3-5 — the headline item)*
- [ ] **A4 — Algorithm 1 internal inconsistencies.** *(R3-3, R3-4)*
- [ ] **A5 — Subgradient at non-differentiable points.** Prop. 4 handles `μ_j(K) = 0` correctly in the text, but Algorithm 2 line 10 uses `sign(·)` unconditionally. State the convention. *(R3-9)*

## Priority B — Contribution & Method Justification (AE point 2 + Reviewer 3)

- [ ] **B1 — State the novel contribution.** *(R3-1)*
- [ ] **B2 — Justify EA over baselines.** *(R3-2)*
- [ ] **B3 — Justify the specific operators.** *(R3-3)*

## Priority C — Problem-Formulation Clarity (Reviewers 5 & 2)

- [ ] **C1 — Unified objective / trade-off** in abstract + intro. *(R5-1)*
- [ ] **C2 — Assumption 1 rigor.** *(R5-2)*
- [ ] **C3 — Coupling of control/comm/actuator dynamics.** *(R2-2)*

## Priority D — Presentation & Framing (Reviewers 2 & 5)

- [ ] **D1 — Introduction + literature review** toward full SOTA, within page budget. *(R2-1)*
- [ ] **D2 — Conclusion** needs a real takeaway. *(R2-3)*
- [ ] **D3 — Grammar pass.** Abstract: "The proposed methods is validated" → "method is". *(R5-3)*

---

# BLOCKER — Resolve before writing any new Section IV text

**Status: unresolved arithmetic inconsistency between the paper's description of `ℓ` and the
observed simulation behavior. Settle this first; the form of the convergence theorems depends
on it.**

## The observation

Iterating `ℓ` alone (masks fixed at `a = s = 1`) already produces substantially sparse
controllers. This rules out the hypothesis that mask flips are the sole driver of sparsity.

## The structural bound that contradicts it

`ℓ` has a hard per-run reachability limit under the algorithm as written in the paper:

- initialization sets **every** individual to `ℓ = ‖K_d‖_0` — zero initial diversity in `ℓ`;
- crossover only *copies* existing `ℓ` values between individuals, never creates new ones;
- mutation shifts by `|δ| ≤ d`.

Therefore the population minimum of `ℓ` decreases by at most `d` per generation:

    ℓ_min ≥ ℓ_0 − d · G_max

**Reachability Lemma (worth stating in the paper regardless of how this resolves).**

For the 7×7 grid: `ℓ_0 = N_u·N_x = 49·98 = 4802`, `d = 5`, `G_max = 150`
⟹ `ℓ` cannot go below **4052**, i.e. at least 84% of entries are retained.

## Why that cannot produce the observed cost

With Section VII weights (`w_a = 0.4`, `w_s = 0.2`, `w_c = 0.05`) on the 7×7 grid:

    J_EA(K_d) ≈ 1 + 0.4·49 + 0.2·98 + 0.05·2352
              = 1 + 19.6 + 19.6 + 117.6  ≈ 157.8

With `a = s = 1` the masks contribute a fixed 39.2, so every remaining reduction must come from
`w_c·N_c`. The `ℓ`-only floor is `(1 + 39.2)/157.8 ≈ 0.255` normalized — which is roughly what
Fig. 1 shows for 7×7, so the *cost* is consistent. But reaching it requires `N_c → 0`, i.e.
retaining only the diagonal blocks `[K]_ii`. With 2 states and 1 input per bus, that is
`49 × 2 = 98` entries, so `ℓ` must travel `4802 → ≈98`, a displacement of **4704 against a
budget of 750**. Off by more than 6×.

## Candidate resolutions (ranked by suspicion) — for CC to check in the code

1. **`ℓ` is not an element count in the implementation.** If it is a retained *fraction*
   (`[0,1]` or a percentage), then `d = 5` means ±5% per generation and 150 generations sweep the
   whole range — everything reconciles. The paper says `ℓ ∈ [1, N_uN_x]` with `sat(·)` clipping
   to that interval.
2. **Mutation on `ℓ` is multiplicative or log-scaled** (`ℓ ← ℓ·(1±δ)`), not additive as written.
3. **Initialization is not `ℓ = ‖K_d‖_0`** but random over the interval. Then the initial
   population already spans a wide `ℓ` range, crossover can jump straight into the sparse region,
   and the `d·G_max` bound does not apply.
4. **`d` in the run configuration is not 5.**

## Diagnostic — three numbers settle it

From any 7×7 run, report:

- `ℓ_0`, and the **distribution** of `ℓ` across the initial population (single point or spread?)
- terminal `ℓ*`
- terminal `N_c(K_s(θ*))`

Interpretation:
- `ℓ*` pinned near `ℓ_0 − d·G_max` ⟹ the run is **rate-limited by `d`, not converged**. (Unlikely:
  Fig. 1's 7×7 curve plateaus around generation 100, which is not the signature of rate limiting.)
- `ℓ*` far below 4052 ⟹ `ℓ` is parameterized differently from the paper's description. Fix the
  paper text to match the code.

Also worth logging per generation, for the same runs: `ℓ*_t`, `N_c`, `N_a`, `N_s` trajectories.
These make the pruning dynamics visible and would strengthen Section VII on their own.

## Why this outranks everything else on the list

If the discrepancy is that **the implemented algorithm differs from the algorithm described in
the paper**, then the analysis is not merely loose — its subject does not exist. That is a
stronger version of R3-5 and must be closed before any other revision work.

It also determines the theory: the shape of the one-step neighborhood `N_d(θ)` (additive `±d`
vs. multiplicative/fractional) fixes the form of `p_el` and hence Theorem 3.

**Unaffected either way:** Theorem 1 (a.s. finite-time convergence) and Theorem 2 + Corollary
(stalled set, mask-global optimality) in `section-iv-new.tex`. They use only finiteness of `Θ`
and positivity of the per-neighbor mutation probability, not the parameterization of `ℓ`.
**Affected:** Theorem 3 and the whole `h*` / `Φ(h)` apparatus.

---

# PRIMARY PLAN — Closing the Algorithm/Theory Gap

## Part 0 — Concrete defects located in the current text

These are specific, checkable, and are almost certainly what R3 saw.

| # | Location | Defect |
|---|---|---|
| P1 | Sec. III, Mutation | `δ ~ Unif{−d, d}` reads as a two-point set {−d, +d}. Cor. 1 and Prop. 2 both use `1/(2d+1)`, i.e. the integer *range*. If literally two-point, `δ = −1` is impossible and **the improvement path in Prop. 2 does not exist**. Fix to `Unif{−d, …, d}` (integers). |
| P2 | Prop. 2 | The improvement path claims the offspring has "identical masks" to the elite. But `θ^c = [θ¹_{1:k}, θ²_{k+1:N_θ}]` with `k ~ U(1, N_θ−1)` means the child **always** inherits at least the last coordinate from parent 2. The offspring can never equal the elite. The stated path is not realizable as written. |
| P3 | Prop. 2 vs. Thm. 1 | Prop. 2 gives `p_imp ≥ (1−p_m)^{N_u+N_x} / (N_p(2d+1)) · P(K_t ∈ S)`. Thm. 1 then **defines** `p_imp := 1/(N_p(2d+1)) · 1{h>h*}`, silently dropping the mask-preservation and stability factors. Thm. 1's bound is therefore *larger* than what Prop. 2 supports — Thm. 1 is not established by its own lemma. |
| P4 | Lemma 3, converse | "when `h_t ≤ h*`, pruning additional link(s) does not decrease `J_EA`" is **false as stated**. `Φ(h)` is an upper bound on the cost increase; `Φ(h) ≥ w_c` means the *certificate* fails, not that pruning hurts. Thm. 1's stagnation conclusion rests on this. |
| P5 | Thm. 1(ii) | `ℓ*_∞ ≤ h^{−1}(h*)` treats `h(ℓ)` as invertible. It is a non-decreasing step function. Use a generalized inverse in the style already used in Cor. 1 (`ℓ_stab := min{ℓ : h(ℓ) ≥ h_stab}`). |
| P6 | Thm. 1 | The `N_p − n_e` offspring are treated as independent; they share the parent pool `P_{t−1}` and the same softmax weights. Independence needs justification or a union/Chernoff-free restatement. |
| P7 | Sec. VII vs. Eq. (6) | Simulations use `τ = 0`, at which the softmax `exp(−J/τ)` is undefined. Needs an explicit convention (deterministic argmin) plus a note that the `≥ 1/N_p` selection bound survives it. |
| P8 | Alg. 1 line 9–10 | Line 9 passes `p_c` but the Crossover operator description never uses it (crossover appears unconditional). Line 10 reads `Mutation(θ, p_m, d)` — argument should be `θ^c`. |
| P9 | Def. 2(a) | Uses `F(θ)`, never defined. Should be `J_EA(θ)`. |
| P10 | Sec. V, Eq. (22) | "Rearranging (21) gives the unit-decrease identity" — (22) does **not** follow from (21); it is the defining Lyapunov equation `V − (A+BK_d)ᵀV(A+BK_d) = I`. Restate. |
| P12 | Def. 2(b) | Asserts `|X_t| ≤ d` for `X_t = ℓ*_{t−1} − ℓ*_t`. False: an incoming crossover elite inherits `ℓ` from an arbitrary parent, so `ℓ*_t` can jump by more than `d` and is not monotone. Elitism makes the cost monotone, not the link count. Use `J_EA` as the potential function. |
| P11 | Prop. 1 / Thm. 1 proofs | Cross-references are wrong ("From Lemma 2, we have …" where the step is just the definition of `J_EA`; "combine Lemmas 2 and 1" should route through Lemma 4). |

## Part 1 — The structural fix: make the algorithm match the theory, not the reverse

**Recommended headline change: add a dedicated elite-mutation offspring to Algorithm 1.**

Currently every non-elite offspring is produced by `Selection → Crossover → Mutation`, and crossover destroys the improvement path (P2). Instead, reserve one offspring slot per generation:

```
θ^elite-child ← Mutation(θ*_{t−1}, p_m, d)        # no crossover
```

This is a two-line change (an embedded (1+λ)-EA inside the GA), costs nothing, and makes the
improvement path **exact rather than approximate**:

    p_elite = (1/(2d+1)) · (1 − p_m)^{N_u+N_x}

with no selection factor, no crossover factor, no independence assumption. The remaining
`N_p − n_e − 1` offspring can only add improvement probability, so `p_elite` is a valid lower
bound on `P_imp(t)`. This single change resolves P2, P3, P6, and most of P7 at once.

**Framing for the paper:** this is not a patch — it is the point. Present it as "Algorithm 1 is
designed so that a per-generation improvement certificate is available," which is precisely the
answer to AE point 2 (what is novel about Algorithm 1 vs. a textbook GA).

## Part 2 — Restate the convergence result as a drift + hitting-time theorem on the true chain

**STATUS: drafted. See `section-iv-new.tex` for compilable LaTeX (revised Algorithm 1 + Thms 1-4).**

Replace "the EA follows a pruning path" with "the EA is a finite-state elitist chain with a
certified drift." Four results, in dependency order:

**Prop. (elite-child hitting probability).** With the line-9 elite-mutation slot, for any
`θ' ∈ N_d(θ*_{t−1})`, `P(θ^el = θ' | F_{t−1}) ≥ q_min := (min{p_m, 1−p_m})^{N_u+N_x}/(2d+1) > 0`.
No selection or crossover factor enters. This is the technical payoff of the Part 1 change.

**Thm 1 (a.s. finite-time stabilization).** `Θ` is finite (`|Θ| = N_uN_x·2^{N_u+N_x}`), so
`J_EA(Θ)` is a finite value set; elitism makes `J_EA(θ*_t)` non-increasing; a non-increasing
sequence in a finite ordered set is eventually constant. Gives **almost-sure, finite-time**
convergence — strictly stronger than the submitted expectation-based statement, and it needs
neither Lemma 3/4 nor `h*`.

**Thm 2 (limit set) + Cor. (mask-global optimality).** Because each entry of `m_a, m_s` is an
independent Bernoulli(`p_m`) flip, the one-step neighborhood `N_d(θ)` constrains only `ℓ`
(to `±d`) and leaves the masks **unconstrained**. Hence `θ*_t → Θ_stall` a.s., and taking
`δ = 0` yields: the returned `(a, s)` is **globally optimal over the entire mask space at the
returned link count**. This is a genuinely new and quotable statement; it also explains the very
sparse actuator selections in Fig. 1.

**Thm 3 (expected certified-phase length).** With `Δ_min := w_c − Φ(h*+1) > 0` and
`p_el := (1−p_m)^{N_u+N_x}/(2d+1)`, `E[N_cert] ≤ J_EA(θ*_0)/(Δ_min · p_el)` via geometric
trials + Wald.

⚠️ **Use `J_EA` as the potential, not `ℓ*_t`.** Submitted Def. 2(b) asserts `|X_t| ≤ d` for
`X_t = ℓ*_{t−1} − ℓ*_t`. False for the submitted algorithm: an incoming crossover elite inherits
`ℓ` from an arbitrary parent, so `ℓ*_t` can jump by more than `d` and is not monotone. Elitism
makes the *cost* monotone, not the link count. **Log as defect P12.**

**Remark (mutation-rate scaling) — the practically valuable by-product.** Thm 3 degrades as
`(1−p_m)^{N_u+N_x}`. At the Section VII settings (`p_m = 0.05`, `d = 5`):
`p_el ≈ 1.9e−3` for the 5×5 grid (`N_u+N_x = 75`) but `≈ 4.8e−5` for the 7×7 grid (`147`),
i.e. ~2e4 generations per certified decrement against `G_max = 150` — **the bound is vacuous at
that scale**, even though the EA works well. Report this honestly: the bound counts a single
improvement path while the true neighborhood is far larger. The useful inference is that *mask
preservation, not link search, is the binding constraint*, which points to the standard
length-scaled rate `p_m = c/(N_u+N_x)`, giving `(1−p_m)^{N_u+N_x} → e^{−c}`, dimension-free.
At `c = 1`: `p_el ≈ 3.3e−2`, ~30 generations per decrement regardless of size.
**Testable prediction: on the 7×7 grid, `p_m = 1/147` should beat `p_m = 0.05`. Run this.**
This is strong material for AE point (ii): the analysis *designed* a parameter rule rather than
merely describing the algorithm after the fact.

**Optional extension (only if pages allow):** define `h*_a`, `h*_s` analogously to `h*` by
balancing `w_a`, `w_s` against the `L_J/2‖·‖²_F` cost of removing an actuator row / sensor
column, giving a certified region in all three coordinates.

## Part 3 — Repair (Algorithm 2): two fixes

**A1 (step size).** The clip is *not* fatal — prove the damped-Polyak version and the clip
becomes a special case. For `η_t = λ_t (R̄ − ρ*)/‖g̃‖²_F`:

    ‖K^{(t+1)} − K*‖²_F ≤ ‖K^{(t)} − K*‖²_F − λ_t(2 − λ_t)(R̄(K^{(t)}) − ρ*)² / ‖g̃^{(t)}‖²_F

which decreases for any `λ_t ∈ (0, 2)`. The implemented clip is exactly
`λ_t = min(1, 0.5‖g̃‖²_F/(R̄ − ρ*)) ∈ (0, 1]`. Since `‖g̃‖_F ≤ ‖B‖_F √N_x` and the iterates stay
in a bounded set (the distance to `K*` is non-increasing), `λ_t ≥ λ_min > 0` explicitly.
Substituting recovers (38) with an extra `λ_min(2 − λ_min)` factor. **Prop. 5 then describes the
code that was actually run.** State `λ_min` in closed form.

**A2 (repair vs. EA convergence).** Delete the unsupported sentence and replace with:
- Define `J^r_EA` explicitly as the objective minimized when repair is on. Note `J^r_EA ≤ J_EA`
  pointwise (finite vs. `+∞` on unstable genes), and `J^r_EA = J_EA` on stable genes — so the
  elite sequence is still monotone and Part 2's drift argument transfers verbatim.
- Flag honestly what does *not* transfer: Lemma 4 bounds `J_LQR(K_t) − J_LQR(K_{t−1})` using
  `∇J(K_d) = 0` and `‖K_{t−1} − K_d‖_F ≤ √(N_uN_x) Υρ^{h}`. A repaired gain `K^r` has moved off
  the pruned-`K_d` family, so this needs the extra term `‖K^r − K_s‖_F ≤ Σ_t η_t ‖g̃^{(t)}‖_F`,
  which the repair analysis already bounds. Add it, or scope Lemma 4 to unrepaired elites.
- Also state the fallback when `G ∩ K_S = ∅` (Prop. 5's hypothesis): Algorithm 2 exits after `T`
  iterations with a possibly unstable gain, which is then correctly assigned infinite cost. And
  add a remark that the Gershgorin condition `R̄(A+BK) < 1` is sufficient but conservative.

## Part 4 — What is actually novel about Algorithm 1 (answers AE point 2 + B1/B2/B3)

Three claims, all supportable from the current paper:

1. **The encoding is the contribution.** `θ = [ℓ, a, s]` replaces a search over `2^{N_uN_x}`
   sparsity patterns with a search over a *totally ordered chain* (`ℓ` indexing the nested family
   `Π_ℓ(K_d)`) crossed with two hypercubes. The chain is ordered by `|K_{d,ij}|`, which by Lemma 1
   is aligned with graph distance — **the theoretical prior of Lemma 1 is compiled into the search
   space.** This is why the drift certificate exists at all; a GA over raw support masks would
   admit no such analysis. Say this explicitly in the intro and next to Algorithm 1.
2. **The repair mechanism.** A sparsity-preserving convex-feasibility projection that converts
   infinite-cost genes into usable ones without touching the support. This is the clearest
   differentiator from any off-the-shelf GA, and Fig. 2 quantifies it (25% → 35%).
3. **The certified drift.** After Parts 1–2, Algorithm 1 is a GA with an explicit per-generation
   improvement lower bound and an expected hitting-time bound tied to plant spatial decay.
   Generic EAs have no such guarantee.

**Why an EA and not something simpler (B2):**
- *Greedy pruning*: `J_EA` is not submodular and not even monotone in the support — removing one
  link can destabilize the closed loop, an infinite jump. No approximation guarantee. Empirically
  this is the "diagonal LQR" / naive-truncation regime the EA already beats by 28–52%.
- *Mixed-integer optimization*: `J_LQR(K_s(θ))` requires a Lyapunov solve per candidate and is
  non-convex, non-smooth, and `+∞` on the unstable set — there is no MIP-representable closed
  form. Big-M + LMI relaxations do not scale to `N_x = 98`.
- *Random search*: **add it as a baseline in Fig. 1.** It reuses the same cost oracle, costs
  almost nothing to run, and converts B2 from an assertion into evidence. Strongly recommended.

**Operator justification (B3):** single-point crossover with `ℓ` at position 1 means the child
always inherits the parent-1 link count — i.e. crossover recombines *architecture* (masks) while
preserving *density*, which matches the product structure of the encoding. Bit-flip mutation is
the natural neighborhood on the hypercube coordinates; bounded additive mutation `±d` is the
natural neighborhood on the chain coordinate, and `d` directly sets the granularity of the
pruning steps the drift bound relies on. Say this — currently it reads as arbitrary.

## Part 5 — Suggested order of work

0. **Resolve the BLOCKER above.** Nothing below is safe to write until the `ℓ` parameterization
   in the code and in the paper agree.
1. Fix P1 and P8–P9 (typo-level, but P1 breaks Prop. 2 outright).
2. Make the Algorithm 1 change in Part 1; re-run simulations to confirm nothing degrades.
3. Rewrite Section IV per Part 2 (drift + hitting time). **Draft exists: `section-iv-new.tex`.**
4. Fix Prop. 5 per Part 3 (damped Polyak) — self-contained, can be done in parallel.
5. Rewrite the Section VI convergence claim per Part 3.
6. Add the random-search baseline to Fig. 1.
7. Run the `p_m = 1/(N_u+N_x)` experiment on the 7x7 grid (Part 2 remark) — cheap, and it turns
   the scaling argument into evidence.
8. Write the contribution paragraph per Part 4 into abstract + intro + Algorithm 1 discussion.
9. C and D items last, fitting within 8 pages.

## Handoff note for Claude Code

When run against the paper source and the simulation repo together, the intended first task is
the BLOCKER section: read the actual `Mutation` and population-initialization code, determine how
`ℓ` is represented and mutated, and report the three diagnostic numbers. Do not start editing the
`.tex` until that is answered — the answer determines whether Section IV needs a new neighborhood
definition.

Second task: apply the Part 0 defect fixes (P1–P12), which are mechanical and independent of the
blocker.

Third task: the Part 1 algorithm change plus a re-run to confirm no regression.

`section-iv-new.tex` holds compilable draft replacements for Algorithm 1 and the convergence
results; Theorems 1–2 there are safe regardless of how the blocker resolves, Theorem 3 is not.

## Working rules for this repo

- **Symbol-before-use check** on every new formula: read the current `.tex`, confirm each symbol is
  defined before the insertion point, and report new symbols with a proposed definition site.
  New symbols introduced by this plan and not yet in the paper: `Θ`, `N(θ)`, `Θ_stall`,
  `p_el`, `q_min`, `λ_t`, `λ_min`, `N_cert`, `Δ_min`, `ℓ(h*)` (generalized inverse), `J^r_EA`,
  `h*_a`, `h*_s`. Definitions for all of these are in `section-iv-new.tex`.
- Keep the *analyzed* algorithm and the *implemented* algorithm identical. Where they must differ,
  say so explicitly in a remark rather than leaving the reader to find it.
- Track camera-ready length continuously against the 8-page cap.
