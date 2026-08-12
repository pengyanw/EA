# CDC 2026 #642 — Revision Issue Register

Living record of (a) what has already been changed in this repo, and (b) every
open problem, indexed by **where in the paper it lives**.

Convention from here on: every proposed paper change states **section → location
→ what exactly changes**.

Last updated: 2026-08-01.

---

# Part 0 — What has already changed

## 0.1 Code

| File | Change | Blast radius |
|---|---|---|
| `EA functions/get_lqr_cost.m` | `‖P‖_F²` → `trace(P·Σ)`, `Σ = I` (matches Sec. II) | **Everything.** EA prunes far harder; vs dense 58/75% → 79/90% |
| ″ | **Marginal-stability guard**: reject `ρ ≥ 1 − 10⁻⁸`; reject `trace(P·Σ) ≤ 0` | Fixes a fatal regression the line above introduced — see the note below |
| `EA functions/ea_lqr_codesign_gershgorin.m` | `K_d` thresholded at `denseTol = 1e-3`; `max_links = nnz(K_d)` | `‖K_d‖₀`: 1250/4802 → 276/734 |
| ″ | `options.initLinks` (default `'paper'` = `ℓ₀ = ‖K_d‖₀`) | Cost curve now starts at the dense controller |
| ″ | `options.gateLinkMutation` (default `false`, i.e. `δ` applied every child) | Matches Sec. III text |
| ″ | `ea_params.mutRange` replaces hard-coded `randi([-5,5])` | none (same default) |
| ″ | Logging: `bestLinks/bestNnz/bestNa/bestNs/bestGene/repairPhase/repairIters` | none (additive) |
| `EA functions/gersgorin_stabilize_K.m` | Pure Polyak step of Eq. (37) (`etaCap` default `Inf`); dead `s_row_new` removed; non-finite guard | Repair failure rate 42.0% → 38.1%; final `J` 21.49 → 21.38 |
| `analysis_perf_bounds.m` | `costBM` from thresholded `K_d`; `K_full_all = K_dense`; IEEE-13 `K_dense_shin` thresholded | Fig. 1 normalization now matches the EA's internal one |
| ″ | Fitted exponential replaced by `Phi_predictor`; legend → `Thm.~1 prediction`; `LB_grid` saved | Fig. 1 yellow line |
| `plot_only/plot_figs.m` | Reads saved `LB_grid`/`LB_sh` instead of recomputing; `κ` labels from `shin_kappa` | removes a drifted duplicate formula |
| `run_multi_seed_unstable.m` | Vector PDF export; algorithm variant pinned explicitly; 2-panel layout (panel (c) removed); legend fixes | Fig. 2 |
| **new** `EA functions/Phi_predictor.m` | Eq. (16)–(21) link-pruning prediction | — |
| **new** `EA functions/Psi_predictor.m` | Actuator/sensor removal certificate + multi-flip expectation bound | — |
| **new** `plot_only/plot_fig2_pdf.m` | Rebuild Fig. 2 from saved data, no EA re-run | — |

Dead code, no callers, untouched: `cost_funcs/get_lqr_cost_new.m`,
`EA functions/local_block_lqr_cost.m`.

### Note — the `trace(P)` change had a fatal secondary effect (found and fixed)

On the IEEE 13-bus, `ρ(A) = 1` **exactly**. At generation 22 the EA reached a gene
(`ℓ = 338`, `N_a = 1`, `N_s = 5`, `nnz = 5`) whose closed loop had `1 − ρ = 2.1×10⁻¹⁵`.
The test `any(abs(eig(Acl)) >= 1)` evaluates `0.999999999999998 < 1` and **passes**,
after which `dlyap` solves an essentially singular system and returns a non-PSD `P`:

| | value | old objective | new objective |
|---|---|---|---|
| `trace(P)` | **−3.99×10¹⁵** | — | **looks like an excellent individual** |
| `‖P‖_F²` | +1.59×10³¹ | huge cost, discarded | — |

`‖P‖_F² ≥ 0` always, so the numerical failure previously surfaced as a *huge* cost
and was discarded. `trace(P)` can go negative, so elitism **locked onto it**:
`bestCost = −1.55×10¹³` for 129 of 150 generations.

After the guard: IEEE-13 `bestCost` 26.63 → 2.448 (normalized 0.0865), no negative
values, terminal controller `N_a = 1`, `N_s = 1` — which reproduces the paper's
"one sensor, one actuator, one communication link". The grid results were checked
and are **unaffected** (79.4% / 90.3% vs dense, unchanged).

General lesson worth a sentence in the paper: the `J_LQR = +∞` test needs a
numerical margin, not a strict `ρ < 1`. See **P-II-4**.

## 0.2 main.tex

| Section | What changed | Status |
|---|---|---|
| III, Mutation | `Unif{−d, d}` → uniform on the integers `{−d,…,d}` | **P-III-1 fixed** |
| III, Alg. 1 | `Mutation(θ,·)` → `Mutation(θ^c,·)` | **P-III-3 fixed** |
| IV, Lemma 3 | Deleted the false converse; added `rem:one_sided` | OK per review |
| IV, Def. 2(a) | `F(θ)` → `J_EA(θ)`, ties broken arbitrarily | **P-IV-11 fixed** |
| IV, Def. 2(b) | Deleted the false `\|X_t\| ≤ d`; noted why it fails | **P-IV-12 fixed** |
| IV, after Def. 2 | Added `h(ℓ)` as a function + its monotonicity | **P-IV-16 fixed** |
| IV, Prop. 1 proof | Fixed both wrong cross-references | **P-IV-14 fixed** |
| IV, Prop. 2 | Stability hypothesis made explicit; `P[K_t ∈ 𝒮]` removed from the bound | **still defective — see P-IV-6** |
| IV, Thm. 1 | `p_imp` matches Prop. 2; `P_imp := p_imp`; now a **conditional** expectation on `F_{t−1}` with a uniform `Δ_min := w_c − Φ(h*+1)`; part (i) labelled | **P-IV-7, P-IV-15 fixed**; **P-IV-3 open** |
| IV, Thm. 1(ii) | Generalized inverse `ℓ(h*)`; `ℓ_stab` floor; `limsup`; proof rewritten via Borel–Cantelli with the `E[J_EA(θ*_0)] < ∞` caveat stated | **P-IV-8 fixed** |
| IV, `eq:EA_sum` | Deterministic horizon made explicit; optional-stopping caveat added | **P-IV-9 fixed** |
| IV, closing ¶ | Rewritten: certificate is one-sided, bound is a lower bound | **contradicts Sec. VII — see P-VII-1** |
| V, `eq:lyap-decrease` | "Rearranging (25)" replaced by the defining Lyapunov equation; `V((A+BK)x)` parenthesised | **P-V-4 fixed** |
| V, `cor:ea-stab` | Added the missing `min{2d+1, ·}` | **P-V-1 partly fixed** (parent-identification and mask hypothesis still open) |
| V | New `cor:elite-child-stab` | **still defective — see P-IV-6** |

### Section VI — rewritten to "surrogate, not certificate", then corrected after review

| Location | What changed |
|---|---|
| Opening ¶ | No longer claims the repaired gain has the same sparsity and is stable; states the three outcomes honestly |
| `rem:G_empty` | New. Structural obstruction; **argument replaced** with the principal-submatrix version, valid for all `K` with no sampling |
| `rem:surrogate` | New. Numbers unified to one plant; explicit limits on what Prop. 5 claims |
| `prop:phase1_conv` | Hypothesis strengthened to `ℱ := {K ∈ 𝒦_𝒮 : R̄ ≤ γ̄} ≠ ∅`; rate replaced by the summability the proof actually gives; unproven `→ 0` deleted; missing lower bound on `‖g̃‖` stated |
| Algorithm 2 | Pure Polyak step; zero-subgradient guard; **Phase 2 deleted** (proved to be an identity map); old Phase 3 renumbered |
| After Alg. 2 | "Fallback" reframed as the operative mechanism, with the measured phase split |
| Alg. 2 → Sec. IV ¶ | Deleted the false "which (38) bounds" claim; both failing ingredients of Lemma 4 now listed; repaired case left open |
| `\rho^*` → `\bar\gamma` | 15 occurrences, resolving the collision with Lemma 1's `ρ` and the spectral radius `ρ(·)` |
| Table I | `ρ(A)`: `<1` → `=1` for all three plants (measured `1.000000000000`) |
| `:923` | Deleted a stray `\zeta` |

---

# Part 1 — Open issues, by paper section

Severity: **C** = invalidates a stated result · **M** = gap in a proof · **m** = presentation.

## Section II — Problem formulation

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-II-1** | **C** | `:61` def. of `N_c`; `:76` def. of `ℓ` | `N_c(K) = ‖𝒦(K)‖₀` is the **subsystem-block** adjacency count with `𝒦ᵢᵢ = 0`; `ℓ ∈ [1, N_uN_x]` is a **scalar entry** count. Prop. 1's identity `ΔJ_EA = w_c X_t − …` requires `ΔN_c = X_t`. Measured gap: `nnz(K_d)` vs block-level `N_c` = 276 vs 131 (5×5), 734 vs 382 (7×7); 18.1% / 13.4% of entries sit in diagonal blocks and contribute **zero** to `N_c`. Removing one entry from a still-nonzero block gives `ΔN_c = 0`, so Prop. 1 certifies moves that strictly increase cost. | Either redefine `N_c := ‖K‖₀` in Eq. (2)/(5) (code already does this; Prop. 1 becomes correct as written), or keep block-level and re-derive `Φ` and `X_t` at block resolution (code must change too). |
| P-II-2 | M | `:83–87` | `J_LQR = tr(PΣ)` but `Σ` is never given a value. Implemented as `Σ = I`. | One sentence stating `Σ = I`. |
| **P-II-4** | M | `:87` | The `J_LQR(K) = +∞` determination is stated as "when the closed loop is unstable", i.e. a strict `ρ < 1` test. Numerically this is not usable: on a plant with `ρ(A) = 1` the EA reaches closed loops with `1 − ρ ≈ 10⁻¹⁵` that pass the test but make the Lyapunov solve meaningless. | State that the cost is taken to be `+∞` whenever `ρ(A+BK) ≥ 1 − ε` for a numerical margin `ε` (code uses `10⁻⁸`), and note the cost diverges as `ρ → 1` so nothing is lost. |
| P-II-3 | M | Eq. (5) | Structural terms are unnormalized counts, so at `K_d` they are **34× (5×5) / 76× (7×7)** the performance scale (which is 1 by construction). The trade-off is decided almost entirely by structure. | Normalizing (`N_a/N_u` etc.) balances it but raises `h*` by `log‖K_d‖₀/(2\|log ρ\|)` ≈ 1.4–1.8. Derivation is unaffected — `w_c` enters in exactly one role, so it is a substitution `w_c → w_c/‖K_d‖₀`. |

## Section III — Algorithm

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-III-1** | **C** | `:134` vs `:378` vs `:825` | Mutation law written three ways: `Unif{−d, d}` (two-point), then used as `1/(2d+1)`, then written `Unif{−d,…,d}`. Under the literal `:134` reading with `d = 5`, `δ = −1` has probability zero and **Prop. 2's improvement path does not exist**. | Fix `:134` to `Unif{−d,…,d}` (integers). One line. |
| **P-III-2** | **C** | `:128` + Table I vs code `:296–297` | Table I reports `τ = 0`, at which the softmax `exp(−J_EA/τ)` is `0/0`. **The code does not use `τ = 0`**: `tau = max(std(validCosts), 1e-6)` is recomputed every generation from the population's cost spread. So the paper reports a parameter value the code never uses, and the operator is adaptive, not fixed. | Either report the adaptive rule in Table I and Section III, or pin `τ` in the code. The `≥ 1/N_p` selection bound in Prop. 2 survives either way (softmax weights are strictly positive whenever `τ > 0`). |
| **P-III-5** | **C** | `:141–148` vs code `:336–345` | **Crossover is a different operator than the paper describes.** Paper: single-point, `θ^c = [θ¹_{1:k}, θ²_{k+1:N_θ}]`, `k ~ U(1, N_θ−1)`. Code: *uniform* crossover — each mask bit is drawn independently from parent 1 w.p. `w1 = softmax(−J_1, −J_2)`, and `ℓ` is a fitness-weighted **convex combination** `round(w1 ℓ¹ + (1−w1) ℓ²)`, not an inherited value. | Rewrite the Crossover paragraph in Section III to the implemented operator. Consequences: (a) P2's objection ("single-point crossover always transmits the last coordinate of parent 2, so the offspring can never equal the elite") does not apply to uniform crossover; (b) Prop. 2's "crossover preserves `ℓ*_{t−1}` with probability 1" is false — it holds only when both parents share that `ℓ`, or via the bypass in P-III-6; (c) the convex combination means `ℓ` can never leave the population's current `ℓ`-range by crossover, which together with `\|δ\| ≤ d` gives the reachability bound `ℓ_min ≥ ℓ_0 − d G_max` (measured: the EA uses only 27 of the 750 available on 5×5, 9 of 750 on 7×7 — it is not rate-limited, it simply does not move `ℓ`). |
| P-III-6 | — | code `:332–334` | **Not a defect — a finding that removes planned work.** With probability `1 − p_c = 0.2` the code skips crossover entirely and sets `child ←` the fitter of the two selected parents. When the elite is selected, this is exactly the "dedicated elite-mutation offspring" the revision plan proposed adding. | The Part 1 algorithm change is unnecessary; describe the existing bypass in Algorithm 1 instead. Gives the improvement path `p ≥ (1−p_c)·P(elite selected)·(mutation factors)` with no crossover factor and no new operator. |
| P-III-3 | m | `:155–156` | Algorithm 1 passes `p_c` to `Crossover` but the operator description never uses it; `Mutation(θ,·)` argument should be `θ^c`. | Two-line edit; also decides whether Prop. 2's "crossover preserves `ℓ` w.p. 1" is right. |
| P-III-4 | m | `:124` | Initialization draws masks from Bernoulli(p) without stating `p` (code uses 0.99). | State it. |
| **P-0-1** | **C** | `\bibliography{cit}` | **No `.bib` file exists anywhere in the repository**, so the bibliography resolves to nothing and all 8 cite keys are unverifiable. Two are also mis-cited: `\cite[Sec.~3.1.2]{shin}` supports a max-of-convex-functions subgradient rule (shin is the spatial-decay reference), and `fazel2018global` was cited for Polyak-step subgradient methods (it is the LQR policy-gradient paper). The latter has been replaced by `polyak1987introduction`, which is equally unverifiable until the `.bib` exists. | Supply `cit.bib`. Nothing else in the register can be closed without it. |

## Section IV — Convergence

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-IV-1** | **C** | inherits P-II-1 | Prop. 1 → Lemma 3 → Thm. 1 strict positivity all rest on `ΔN_c = X_t`. | resolved by P-II-1. |
| **P-IV-2** | **C** | `:199` def. of `𝒮` | The sublevel constant `c` is **never fixed**. `L_J` is defined per sublevel set, `Φ` and `h*` depend on `L_J`. Enlarging `c` to make membership easy enlarges `h*` and shrinks the certified region. Everything plotted from `eq:plotted_bound` inherits the ambiguity. | Fix `c` explicitly (e.g. from the `(Ω,β)` bound), then `L_J`, `Φ`, `h*`, `σ_crit`, `h_stab`, `ℓ_stab` become determinate. |
| **P-IV-3** | **C** | `:212` vs `eq:p_imp` | `:212` assumes `a = 1, s = 1`. If the masks are pinned, the factor `(1−p_m)^{N_u+N_x}` ("masks not flipped") should be **1**. The section both fixes the masks and pays for not moving them. *Introduced by our edit restoring that factor.* | Decide: (i) genuinely restrict to `a = s = 1` → drop the factor, and state that this is not the implemented algorithm; or (ii) carry general masks → use the a/s certificate (Part 3). |
| **P-IV-4** | **C** | `:380–382` | Prop. 2's path claims the offspring has "identical masks", but crossover `θ^c = [θ¹_{1:k}, θ²_{k+1:N_θ}]` with `k ~ U(1,N_θ−1)` always inherits a tail from parent 2. Only trivial under `a = s = 1`. | Ties to P-IV-3 and P-III-3. |
| P-IV-5 | M | `:291–308` (Lemma 4) | Uses block-level `h_{t−1}` to bound **scalar** entry magnitudes: "all blocks within `h` are nonzero" does not imply "all entries in them are retained". Also asserts `supp(K_t) ⊂ supp(K_{t−1})`, false when `θ*_t` is a crossover offspring; and needs the segment `[K_{t−1},K_t] ⊂ 𝒮`. | Restate Lemma 4 for `K_t = K_s([ℓ*_{t−1}−X_t, a*, s*])` and fix the resolution mismatch. |
| **P-IV-6** | **C** | `cor:elite-child-stab` | Concludes `K_s(·) ∈ 𝒮` but the proof only supplies `(Ω,β)`-stability via `thm:stab-trunc`(i). Stability gives `J_LQR < ∞`, not `J_LQR ≤ c`. **Prop. 2's hypothesis is still not discharged.** *Our edit.* | Add the quantitative step `J_LQR(K_s) ≤ Ω²(‖Q‖+L²‖R‖)‖Σ‖/(1−β²)` and define `c` to be at least that (couples to P-IV-2). |
| P-IV-7 | M | `eq:EA_step` | `h_{t−1}` and `P_imp(t)` are random; LHS is a number. Type error. *Our edit did not fix this.* | Write `E[· \| F_{t−1}]` on the event `{h_{t−1} > h*, ℓ*_{t−1} > ℓ_stab}`; both are `F_{t−1}`-measurable. |
| P-IV-8 | M | Thm. 1(ii), `eq:ell_inf` | `ℓ*_∞` never defined and `ℓ*_t` is not monotone → use `limsup`. Our tail argument ("cannot persist indefinitely") only gives *infinitely many* uncertified generations; the tail claim needs *finitely many certified* ones. *Our edit.* | Use `Δ_min := w_c − Φ(h*+1) > 0` uniformly on the certified event ⟹ `Σ_t P(certified) ≤ J_EA(θ*_0)/(Δ_min p_el) < ∞` ⟹ Borel–Cantelli. Also needs `J_EA(θ*_0) < ∞`. |
| P-IV-9 | M | `eq:EA_sum` | Sum over "the first `T` certified generations" while the index runs over all generations, with a random RHS. If `T` is the count of certified generations it is a stopping time and needs Wald. | Decide deterministic vs random `T`; the hitting-time use wants the random version. |
| P-IV-10 | M | `:212` | The section drops `w_a N_a + w_s N_s` but still writes `J_EA`. The omission is justifiable (both are non-increasing under support shrinkage) but the argument is never made. Also, even at `a = 1`, `N_a` is not constant in `ℓ` — pruning can empty a row. | Name the analyzed function `J̃_EA`, or make the monotonicity argument. |
| P-IV-11 | m | `:219` | `F(θ)` undefined; should be `J_EA`. Argmin ties unbroken. | one line |
| P-IV-12 | m | Def. 2(b) | `\|X_t\| ≤ d` is false (crossover can bring an arbitrary `ℓ`). Nothing depends on it now. | delete |
| P-IV-13 | m | `:248` | Lemma 3's proof is a forward reference to Prop. 1; renders as a bare number. | merge, or make Lemma 3 a definition + remark |
| P-IV-14 | m | `:337`, `:348` | Wrong cross-references (the `J_EA` identity is the definition, not Lemma 2; `Φ(h)<w_c` follows from Eq. (18), not Def. 2). | fix refs |
| P-IV-15 | m | `:449` | A part (ii) with no labelled part (i). | label |
| P-IV-16 | m | — | `h(·)` is used as a function (`h(ℓ)`) and its monotonicity is relied on, but it is only defined per-generation in Def. 2(c). | one sentence: `Π_ℓ` is a nested family |

## Section V — Stability

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-V-1** | **C** | `cor:ea-stab`, `eq:p_stab` | `[ℓ*_t − ℓ_stab + d + 1]_+/(2d+1)` lacks `min(2d+1, ·)` → exceeds 1 whenever `ℓ*_t − ℓ_stab > d` (the normal early case). Proof also says the offspring's `ℓ` is `ℓ*_t + δ`; in fact it comes from the crossover parent drawn by Selection, not the elite. And it invokes `thm:stab-trunc`(i) without assuming `a = s = 1`. | Add the `min`; restate for the elite-child slot; add the mask hypothesis. Sec. IV uses this only for the *definition* of `ℓ_stab`, so Thm. 1 does not inherit the defect. |
| P-V-2 | M | `eq:trunc-bound` | `N_Δ(r)` counts **subsystem pairs** but bounds a Frobenius norm over **entries** — missing a block-size factor (1×2 here). `h_stab`, hence `ℓ_stab`, hence Thm. 1 inherit it. Also `K^⋆` undefined (should be `K_d`). | multiply by the block size |
| P-V-3 | M | `thm:stab-trunc`(ii) | No bound on `‖ΔK_as‖_F` is ever supplied, so part (ii) is unusable. In practice zeroing a row makes it far larger than `σ_crit`. | supply the bound, or state that the general-mask regime is uncertified |
| P-V-4 | m | `:552` | "Rearranging (25) gives the unit-decrease identity" — (26) is the *defining* Lyapunov equation, not a consequence. Also `V(A+BK_d x)` should be `V((A+BK_d)x)`. | restate |

## Section VI — Repair

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-VI-1** | **C** | Prop. 5 hypothesis | `G = {K : R̄(A+BK) < 1} = ∅` on **every plant in the paper**: half the rows of `B` are identically zero, and for those `R_i(A+BK) = R_i(A) = 1.20` (Fig. 1 plant) / `1.32` (Fig. 2 plant) for every `K`. Prop. 5's hypothesis can never hold. Cross-check: `K_d` itself has `ρ = 0.79` but `R̄ = 1.42`; along the pruning path 38/40 sampled controllers are Schur stable and **0/40** are Gershgorin stable. | State the structural necessary condition (`B(i,:) = 0` and `R_i(A) ≥ 1` ⟹ `G = ∅`) and scope Prop. 5 honestly. |
| **P-VI-2** | **C** | Algorithm 2 vs code | Code has fallback Phases 2 and 3 absent from Algorithm 2. **Phase 3 interpolates toward the *dense* `K_ref` then re-sparsifies by magnitude, changing the support** — violating `K_r ∈ 𝒦_𝒮(θ)`, the premise of all of Sec. VI. Measured over 7,669 repair calls: Algorithm 2 as written = 6.9% of exits, Phase 3 = 57.8%; of *successful* repairs ≈ **89% come from Phase 3**. | Either add the fallback to Algorithm 2 with a remark that Phase 3 leaves `𝒦_𝒮` and fires only when `G ∩ 𝒦_𝒮 = ∅` (exactly Prop. 5's excluded case), or remove it and re-run. |
| P-VI-3 | M | Sec. VI prose | "This change does not affect the … convergence properties of the original algorithm" — no argument. Lemma 4 needs `K_{t−1}` to equal `K_d` entrywise on its support; a repaired gain has moved off the `Π_ℓ(K_d)` family. **Phase 1 alone triggers this**, independent of Phase 3. | Add `‖K^r − K_s(θ)‖_F ≤ √T·‖K_s(θ)−K*‖_F/√(λ_min(2−λ_min))`, or scope Lemma 4 to unrepaired elites. |
| P-VI-4 | m | Algorithm 2 line 11 | Still shows `min(·, 0.5)`; Eq. (37) is the pure Polyak step. Code now implements pure Polyak. | delete the `min`. Note: **not** a no-op — measured 42.0% → 38.1% repair failure rate. |
| P-VI-5 | m | Prop. 4 | Subgradient convention at `μ_j(K) = 0` not stated; code's `sign(0) = 0` is a legal choice. | one sentence |
| P-VI-6 | m | Table I | `T` (max repair iterations, code 80) and `γ = 0.95` not listed. | add |

## Section VII — Simulations

| ID | Sev | Location | Problem | What a fix changes |
|---|---|---|---|---|
| **P-VII-1** | **C** | `:1086` | Claims the Eq.-(21) predictions "approximate true EA behavior quite well, particularly in later generations". **Directly contradicts our rewritten Sec. IV closing paragraph.** Numerically the certified per-generation decrease is `~1e−7` (7×7) against an observed normalized drop of order 0.5–0.75. | Rewrite to report the looseness factor honestly. |
| **P-VII-2** | **C** | `:1085` | "47–72% over dense LQR and 28–52% over diagonal LQR" — **obtained with the old `‖P‖_F²` objective**. With `trace(PΣ)`: 79–90% and 70–84%. | Restate with the new numbers. |
| P-VII-3 | M | Fig. 1 caption | The yellow line's certified region is **empty** under Def. 2(c) + the 1e-3 threshold: `h_t ≡ 1` while `h* ≈ 2.0–2.2`, so 0/150 generations are certified. (With the *unthresholded* `K_d`, `h_max = 12 > h*` and 150/150 are certified — but then the curve does not start at the dense baseline.) | Decide what the yellow line shows; see Part 3. |
| P-VII-4 | m | Footnote 1 | Says truncation baselines were "omitted for simplicity's sake", but the code plots a `κ=1` truncation line. | remove the line or restore the text |
| P-VII-5 | m | Fig. 2 caption | Panel (a)'s legend was removed; the dashed grey line at 1 now has no label. | say in the caption that it is the dense-LQR reference |

---

# Part 2 — Code issues not yet reconciled

| ID | File | Problem |
|---|---|---|
| **C-1** | `ea_lqr_codesign_gershgorin.m:176-178` | Repair writes `ℓ/a/s` back into the gene. Not in the paper (repair is defined as an *alternative cost evaluation*). The repaired numerical values are then discarded next generation by the top-`ℓ` decode, so only the overwrite survives. |
| C-2 | same, crossover | Fitness-weighted arithmetic interpolation on `ℓ` + per-bit biased uniform on masks. Paper says single-point crossover. |
| C-3 | same, selection | `tau = max(std(costs), 1e-6)` adaptive. Paper says `τ = 0`. |
| C-4 | `cost_funcs/cost_EA.m:57` | `comm_cost = nnz(K_sparse)` — element count. Matches Prop. 1's implicit use, contradicts Sec. II's stated `N_c`. See P-II-1. |
| C-5 | `cost_funcs/cost_EA.m` signature | 8th argument `max_links`, documented as "Normalization term for nnz(K)", is `~` (ignored). Suggests the original design intent was P-II-3's normalization. |

---

# Part 3 — Decisions that unblock several items at once

1. **`N_c` definition (P-II-1 / P-IV-1 / C-4).** Element count → code unchanged, Prop. 1 correct as written, Sec. II Eq. (2) changes. Block-level → `Φ` and `X_t` re-derived, code changes.
2. **Fix `c` in `𝒮` (P-IV-2).** Unblocks P-IV-6 and makes `h*` determinate, hence makes the Fig. 1 yellow line falsifiable.
3. **`a = s = 1` or not (P-IV-3 / P-IV-4).** Restricting is honest but analyses a different algorithm. Carrying general masks needs the a/s certificate below.
4. **What Fig. 1's yellow line shows (P-VII-3).** Three candidates measured on the current code:

| candidate | 5×5 | 7×7 | active generations | hypothesis holds? |
|---|---|---|---|---|
| link certificate (`Φ`, Sec. IV as written) | vacuous | vacuous | 0/150 | — |
| link certificate with `ε_ℓ` in place of `Υρ^h` | 69× loose | n/a (`ℓ` never moves) | 150/150 | yes |
| a/s certificate, `σ = σ_crit` | ~~8.8×~~ | ~~9.1×~~ | 150/150 | **NO — 0/150** |
| a/s certificate, `σ = ‖K−K_d‖_F` | 94× | 1101× | 150/150 | yes, but `L_J` invalid |
| a/s certificate, exact `‖∇J‖_F` | 260× | 1.5e4× | 126/150 | `L_J` still invalid |
| **oracle** (exact LQR increase — ceiling for *any* certificate) | **17.8×** | **196×** | — | — |

   **The `7.05×/6.96×` previously recorded here is void.** `Psi_predictor.m` bounds
   `‖∇J(K)‖ ≤ L_J‖K−K_d‖ ≤ L_J σ_crit`; measured `‖K_s(θ*_t) − K_d‖_F` is 0.79→3.19
   (5×5) and 0.98→4.70 (7×7) against `σ_crit` = 0.0196/0.0169 — the hypothesis fails in
   **0/150 generations**, violated by 163×/279×. Separately, `L_J ≈ 230` was fitted on the
   `a=s=1` pruning family; on the masked family the descent-lemma constant actually needed
   has median 675/191 and **max 3.2e4 / 3.0e7**, so no uniform-`L_J` version is sound
   either.

   **The looseness is structural, not a slack constant.** With the *exact* LQR increase
   substituted for `Ψ` (the best any certificate can do), 7×7 is still 196× loose. The
   residual sits in the probability factor: for the elite child, the flipped-bit set `F`
   must contain no bit outside the certified set `C`, or `gain(F)` cannot be lower-bounded
   — and it cannot, because `J_EA = +∞` on unstable genes, so the harm of an arbitrary flip
   is unbounded. Conditioning costs `(1−p_m)^{n−c}`; at `n=147, c≈32` that is `0.95^115 =
   2.8e−3`, two orders of magnitude on its own. This is specific to the `+∞` objective, not
   generic EA slack, and is worth saying in the paper.

5. **Which mask coordinate to build on — measured: actuators, decisively.**
   Counterfactual restore at the final elite (hold `ℓ*`, restore one mask to all-ones):

| | `J_EA` neither pruned | actuators pruned only | sensors pruned only | both (actual) |
|---|---|---|---|---|
| 5×5 | 34.07 | **12.20** (80.5% of joint) | **Inf** | 6.90 |
| 7×7 | 74.87 | **15.16** (91.0% of joint) | **Inf** | 9.27 |

   Actuator pruning carries 80–91% of the improvement *and* is independently viable; the
   final sensor mask with all actuators retained is unstable, so sensors have no standalone
   contribution to attribute.

   Submodularity diagnostics agree. With `f(S) = −J_LQR(K(S))/J_LQR(K_d)` on the retained
   set, sampling nested `S ⊂ T`, `v ∉ T` and testing `f(S+v)−f(S) ≥ f(T+v)−f(T)`:

| | valid triples | DR violations | ratio 5%ile | median | empirical `γ` |
|---|---|---|---|---|---|
| actuators 5×5 | 250/250 | 16.0% | 0.986 | 1.001 | ~~0.986~~ **see below** |
| actuators 7×7 | 250/250 | 15.2% | 0.965 | 1.001 | ~~0.965~~ **see below** |
| sensors 5×5 | **24**/250 | 29.2% | −1.100 | 1.029 | unusable |
| sensors 7×7 | **10**/250 | 10.0% | 0.703 | 1.023 | too few samples |

   Most sensor subsets are unstable (`f = −∞`), so that coordinate has no usable set-function
   structure.

   ⚠️ **The `γ` column above is superseded.** It is the 5th percentile of the *pairwise*
   diminishing-returns ratio, whose denominator is a single marginal gain — it diverges and
   changes sign when that marginal is near zero, returning 0.986 on one seed and **−16.96** on
   another for the same plant family. It was also measured on one seed only. Sampling the
   **Das–Kempe** quotient `Σ_{u∈U}[f(S+u)−f(S)] / [f(S∪U)−f(S)]` instead (aggregate
   denominator, far better conditioned), over 150 pairs per plant and five seeds per grid:

| | per-seed `γ_min` | worst per-seed 5th pct | median of per-seed medians |
|---|---|---|---|
| 5×5 | 0.948, **−0.208**, 0.979, 0.995, 0.934 | **0.985** | **1.0000** |
| 7×7 | 0.905, 0.934, 0.881, 0.906, 0.978 | **0.956** | **1.0000** |

   The per-sample quotient is **exactly 1.0000 at the median on every seed** — `f` is submodular
   to numerical precision in the typical case. One outlier (5×5 seed 10) gives a single
   negative sample, so the paper must quote a percentile, not the minimum, and say so.
   `EA functions/gap_predictor.m` implements the corrected estimator; script
   `scratchpad/probe_gamma.m`.

6. **The rate that follows, and the first bound whose magnitude matches the run budget.**
   For approximately submodular maximization with submodularity ratio `γ`, GSEMO / `(1+1)`-EA
   reaches a `(1 − e^{−γ})` approximation in expected `O(n²k)` evaluations (Friedrich &
   Neumann 2015; Qian, Yu & Zhou 2015; `γ` per Das & Kempe). With `γ = 0.97`,
   `1 − e^{−γ} = 0.621`:

   - 5×5: `N_u²k = 25²·4 = 2500` vs budget `20·150 = 3000` — **enough**;
   - 7×7: `49²·4 = 9604` vs 3000 — **3.2× short**.

   Testable prediction, confirmed on both grids. Extending `G_max` 150 → 500 at the
   theory's mutation rate `p_m = 1/n` (5 seeds):

| | budget vs `N_u²k` | gain @ `p_m = 1/n` | gain @ `p_m = 0.05` |
|---|---|---|---|
| 5×5 | 3000 vs 2500 — sufficient | **+0.00%** (6.7772 → 6.7772, identical to 4 d.p.) | +3.19% |
| 7×7 | 3000 vs 9604 — 3.2× short | **+5.15%** (8.0684 → 7.6527), 10-seed paired, 7/10, t=2.34, **p=0.044** | +11.57% |

   The bound's qualitative content is exactly this split: no headroom left where the budget
   covers `O(N_u²k)`, significant headroom where it does not. Three of the ten 7×7 seeds gained
   exactly zero — already converged inside 150 generations — which is the right mixture for a
   bound rather than a threshold. Needs **no** `L_J`, `σ_crit`, `Υ`, `ρ` — only `γ`, which is
   measurable. The full framework and proofs are in **`proof_IV.md`**.

   Two honesty caveats for the write-up: (i) 15–16% of DR inequalities are violated, so this
   is *approximately* submodular and `γ = 0.97` is the 5th percentile, not the Das–Kempe
   minimum (which is 0.915 / −3.34, the negative value indicating occasional monotonicity
   violations too); (ii) `J_EA` is not monotone in the mask — over-removal destabilizes to
   `+∞` — so the theorem must be stated on the stabilizable region under a cardinality
   constraint.

7. **Mutation-rate scaling — the same finding from the other side.** The EA-submodularity
   results are proved for standard bit mutation `p = 1/n`. Independently, a 10-seed paired
   comparison of `p_m = 1/n` against the paper's `p_m = 0.05`:

| | `(1−p_m)^n` at fixed rate | improvement | wins | t | p |
|---|---|---|---|---|---|
| 5×5 (`n=75`) | 0.0213 | +1.56% | 5/10 | 1.44 | 0.183 (n.s.) |
| 7×7 (`n=147`) | **5.31e−4** | **+11.36%** | **9/10** | 3.37 | **0.0083** |

   The predicted signature is not "scaled is better" but "the advantage grows with `n`",
   since `(1−p_m)^n` at fixed `p_m` decays 40× from `n=75` to `n=147` while `(1−c/n)^n → e^{−c}`
   is dimension-free (measured 0.3654 vs 0.3666). Null on the small grid, significant on the
   large one — which is the signature.

8. **B2 ("why an EA rather than something simpler") — the planned argument is refuted by our
   own measurement.** The revision plan justified the EA by asserting that `J_EA` "is not
   submodular and not even monotone in the support, so greedy has no approximation guarantee".
   Item 5 above measures the opposite on the actuator coordinate (`γ ≈ 0.96–1.0`, monotonicity
   violated in <0.4% of sampled pairs) — precisely the regime where greedy is near-optimal.
   `greedy_baseline.m` (new, standalone; does not touch the main line) tests it directly:
   same gene, same cost oracle, same plants, matched **evaluation** budget.

| plant | best greedy | EA @3000 evals | difference | greedy evals |
|---|---|---|---|---|
| 5×5 | **G4 6.232** | 7.044 | **−11.5%** | 652 (0.22× budget) |
| 7×7 | **G4 6.898** | 9.756 | **−29.3%** | 1589 (0.53×) |
| IEEE 13-bus | G2 **2.448** / G4 2.551 | 2.541 | −3.6% / +0.4% | 1483 / 403 |
| open-loop-unstable 5×5 (Fig. 2) | **G4 15.424** | 15.579 ± 0.268 | **−1.0%** | 494 (0.16×) |

   G4 = round-robin greedy over `{ℓ, a, s}`, each coordinate run to a 1-flip fixed point. It is
   **better and cheaper on both grids, and the margin grows with size** — the opposite of the
   claim at `main.tex:1200` that "the EA performs better compared to baseline for larger
   systems".

   **Extra budget does not rescue the EA.** Per-seed on 7×7, against the `G_max = 1000` runs
   (20 000 evaluations, 12× greedy's):

| seed | G4 (converged) | EA @3000 | EA @20 000 |
|---|---|---|---|
| 1 | **5.145** (1475 evals) | 7.646 | 5.561 |
| 10 | **7.446** (1573) | 10.437 | 8.345 |
| 15 | **8.104** (1720) | 11.185 | 8.803 |

   Greedy wins all three even at 12× the budget.

   ⚠️ **This does not affect Section IV.** Theorems A / B / C′ describe where a search stops and
   how far that is from optimal; they hold for greedy too — Theorem C in fact characterises
   exactly the fixed points of variant G1. What it affects is only the *contribution* claim.

   **The unstable plant does not rescue it either — option (ii) is closed.** That was the
   regime where the repair mechanism is the EA's claimed differentiator (~half the population
   infeasible), and the hypothesis was that a 1-flip local search would stall at the stability
   boundary. It does not: G4 reaches `J = 15.424` with architecture `N_a/N_s/N_c = 13/26/26`,
   **exactly the value and architecture of the EA's best of three seeds**, in 494 evaluations
   against 3000. It crosses the infeasible set freely — 202 of those 494 candidates scored
   `Inf`, repaired by the same `gersgorin_stabilize_K`. (The single-coordinate variant G1 does
   stall there, at 26.174; that is a property of searching one coordinate, not of greedy.)

   **Independently re-verified.** The rows above came from the subagent's `greedy_baseline.m`.
   They were re-tested with a greedy written from scratch and a cost oracle transcribed line by
   line from `ea_lqr_codesign_gershgorin.m:206–243` (`scratchpad/verify_unstable.m`,
   `scratchpad/verify_lscan.m`), varying only the `ℓ` phase:

| plant | greedy, `ℓ` local (±1/±5 — the EA's own link neighbourhood) | greedy, `ℓ` global (64-pt scan) | EA, 3000 evals |
|---|---|---|---|
| unstable 5×5 (repair on) | 16.740 (491 evals) | **15.4238** — arch. 13/26/26, 497 evals | 15.5787 ± 0.2682 |
| stable 7×7 | **5.1450** (1812) | **5.1450** (1502) | 7.6457 |

   Two things this settles:

   - **The subagent's unstable-plant number is confirmed**: 15.4238 with architecture
     `N_a/N_s/N_c = 13/26/26` in 497 evaluations, reproducing its 15.424 / 13/26/26 / 494. A
     first attempt with a weaker `ℓ` phase got 16.740 and did *not* reproduce it — the `ℓ`
     search strategy is what differs, not the rest of the method.
   - **The "EA's `±d` link mutation is the weak operator" hypothesis is FALSE.** On the stable
     7×7, local `±5` descent and a global scan return the *identical* optimum (5.1450, arch.
     2/4/4). Greedy beats the EA there by 33% **using exactly the link neighbourhood the EA
     has**. The global scan matters only on the unstable plant, whose `ℓ` landscape has a local
     minimum.

   **Route (i) is closed.** The unstable plant was the last candidate regime — repair is the
   EA's claimed differentiator there. Greedy calls the *same* `gersgorin_stabilize_K`, spends
   44% of its probes (215/491) on infinite-cost candidates versus 4–22 on stable plants, and
   still ties the EA: 15.4238 against a mean of 15.5787, exactly equalling the EA's best of
   three seeds. A tie, not a win, at 17% of the budget.

   **Route (ii) — reframe the contribution as the analysis** (Theorems A/C′, the computable
   certified gap, the `p_m = c/n` rule, the stopping criterion), report greedy as a baseline in
   Fig. 1, and keep the two claims the data does support: naive truncation is 2–4× worse
   (G3s/G3d), and decoupling the coordinates is what kills a search (G1 stalls at 11–15 while
   G4 reaches 5–7). Both of those argue for the `[ℓ,a,s]` **encoding**, which is defensible as
   the contribution; neither argues for the evolutionary **operators**, which is what the
   paper currently claims.

## Sketch of the a/s certificate (not yet in the paper)

For gene `θ` with `K = K_s(θ)` and `‖K_d − K‖_F ≤ σ_crit`:

```
Δ_a(u) := w_a + w_c·|supp(K(u,:))|                                  exact saving
Ψ_a(u) := (L_J/J_LQR(K_d))·( σ_crit·‖K(u,:)‖_F + ½‖K(u,:)‖²_F )     LQR cost bound
```

`Δ_a(u) > Ψ_a(u)` ⟹ flipping `a_u → 0` strictly decreases `J_EA`. Sensors
analogous with `w_s` and column norms. Measured: the certified set is never empty
(150/150 generations on both grids), and `c_t` falls 73→10 of 75 (5×5) and 144→15
of 147 (7×7), correctly tracking the architecture becoming irreducible.

⚠️ **The `‖K_d − K‖_F ≤ σ_crit` hypothesis above is false on both plants** (0/150
generations, violated 163×/279×), and `L_J` fitted on the `a=s=1` family is not valid
on the masked family. See Part 3 item 4. The certificate's *shape* is right — it is the
only non-empty one measured — but its constants are not, and the oracle test shows
tightening them cannot recover more than one of the three missing orders of magnitude.
Part 3 items 5–6 give the replacement route (actuators + approximate submodularity),
which needs none of these constants.

For the multi-flip bound, claim the link saving for rows only (rows have pairwise
disjoint supports, as do columns; only a row and a column can share an entry), so
per-bit gains are exactly additive, and `Ψ` is sub-additive since
`‖ΔK_S‖²_F ≤ Σ_v‖v‖²` and `√(Σx²) ≤ Σx`. Then

```
E[ΔJ_EA] ≥ (1−p_m)^(n−c') · p_m · Σ_{v∈C'} g(v),     n = N_u+N_x
```

Take the max with the single-bit bound each generation; the multi-flip term wins
early (large `c'`), the single-bit term late.

---

# Part 4 — Measured constants (current code)

| | 5×5 | 7×7 |
|---|---|---|
| `N_x, N_u` | 50, 25 | 98, 49 |
| `‖K_d‖₀` (1e-3 threshold) | 276 | 734 |
| block-level `N_c(K_d)` | 131 | 382 |
| `J_LQR(K_d) = tr(P)` | 226.8 | 447.4 |
| `Υ`, `ρ` | 1.836, 0.1421 | 1.677, 0.1518 |
| `L_J` (sublevel `5·J_LQR(K_d)`) | 235.4 | 235.3 |
| `σ_crit` | 0.0196 | 0.0169 |
| `h*` | 2.01 | 2.21 |
| `h_max` reachable | **1** | **1** |
| EA: `ℓ*` | 276 → 255 | 734 → 718 |
| EA: `N_a` | 24 → 4 | 48 → 3 |
| EA: `N_s` | 49 → 7 | 95 → 9 |
| EA vs dense | 79.4% | 90.3% |
| EA vs diagonal | 69.6% | 83.6% |

### Constants for the actuator / submodularity route (seed 1, `p_m = 0.05`, `G = 150`)

| | 5×5 | 7×7 |
|---|---|---|
| `ℓ*` trajectory | 276 → 259 (min 249) | 734 → **734** (min 725) |
| `N_c` trajectory | 258 → 12 | 692 → 16 |
| `N_a`, `N_s` | 24→4, 49→10 | 47→4, 97→11 |
| `J_EA` | 33.41 → 7.27 | 73.88 → 8.74 |
| share of drop: `w_cΔN_c` / `w_aΔN_a` / `w_sΔN_s` / LQR | 47.1 / 30.6 / 29.8 / −7.5 % | 51.9 / 26.4 / 26.4 / −4.7 % |
| actuator share by counterfactual restore | **80.5%** | **91.0%** |
| Das–Kempe `γ`, worst per-seed 5th pct over 5 seeds | **0.985** | **0.956** |
| Das–Kempe `γ`, per-seed minimum (range) | 0.934–0.995, one outlier −0.208 | 0.881–0.978 |
| Das–Kempe quotient, median | **1.0000** | **1.0000** |
| `N_u²k` (k=4) vs budget `N_p G_max` | 2500 vs 3000 | 9604 vs 3000 |
| `G` 150→500 gain at `p_m = 1/n` | **+0.00%** (exact) | **+5.15%** (p = 0.044) |
| computable optimality gap, Cor. C.1′ of `proof_IV.md` (per seed) | 2.100 / 1.450 / 0.975 (31 / 19 / 15%) | 2.200 / 2.248 / 7.199 (29 / 22 / 64%) |
| `γ_f` measured on the band the search occupies | **1.0000** | **1.0000** |
| same bound evaluated at the **greedy** fixed point | 2.023 / 1.504 / 1.410 (30 / 24 / 25%) | 1.184 / 1.906 / 2.527 (23 / 26 / 31%) |
| Thm. C hypothesis holds (generations, end) | 150 ✔ / 63 ✔ / 24 ✘ | 144 ✔ / 95 ✘ / 22 ✘ |

IEEE 13-bus, same run: `γ = 1.0000`, `N_a(S*) = 1`, `M(S*) = 0.450`, gap **0.450 = 18%** of
`J_EA(S*) = 2.448`. Fig. 1 (all three panels) now plots this bound; `Phi_predictor` is no
longer called by `analysis_perf_bounds.m`.
| `G` 150→500 gain at `p_m = 0.05` | +3.19% | +11.57% |
| effective descent-lemma constant on the masked family (median / max) | 675 / 3.2e4 | 191 / **3.0e7** |
| `‖K_s(θ*_t) − K_d‖_F` vs `σ_crit` | 0.79→3.19 vs 0.0196 | 0.98→4.70 vs 0.0169 |
| exact `‖∇J(K_s(θ*_t))‖_F` | 167 → 9638 (max 1.5e5) | 520 → 5.9e4 (max 7.4e4) |

`ℓ` is **inert**: it moves 17 of an available 750 on 5×5 and 0 on 7×7, while `N_c` collapses
258→12 and 692→16. All of that collapse is mask-driven. This closes the CLAUDE.md BLOCKER:
the paper's description of `ℓ` is accurate; `ℓ` simply is not the mechanism that produces
sparsity, so no re-parameterisation of it is needed or useful.
