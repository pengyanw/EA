# Greedy baseline — design note

**Question being answered.** "Why an EA rather than something simpler, e.g. greedy pruning?"
The paper currently *asserts* that greedy fails because `J_EA` is non-submodular, non-monotone
in the support, and `+inf` on the unstable set. This experiment tests that assertion. It is
deliberately built so that a greedy win is visible rather than hidden.

## Variants

All operate on the same gene `theta = [ell, a, s]`, start from `theta_0 = [nnz(K_d), 1, 1]`
(which decodes to `K_d` itself), and are scored by the same oracle.

| ID | Method | Why it is here |
|----|--------|----------------|
| `DENSE` | reference row, `theta_0` | shows where every descent curve starts |
| `G1` | 1-flip local search on `a` only | the actuator-selection problem in isolation; its fixed points are exactly the 1-flip actuator local optima that Theorem C / `gap_predictor.m` characterise |
| `G2` | 1-flip local search on `a` and `s` jointly | required baseline #1 — steepest-descent over the mask coordinates |
| `G3s` | `ell` scan (coarse grid + nested refinement), masks all-ones | required baseline #2 — "just threshold the LQR gain" |
| `G3d` | `ell <- ell - d` descent, `d = 5` (= EA `mutRange`), stop at first non-improving step | the textbook-greedy reading of "prune links greedily"; no patience knob, because any patience value would be chosen knowing the answer |
| `G4` | round-robin `link -> act -> sens`, each phase to its own fixed point, until a full cycle gains nothing | required baseline #3, and the **strongest reasonable greedy**: it is the only variant that searches all three gene coordinates, i.e. the only one facing the same search space as the EA. This is the row to quote head-to-head. |
| `EA` | `ea_lqr_codesign_gershgorin` with the exact `ea_params`/`opts`/`rng` of `analysis_perf_bounds.m` | reference |

Flips are evaluated in **both** directions (removal *and* addition), so a greedy fixed point is a
genuine 1-flip local optimum, not merely "nothing left to delete".

## Fairness controls

1. **One cost oracle.** Everything (greedy, EA, dense) is scored by `cost_EA(...)`, which calls
   `get_lqr_cost`. No guard bypassed, no separate code path for any method.
2. **Verbatim decode.** `K_d = -dlqr(A,B,Q,R)` thresholded at `1e-3`; `sortIdx` from
   `sort(abs(K_d(:)),'descend')`; keep top `ell`; then zero rows outside `a`, columns outside `s`.
   Copied line-for-line from the EA, not re-derived. The EA's returned `J` is re-scored through
   this oracle and a mismatch with `history.bestCost(end)` raises a warning.
3. **Evaluation budget is the currency.** EA = `20 x 150 = 3000`. Greedy counts one evaluation per
   candidate *considered*, including ones rejected by the spectral pre-check with no Lyapunov
   solve — that over-charges greedy, the conservative direction. **No memoisation** (the EA
   re-evaluates its elites every generation, so caching greedy's repeats would be an asymmetric
   discount). Incumbent cost is reported at 500 / 1000 / 3000 evaluations for every method.
4. **Instability is scored, not skipped** — a destabilising move costs a full evaluation and gets
   `J = Inf`. Skipping would refund greedy for probing the infeasible set.
5. **Repair parity.** Grids: repair off for both (mirrors `analysis_perf_bounds.m`).
   IEEE 13-bus: repair on for both. Repair calls counted separately from evaluations.
6. **Determinism.** Greedy has no seed dependence; the 5x5/7x7 spread is *plant* variation only,
   and IEEE 13-bus is seed-independent, so greedy is run once there (seed printed as `--`)
   rather than three times to manufacture a zero-variance row.
7. **1-flip audit [D11].** Every returned architecture, including the EA's, is re-tested against
   all `Nu + Nx` single flips in both directions. Its cost is reported but charged to nobody.

## What I expect, and what would falsify it

**Expectation (the paper's claim):** greedy stalls at a materially worse 1-flip local optimum,
and/or burns its budget on infinite-cost dead ends.

**What would falsify it:** `G4` reaching `J_EA` at or below the EA's, at or below 3000
evaluations. Recent measurements on this objective put the Das–Kempe submodularity ratio at
~0.96–1.0 with monotonicity violated in <0.4% of sampled pairs — i.e. squarely in the regime
where greedy carries strong guarantees. **Greedy may well win, and if it does that is the
result.** The honest write-up in that case is that the EA's value is not raw actuator selection
but (a) the joint `(ell, a, s)` search and (b) the repair mechanism, and the paper should say so.

**Residual biases, both directions.** (+greedy) it starts from the best-known feasible point and
`G3s` exploits the total order on `ell` that the EA can only walk by `+/-d`. (+EA) `G1`/`G2` hold
`ell` fixed, so joint moves are invisible to them — `G4` exists to close that, which is why `G4`
is the row that must be quoted.
