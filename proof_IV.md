# proof_IV.md — Proof repository for Section IV (actuator / submodularity route)

Working document for the rewritten convergence analysis. Every statement is tagged:

- **[PROVED]** — complete argument given here, no unverified hypothesis.
- **[ASSUMPTION]** — needed, stated explicitly, with the measured violation rate.
- **[MEASURED]** — a number from the simulation, with the script that produced it.
- **[NOT CLAIMED]** — deliberately withheld, with the reason.

Companion: `REVISION_ISSUES.md` Part 3 items 4–7 holds the measurements that force this
framework; this file holds the mathematics. A compilable LaTeX form of the mathematics
alone is `paper/proof_IV.tex`; notation here matches it and `paper/main.tex`.

---

## 0. Why this framework and not the submitted one

Three measurements, in the order they close off alternatives.

**(a) The link coordinate is inert.** Over 150 generations the elite's $\ell$ moves
$276 \to 259$ (5×5) and $734 \to \mathbf{734}$ (7×7), while $N_c$ collapses $258 \to 12$ and
$692 \to 16$. All of the sparsity is produced by the masks. The submitted Section IV analyses
$\ell$ under $\mathbf{a} = \mathbf{s} = \mathbf{1}$, i.e. the one coordinate the algorithm does
not use.

**(b) The link certificate is empty, not merely loose.** $h_t \equiv 1$ against
$h^\star \approx 2.0\text{–}2.2$, so $\{h_{t-1} > h^\star\}$ occurs in **0 of 150**
generations. Theorem 1 as submitted certifies nothing about the runs in Fig. 1. Replacing the
hop-distance proxy $\Upsilon\rho^{h}$ with the exact order statistic $\varepsilon_\ell$ makes
it non-empty (150/150) but does not help, by (a).

**(c) The drift route cannot be made tight, at any constant.** Substituting the *exact* LQR
increase for the cost bound — the ceiling for any certificate — still leaves the bound
$17.8\times$ (5×5) and $196\times$ (7×7) below the measured decrease. The residual is the
probability factor $(1-p_m)^{n-c}$ with $n = N_u + N_x$, which is forced: for the elite child
the flipped-bit set must contain no bit outside the certified set, because
$J_{\mathrm{EA}} = +\infty$ on unstable genes and the harm of an arbitrary flip therefore has
no bound. At $n = 147$, $c \approx 32$ that factor is $0.95^{115} = 2.8\times10^{-3}$ on its
own.

**Conclusion.** Abandon *rate* claims obtained by counting one improvement path. Prove instead
(i) where the algorithm stops, and (ii) how far that is from optimal. Both are achievable
without $L_J$, $\sigma_{\mathrm{crit}}$, $\Upsilon$, $\rho$.

**(d) Which mask coordinate.** Counterfactual restore at the final elite: actuator pruning
alone accounts for 80.5% (5×5) and 91.0% (7×7) of the improvement and is independently
viable; the final *sensor* mask with all actuators retained is **unstable**, so sensors have
no standalone contribution. Submodularity diagnostics agree: 250/250 sampled triples valid on
the actuator coordinate, only 24/250 and 10/250 on the sensor coordinate (the rest give
$J_{\mathrm{LQR}} = \infty$). **The framework is built on actuators.**

---

## 1. Setup

Fix a generation and hold the link count $\ell$ and the sensor mask $\mathbf{s}$ at their
current values. Write $\mathcal{R} := \{1,\dots,N_u\}$ for the actuator index set and, for a
retained set $S \subseteq \mathcal{R}$, let $K(S)$ denote $\Pi_\ell(K_{\mathrm d})$ with rows
outside $S$ and columns outside $\mathbf{s}$ zeroed, so that
$K(\mathcal{R}) = K_{\mathrm s}([\ell,\mathbf{1},\mathbf{s}])$.

**Value function.**

$$
f(S) \;:=\; -\,\frac{J_{\mathrm{LQR}}(K(S))}{J_{\mathrm{LQR}}(K_{\mathrm d})}
\;\in\; [-\infty,\,0)
$$

$f$ is the negated normalised LQR cost, so larger is better and $f = -\infty$ exactly on the
destabilising sets.

**Structural cost.** Let
$m_u := \bigl|\operatorname{supp}\bigl([\Pi_\ell(K_{\mathrm d})]_{u,:}\bigr) \cap \mathbf{s}\bigr|$
be the number of communication links actuator $u$ contributes.

**Lemma 1 ($M$ is exactly modular). [PROVED]**
Distinct rows of $K$ have disjoint supports, so $N_c(K(S)) = \sum_{u \in S} m_u$ exactly. An
actuator whose row is emptied by the $\ell$-truncation carries no controller entries and is
not charged $w_a$ by $J_{\mathrm{EA}}$, so the actuator term counts only the rows that
survive:

$$
M(S) \;:=\; \sum_{u \in S}\Bigl( w_a\,\mathbf{1}\{m_u > 0\} \;+\; w_c\,m_u \Bigr)
$$

Each summand is a weight fixed by $(\ell,\mathbf{s})$ alone, so $M$ is a sum of per-element
weights, i.e. modular. $\blacksquare$

> ⚠️ **$M$ must be charged on the decoded controller, not on the mask.** Writing $w_a|S|$
> instead of $w_a\,\#\{u \in S : m_u > 0\}$ overcharges every masked-on actuator whose row the
> truncation emptied. The two coincide whenever $\ell$ is large enough that every masked row
> survives — which is why the EA never exposed the difference ($\ell$ sits at 90–99% of its
> maximum) and the greedy baseline immediately did ($\ell$ at 5–12%): at one greedy fixed
> point, 16 actuators were masked on but only 4 had entries, inflating $M(S^\star)$ from $2.0$
> to $6.8$ and the reported gap from $2.02$ to $6.82$. Fixed in
> `EA functions/gap_predictor.m`. The modularity argument is unaffected, so Theorems C / C′
> and their proofs stand as written.

> This matters: **all** of the curvature in the problem sits in $f$. No approximation is
> incurred anywhere in the structural terms.

**Objective.** With $\ell$ and $\mathbf{s}$ held fixed, the sensor and constant terms are the
same for all $S$, so minimising $J_{\mathrm{EA}}$ over the actuator mask is equivalent to
maximising

$$
U(S) \;:=\; f(S) \;-\; M(S)
\tag{IV.1}
$$

and $J_{\mathrm{EA}}(S) = -U(S) + w_s N_s + \text{const}$.

---

## 2. Structure of $U$

**Definition (submodularity ratio, Das & Kempe).** For monotone $f$,

$$
\gamma_{\mathrm f} \;:=\; \min_{S,\,T \subseteq \mathcal{R}}\;
\frac{\sum_{u \in T \setminus S}\bigl[\, f(S \cup u) - f(S) \,\bigr]}
     {f(S \cup T) - f(S)}
$$

over pairs with $f(S \cup T) > f(S)$. $f$ is submodular iff $\gamma_{\mathrm f} \geq 1$;
$\gamma_{\mathrm f} \in (0,1]$ measures how far it is from submodular.

**Assumption A1 (approximate submodularity of $f$). [ASSUMPTION]**
Let $\mathcal{B}$ be the range of architecture sizes the search occupies. For all $S$ with
$|S| \in \mathcal{B}$ and $A + BK(S)$ Schur stable, $f$ has submodularity ratio
$\gamma_{\mathrm f} > 0$.

> **[MEASURED] The band $\mathcal{B}$ is not cosmetic.** Theorem C is invoked at $S = S^\star$,
> whose size is 8–30% of $N_u$; the estimator originally sampled $[0.45 N_u,\, 0.80 N_u]$,
> a region neither search visits. Stratified sampling, 250 draws per cell, $|T| \geq 2$
> (with $|T| = 1$ the quotient is identically 1 and carries no information):
>
> | $\lvert S\rvert/N_u$ | $0.05$–$0.15$ | $0.15$–$0.30$ | $0.30$–$0.45$ | $0.45$–$0.80$ |
> |---|---|---|---|---|
> | 5×5 s1 | $1.000$ | $0.980$ | $0.981$ | $0.911$ |
> | 5×5 s10 | $1.000$ | $0.999$ | $0.871$ | $0.879$ |
> | 7×7 s1 | $0.976$ | $\mathbf{-0.119}$ | $\mathbf{-2.012}$ | $0.717$ |
> | 7×7 s10 | $1.000$ | $0.986$ | $0.652$ | $\mathbf{0.110}$ |
>
> Entries are $\min$ over samples. The searches return $|S^\star|$ of order $10\%$ of $N_u$ —
> the leftmost column, where $f$ is submodular. The negative values sit in the middle bands,
> where the aggregate gain $f(S \cup T) - f(S)$ is near zero and the quotient is dominated by
> its denominator; 18–26% of draws there are discarded for that reason. Estimating
> $\gamma_{\mathrm f}$ on the occupied band gives $\gamma_{\mathrm f} = 1.0000$ on all three
> plants. Script: `scratchpad/probe_gamma_strata.m`.
>
> ⚠️ **Two earlier estimators were wrong.** (i) The *pairwise* diminishing-returns ratio
> $[f(S{+}v) - f(S)]/[f(T{+}v) - f(T)]$ has a single marginal in the denominator; it diverges
> and changes sign when that marginal is near zero, returning $0.986$ on one seed and
> $\mathbf{-16.96}$ on another for the same plant family. (ii) The Das–Kempe form fixed that
> but was sampled on the wrong band, and included the vacuous $|T| = 1$ draws. Both are fixed
> in `EA functions/gap_predictor.m`, which now takes the band from the last quarter of the
> trajectory.

**Assumption A2 (monotonicity of $f$). [ASSUMPTION]**
$S \subseteq T$, both feasible $\Rightarrow f(S) \leq f(T)$: restoring an actuator row of
$K_{\mathrm d}$ does not increase the LQR cost.

> This is **not** automatic. $K(S)$ is a truncation of $K_{\mathrm d}$, not the optimal gain
> for the support $S$, so restoring a row can in principle hurt.
>
> **[MEASURED] A2 holds to 99.67% or better.** Over 600 sampled $(S,u)$ with $S$ feasible:
>
> | | violations $f(S{+}u) < f(S)$ | adding $u$ destabilised | relative size |
> |---|---|---|---|
> | 5×5 | **0 / 600 (0.00%)** | 0 | — |
> | 7×7 | **2 / 600 (0.33%)** | 0 | median $0.46\%$, max $0.76\%$ |
>
> No sampled addition destabilised, so the feasible family is empirically upward-closed as
> well. Script: `scratchpad/probe_A2.m`.
>
> **The two assumptions are unequal, and this matters for presentation.** A2 is essentially
> satisfied; A1 needs its band. A monotonicity failure would break the step
> $f(S_{\mathrm{opt}}) \leq f(S^\star \cup S_{\mathrm{opt}})$ outright and is absorbed
> nowhere, whereas departure from submodularity is exactly what $\gamma_{\mathrm f}$ is there
> to absorb, and its cost is visible and bounded in (IV.2).

**Lemma 2 (curvature is inherited from $f$). [PROVED]**
$U = f - M$ with $M$ modular, so $U$ is submodular iff $f$ is, and the submodularity ratio of
$U$'s curvature part equals that of $f$. $\blacksquare$

**Lemma 3 (equivalent removal form). [PROVED]**
For $R := \mathcal{R} \setminus S$,

$$
H(R) \;:=\; J_{\mathrm{EA}}(\mathcal{R}) - J_{\mathrm{EA}}(\mathcal{R} \setminus R)
\;=\; M(R) + \bigl[\, f(\mathcal{R} \setminus R) - f(\mathcal{R}) \,\bigr]
$$

The map $R \mapsto f(\mathcal{R} \setminus R)$ is submodular whenever $f$ is (the complement
of a submodular function is submodular), so $H$ is submodular and $H(\emptyset) = 0$.
$\blacksquare$

> Lemma 3 is the form to use if the paper prefers to speak of *removing* actuators; it is
> equivalent to (IV.1) and is recorded here so the two framings are not re-derived.

---

## 3. Theorem A — almost-sure finite-time convergence to a 1-flip local optimum

This is the qualitative backbone. **It uses no constants at all.**

**Definition (1-flip local optimum).** $S^\star \subseteq \mathcal{R}$ is a 1-flip local
optimum of $U$ if

$$
U(S^\star \cup \{u\}) \le U(S^\star)\ \ \forall\, u \notin S^\star
\qquad\text{and}\qquad
U(S^\star \setminus \{u\}) \le U(S^\star)\ \ \forall\, u \in S^\star
$$

**Proposition A.1 (elite-child hitting probability). [PROVED]**
Algorithm 1 produces, with probability $1 - p_c$ per offspring slot, a child equal to the
fitter of the two selected parents, which is then mutated without crossover
(`ea_lqr_codesign_gershgorin.m:332–334`). The elite has the strictly smallest cost, so
whenever it is selected as either parent it is the fitter one. Softmax selection with
$\tau > 0$ gives it weight at least $1/N_p$. Mutation flips exactly the bit $u$ and no other
with probability $p_m(1-p_m)^{n-1}$, and leaves $\ell$ fixed with probability $1/(2d+1)$.
Hence for every single-bit neighbour $\theta'$ of the elite,

$$
\mathbb{P}\bigl(\,\theta^{\mathrm c} = \theta' \;\big|\; \mathcal{H}_{t-1}\,\bigr)
\;\ge\; q \;:=\; \frac{(1 - p_c)\,p_m(1-p_m)^{n-1}}{N_p\,(2d+1)} \;>\; 0
$$

uniformly in $t$. $\blacksquare$

> **No algorithm change is required.** The revision plan had proposed adding a dedicated
> elite-mutation slot; the $1 - p_c$ bypass already is one. Recorded as P-III-6.

**Theorem A. [PROVED]**
Assume $J_{\mathrm{EA}}(\theta^\star_0) < \infty$. Then almost surely there is a finite
(random) $T$ such that $J_{\mathrm{EA}}(\theta^\star_t) = J_{\mathrm{EA}}(\theta^\star_T)$ for
all $t \ge T$, and the limiting mask is a 1-flip local optimum of $U$ in both the actuator and
sensor coordinates.

*Proof.* $\Theta = \{1,\dots,N_uN_x\} \times \{0,1\}^{N_u} \times \{0,1\}^{N_x}$ is finite, so
$J_{\mathrm{EA}}(\Theta)$ is a finite subset of $\mathbb{R} \cup \{\infty\}$. Elitism copies
the $n_e$ best genes forward unchanged (`ea_lqr_codesign_gershgorin.m:292`), so
$J_{\mathrm{EA}}(\theta^\star_t)$ is non-increasing. A non-increasing sequence taking values in
a finite ordered set is eventually constant; let $T$ be the first index at which it attains its
limit $J_\infty$.

For $t > T$ the elite cost is constant. By Proposition A.1, each generation independently
produces any prescribed single-bit neighbour $\theta'$ with probability at least $q > 0$. Since
$\sum_{t>T} q = \infty$ and the trials are conditionally independent given the past, the second
Borel–Cantelli lemma gives that $\theta'$ is produced infinitely often, almost surely. If
$J_{\mathrm{EA}}(\theta') < J_\infty$ for any such $\theta'$, elitism would carry it into the
next population and the elite cost would drop below $J_\infty$, contradicting the definition of
$J_\infty$. Hence $J_{\mathrm{EA}}(\theta') \ge J_\infty$ for every single-bit neighbour, which
is 1-flip local optimality. $\blacksquare$

**Remark A.2. [NOT CLAIMED]** Theorem A gives no rate. The chain on $\Theta$ is irreducible
(masks are reachable in one step; $\ell$ in at most $\lceil N_uN_x/d \rceil$ steps), so the
classical result for elitist EAs with full-support mutation would give almost-sure convergence
to the *global* optimum as $t \to \infty$. We do not state that, for two reasons: it is
asymptotic with no rate and therefore says nothing at $G_{\max} = 150$, and it holds for
**any** elitist GA with full-support mutation, so it cannot answer what is specific to
Algorithm 1.

---

## 4. Proposition B — exact characterisation of the returned architecture

**Proposition B. [PROVED]**
At the limiting elite of Theorem A, with $S^\star$ its retained actuator set,

$$
\begin{aligned}
u \in S^\star &\;\Longrightarrow\; f(S^\star) - f(S^\star \setminus u) \;\ge\; w_a + w_c m_u
&&\text{(each retained actuator pays for itself)}\\
u \notin S^\star &\;\Longrightarrow\; f(S^\star \cup u) - f(S^\star) \;\le\; w_a + w_c m_u
&&\text{(each discarded actuator does not)}
\end{aligned}
$$

*Proof.* Immediate from the definition of 1-flip local optimality applied to $U = f - M$ and
Lemma 1, which gives $M(S \cup u) - M(S) = w_a\mathbf{1}\{m_u>0\} + w_c m_u$. $\blacksquare$

> This is the precise sense in which the returned architecture is *irreducible*: it is exactly
> the set on which the marginal LQR value of an actuator equals its marginal structural price.
> It is checkable from a run — no constants — at a cost of $N_u - |S^\star|$ Lyapunov solves,
> and is the honest replacement for the submitted claim about pruning depth.

---

## 5. Theorem C — local optimality implies a computable optimality gap

This is the quantitative result. It replaces the drift bound.

**Theorem C. [PROVED under A1, A2]**
Let $S^\star$ be a 1-flip local optimum of $U$ and let $S_{\mathrm{opt}}$ be any feasible
competitor (in particular the global maximiser). Then

$$
U(S_{\mathrm{opt}}) - U(S^\star)
\;\le\; M(S^\star \setminus S_{\mathrm{opt}})
\;+\; \bigl(\gamma_{\mathrm f}^{-1} - 1\bigr)\,M(S_{\mathrm{opt}} \setminus S^\star)
\tag{IV.2}
$$

In particular, at $\gamma_{\mathrm f} = 1$ (exactly submodular),

$$
U(S_{\mathrm{opt}}) - U(S^\star) \;\le\; M(S^\star \setminus S_{\mathrm{opt}})
$$

and if $S^\star \subseteq S_{\mathrm{opt}}$ the local optimum is global.

*Proof.* By the definition of the submodularity ratio with $S = S^\star$ and
$T = S_{\mathrm{opt}} \setminus S^\star$,

$$
f(S^\star \cup S_{\mathrm{opt}}) - f(S^\star)
\;\le\; \gamma_{\mathrm f}^{-1}\!\!\!\sum_{u \in S_{\mathrm{opt}} \setminus S^\star}\!\!\!
\bigl[\,f(S^\star \cup u) - f(S^\star)\,\bigr]
$$

By 1-flip local optimality (the *addition* half only), each bracket is at most
$w_a + w_c m_u$, so the right-hand side is at most
$\gamma_{\mathrm f}^{-1} M(S_{\mathrm{opt}} \setminus S^\star)$. By A2,
$f(S_{\mathrm{opt}}) \le f(S^\star \cup S_{\mathrm{opt}})$, hence

$$
f(S_{\mathrm{opt}}) \;\le\; f(S^\star)
+ \gamma_{\mathrm f}^{-1} M(S_{\mathrm{opt}} \setminus S^\star)
$$

Subtract $M(S_{\mathrm{opt}})$ and use modularity of $M$, i.e.
$M(S_{\mathrm{opt}}) = M(S_{\mathrm{opt}} \setminus S^\star) + M(S_{\mathrm{opt}} \cap S^\star)$:

$$
U(S_{\mathrm{opt}}) \;\le\; f(S^\star)
+ \bigl(\gamma_{\mathrm f}^{-1} - 1\bigr) M(S_{\mathrm{opt}} \setminus S^\star)
- M(S_{\mathrm{opt}} \cap S^\star)
$$

Finally
$f(S^\star) = U(S^\star) + M(S^\star) = U(S^\star) + M(S^\star \cap S_{\mathrm{opt}}) + M(S^\star \setminus S_{\mathrm{opt}})$,
giving (IV.2). $\blacksquare$

> Only the *addition* half of local optimality is used. The removal half is spare; see open
> item 2.

**Corollary C.1 (computable bound). [PROVED under A1, A2]**
$S_{\mathrm{opt}}$ is unknown, but $S^\star \setminus S_{\mathrm{opt}} \subseteq S^\star$ and
$S_{\mathrm{opt}} \setminus S^\star \subseteq \mathcal{R} \setminus S^\star$, and $M$ is
non-negative and monotone, so

$$
\min_{\mathbf a} J_{\mathrm{EA}} \;\ge\; J_{\mathrm{EA}}(S^\star)
\;-\; \Bigl[\, M(S^\star)
+ \bigl(\gamma_{\mathrm f}^{-1} - 1\bigr) M(\mathcal{R} \setminus S^\star) \,\Bigr]
\tag{IV.3}
$$

Every quantity on the right is available from the run.

### 5.1 Theorem C′ — carrying the local-optimality slack, valid at every generation

Theorem C assumes $S^\star$ is a 1-flip local optimum. **Measured, the elite usually is not**:
over 150 generations it is locally optimal in 150/63/24 of them on the three 5×5 seeds and
144/95/22 on the 7×7 seeds, and at the *final* generation on 2 of 3 in each case. Requiring
the hypothesis therefore leaves the 7×7 panel with nothing to plot.

Rather than gate on it, carry the violation. For $u \notin S^\star$ define the **profit of
adding $u$**

$$
\xi_u \;:=\; \bigl[\, f(S^\star \cup u) - f(S^\star) - w_a - w_c m_u \,\bigr]_+ \;\ge\; 0,
\qquad
\Xi(T) \;:=\; \sum_{u \in T} \xi_u
$$

so $\xi_u = 0$ exactly when $u$ cannot be added profitably (in particular when adding it
destabilises, since then $f(S^\star \cup u) = -\infty$). By construction

$$
f(S^\star \cup u) - f(S^\star) \;\le\; w_a + w_c m_u + \xi_u
\qquad \text{for every } u \notin S^\star
\tag{IV.4}
$$

with no hypothesis at all.

**Theorem C′. [PROVED under A1, A2 — no local-optimality hypothesis]**
For any feasible $S_{\mathrm{opt}}$,

$$
U(S_{\mathrm{opt}}) - U(S^\star)
\;\le\; M(S^\star \setminus S_{\mathrm{opt}})
+ \bigl(\gamma_{\mathrm f}^{-1} - 1\bigr) M(S_{\mathrm{opt}} \setminus S^\star)
+ \gamma_{\mathrm f}^{-1}\,\Xi(S_{\mathrm{opt}} \setminus S^\star)
\tag{IV.5}
$$

*Proof.* Identical to Theorem C, with (IV.4) in place of the local-optimality inequality at the
step that bounds each marginal. When $S^\star$ is a 1-flip local optimum every $\xi_u = 0$ and
(IV.5) is exactly (IV.2). $\blacksquare$

**Corollary C.1′ (computable, unconditional).**

$$
\min_{\mathbf a} J_{\mathrm{EA}} \;\ge\; J_{\mathrm{EA}}(S^\star) - \mathrm{gap},
\qquad
\mathrm{gap} \;:=\; M(S^\star)
+ \bigl(\gamma_{\mathrm f}^{-1} - 1\bigr) M(\mathcal{R} \setminus S^\star)
+ \gamma_{\mathrm f}^{-1}\,\Xi(\mathcal{R} \setminus S^\star)
\tag{IV.6}
$$

Every term is available from the run, and the $\xi_u$ come from the same $N_u - |S^\star|$
Lyapunov solves the local-optimality check would need — so the slack costs nothing beyond the
check it replaces.

**Why this is the right form.** An unconverged run gets a *weaker* certificate, not a missing
one: $\Xi(\mathcal{R} \setminus S^\star)$ is exactly the price of not having converged, and it
goes to zero as the elite reaches a local optimum (which Theorem A guarantees it does, almost
surely, in finite time). The figure becomes one continuous band whose width reports both the
approximation quality and the convergence state.

**[MEASURED]** Per seed of the integrated run (`analysis_perf_bounds.m`, $\gamma_{\mathrm f}$
estimated on the occupied band). $m_u$ is computed at the **terminal** $\ell$ and $\mathbf{s}$
— using $N_c$ from generation 1 (a different sensor mask) inflates
$M(\mathcal{R}\setminus S^\star)$ by roughly $20\times$ and must not be done.

| plant / seed | $\gamma_{\mathrm f}$ | gap | $J_{\mathrm{EA}}(S^\star)$ | $\Rightarrow \min_{\mathbf a} J_{\mathrm{EA}} \ge$ | gap as % | Thm. C hypothesis |
|---|---|---|---|---|---|---|
| 5×5 s1 | $1.0000$ | $2.100$ | $6.775$ | $4.675$ | **31%** | 150/150, end ✔ |
| 5×5 s10 | $1.0000$ | $1.450$ | $7.811$ | $6.361$ | **19%** | 63/150, end ✔ |
| 5×5 s15 | $1.0000$ | $0.975$ | $6.547$ | $5.572$ | **15%** | 24/150, end ✘ |
| 7×7 s1 | $1.0000$ | $2.200$ | $7.646$ | $5.446$ | **29%** | 144/150, end ✔ |
| 7×7 s10 | $1.0000$ | $2.248$ | $10.437$ | $8.189$ | **22%** | 95/150, end ✘ |
| 7×7 s15 | $1.0000$ | $7.199$ | $11.185$ | $3.986$ | **64%** | 22/150, end ✘ |
| IEEE 13-bus | $1.0000$ | $0.450$ | $2.448$ | $1.998$ | **18%** | 149/150, end ✔ |

**$\gamma_{\mathrm f} = 1$ on every instance, so the middle term of (IV.6) vanishes and the
bound collapses to**

$$
\mathrm{gap} \;=\; M(S^\star) \;+\; \Xi(\mathcal{R} \setminus S^\star)
$$

— the structural cost of what was kept, plus what was left on the table by not converging.
Where $\Xi = 0$ the gap *is* $M(S^\star)$: 18–31%. The one outlier, 7×7 seed 15 at 64%, is
slack-dominated ($\Xi \approx 5.2$), and that seed is 1-flip locally optimal in only 22 of 150
generations — the bound correctly reports an unfinished run as a wider band rather than as a
tighter but unjustified one.

This is the first bound in this problem that is (i) computable from the run, (ii) free of
$L_J$, $\sigma_{\mathrm{crit}}$, $\Upsilon$, $\rho$, and (iii) finite and non-vacuous.

> ⚠️ **These numbers replace an earlier set produced by a bug in `analysis_perf_bounds.m`.**
> `A_all`, `B_all`, `K_full_all` were `cell(nGrid,1)` — indexed by grid but assigned inside the
> seed loop — so after the loop they held only the *last* seed's plant, which was then paired
> with each seed's elite genes. Every seed generates its own topology and actuator placement,
> so seeds 1 and 2 were certified against the wrong system; the 7×7 seed-1 gap read $196.9$
> instead of $2.854$. Fixed by indexing all three caches by `(grid, seed)`. The bug predates
> the Theorem C work — `Phi_predictor` was called the same way, so the previously published
> Fig. 1 yellow line was affected too.

**Use in Fig. 1.** $J_{\mathrm{EA}}(\theta^\star_t)$ is plotted with the lower envelope
$J_{\mathrm{EA}}(\theta^\star_t) - \mathrm{gap}(t)$; the shaded region between them is the
certified optimality gap. Read together the two curves *bracket* the optimum,

$$
J_{\mathrm{EA}}(\theta^\star_t) - \mathrm{gap}(t)
\;\;\le\;\; \min_{\mathbf a} J_{\mathrm{EA}}
\;\;\le\;\; J_{\mathrm{EA}}(\theta^\star_t)
$$

the upper bound being constructive ($\theta^\star_t$ is feasible) and the lower one certified —
the same reading as a branch-and-bound gap plot. The band is wide early (at full support
$M(S^\star)$ is nearly the whole cost) and narrows as the EA prunes: **the certificate's
strength and the architecture's sparsity move together**, the opposite of the submitted
Theorem 1, whose certified region was empty throughout.

### 5.2 What the greedy baseline says about the scope of Theorem C′

The greedy search of `EA functions/greedy_codesign.m` (Algorithm 3 in
`paper/greedy_baseline.tex`) supplies independent witnesses against which (IV.6) can be
checked. It searches the same gene, decodes it identically, and scores it through the same
oracle, so its returned cost is an upper bound on the true optimum over *all* of $\Theta$.

**[MEASURED]** Certified lower bound versus the greedy witness:

| plant / seed | Cor. C.1′ lower bound | greedy achieves | consistent? |
|---|---|---|---|
| 5×5 s1 | $4.610$ | $6.675$ | ✔ |
| 5×5 s10 | $6.353$ | $\mathbf{6.329}$ | ✘ |
| 5×5 s15 | $5.571$ | $5.693$ | ✔ |
| 7×7 s1 | $4.792$ | $5.145$ | ✔ |
| 7×7 s10 | $7.662$ | $\mathbf{7.446}$ | ✘ |
| 7×7 s15 | $3.507$ | $8.104$ | ✔ |

**The two ✘ rows are not counterexamples.** (IV.6) bounds $\min_{\mathbf a} J_{\mathrm{EA}}$ at
the elite's *current* $\ell$ and $\mathbf{s}$; the greedy solutions change both ($\ell$:
$219 \to 21$ and $742 \to 33$). They therefore lie outside the set the bound quantifies over.

**But they measure what the restriction costs.** Relaxing the two frozen coordinates buys more
than the *entire* certified actuator-only gap on those runs. The fixed-$(\ell,\mathbf{s})$
scope must be stated in the theorem and in the figure legend — not left as a technicality —
because a reader who takes the band as a bound on the global optimum will be wrong by more than
its width.

**The greedy fixed point is certifiably near-optimal, and more tightly than the EA's.**
Evaluating (IV.6) at the greedy architecture instead of the elite
(`scratchpad/probe_greedy_gap.m`):

| plant / seed | greedy $J$ | $\Xi$ | gap | $\Rightarrow \min_{\mathbf a} J_{\mathrm{EA}} \ge$ | gap as % |
|---|---|---|---|---|---|
| 5×5 s1 | $6.675$ | $0$ | $2.023$ | $4.652$ | **30%** |
| 5×5 s10 | $6.329$ | $0$ | $1.504$ | $4.826$ | **24%** |
| 5×5 s15 | $5.693$ | $0.409$ | $1.410$ | $4.283$ | **25%** |
| 7×7 s1 | $5.145$ | $0$ | $1.184$ | $3.961$ | **23%** |
| 7×7 s10 | $7.446$ | $0.358$ | $1.906$ | $5.539$ | **26%** |
| 7×7 s15 | $8.104$ | $0$ | $2.527$ | $5.577$ | **31%** |

Tight (23–31%) and uniform, against 15–64% at the EA's point, and $\Xi = 0$ on four of six —
the greedy search terminates at genuine 1-flip local optima, so the bound is carried by $M$
alone rather than by the slack term.

> The two rows with $\Xi > 0$ are an artefact: `gap_predictor`'s $\xi$ test charges
> $w_a + w_c m_u$ for adding $u$, but does not account for the extra $w_s$ that adding an
> actuator can trigger by activating a previously empty sensor column. Greedy's own oracle
> recounts $N_s$ from the decoded $K$ and correctly finds no improvement. The discrepancy makes
> the bound *conservative*, so it is safe, but the two tests should agree; see open item 8.

**The descent path never leaves the stable set. [MEASURED]** Greedy starts at $K_{\mathrm d}$,
which stabilises by construction, and accepts only strict improvements, so its incumbent is
finite at every step; measured, the accepted trajectory is feasible on 6/6 runs while 4–10
infeasible candidates are probed and rejected. Hence there exists a path from the dense
architecture to a certifiably near-optimal one along which $J_{\mathrm{EA}}$ decreases
monotonically and every intermediate closed loop is Schur stable. **This contradicts the
paper's B2 argument** that removing a link risks "an infinite jump greedy cannot back out of":
on these plants one never has to pass through the infeasible set at all.

**Third independent confirmation that $\ell$ is inert.** The EA leaves $\ell$ at 90–99% of its
maximum while greedy prunes it to 5–12%, and the decoded controllers are of comparable size
either way ($N_c$: EA 5–12, greedy 4–10) — once the masks retain a handful of actuators and
sensors, $\ell$ no longer binds. This agrees with the trajectory measurement (§0(a)) and the
cost decomposition, from a method that had no reason to reproduce them.

**Distinct local optima, not one optimum approached two ways.** Of six grid runs only one
returns identical masks; actuator-set Jaccard elsewhere is 0.00–0.29, and on 5×5 seed 15 the
two actuator sets are *disjoint* ($\{3,16\}$ vs $\{9,11\}$) at costs $6.55$ and $5.69$.
Theorems A / C′ describe where a search stops and how far that is from optimal; this says the
landscape has many such stopping points at similar cost, which is worth one sentence in the
paper.

**By-product — a stopping criterion.** The same quantities give a computable convergence test:
a run with profitable single-actuator additions remaining has not converged, `out.nViol` counts
them and `out.epsSum` measures how much is left on the table. On the 7×7 grid at
$G_{\max} = 150$ this fires on every seed — consistent with $N_u^2 k = 9604$ evaluations
against a budget of $3000$ (§6) and with the $+5.15\%$ still available at $G_{\max} = 500$. The
certificate and the budget analysis agree on the same diagnosis.

---

## 6. Time scale — what can and cannot be proved

**[NOT CLAIMED]** A tight bound on the number of generations to reach $S^\star$.

*Reason.* The same single-path obstruction as in §0(c). From Proposition A.1 the
per-generation probability of a *specific* improving flip is

$$
q \;=\; \frac{(1-p_c)\,p_m(1-p_m)^{n-1}}{N_p(2d+1)}
$$

At $p_m = 1/n$ and $n = 147$ this is $\approx 2.3\times10^{-6}$, so the expected wait for one
improving step is $\approx 4.4\times10^{5}$ generations against $G_{\max} = 150$. Multiplying
by the number of improving steps gives a bound of order $10^{7}$. Counting one path cannot
produce a useful time bound here, and this is a property of the argument, not of the algorithm.

**[MEASURED, reported as an empirical observation only]** The known bound for the idealised
$(1{+}1)$-EA / GSEMO on approximately submodular maximisation is $O(n^2 k)$ *evaluations* to a
$(1 - e^{-\gamma_{\mathrm f}})$ approximation (Friedrich & Neumann 2015; Qian, Yu & Zhou 2015).
Evaluated for our instances with $n = N_u$, $k = 4$:

| | $N_u^2 k$ | budget $N_p G_{\max}$ | $G$: 150→500 gain at $p_m = 1/n$ | significance |
|---|---|---|---|---|
| 5×5 | $2{,}500$ | $3{,}000$ — sufficient | **$+0.00\%$** ($6.7772 \to 6.7772$, identical to 4 d.p.) | exact, 5 seeds |
| 7×7 | $9{,}604$ | $3{,}000$ — $3.2\times$ short | **$+5.15\%$** ($8.0684 \to 7.6527$) | 10-seed paired, 7/10 improved, $t = 2.34$, **$p = 0.044$** |

The split is exactly where the bound places it: no headroom where the budget covers
$O(N_u^2 k)$, significant headroom where it does not. Three of the ten 7×7 seeds gained
*exactly* zero — they had already converged inside 150 generations — which is the right mixture
for a bound rather than a threshold. **This is not a theorem about Algorithm 1** — GSEMO
maintains a Pareto front over cardinalities, which our EA does not — and must be labelled as a
comparison, not a guarantee. Scripts: `scratchpad/probe_horizon.m`,
`scratchpad/probe_horizon2.m`.

---

## 7. Corollary — the mutation-rate rule

**Corollary D. [PROVED]** In Proposition A.1 the factor $p_m(1-p_m)^{n-1}$ is maximised at
$p_m = 1/n$, and for $p_m = c/n$,

$$
(1-p_m)^{n} \;\longrightarrow\; e^{-c}
\qquad (n \to \infty)
$$

so the hitting probability becomes dimension-free. At fixed $p_m$ it decays exponentially in
$n$.

> **[MEASURED]** $(1-1/n)^n = 0.3654$ at $n = 75$ and $0.3666$ at $n = 147$ — dimension-free to
> three digits. At the paper's $p_m = 0.05$, $(1-p_m)^n$ falls from $0.0213$ to
> $5.31\times10^{-4}$, a factor of 40, over the same range.
>
> **Prediction and test.** The corollary predicts not "scaled is better" but "the advantage
> grows with $n$". Ten-seed paired comparison of $p_m = 1/n$ against $p_m = 0.05$:
>
> | | improvement | wins | $t$ | $p$ |
> |---|---|---|---|---|
> | 5×5 ($n = 75$) | $+1.56\%$ | 5/10 | $1.44$ | $0.183$ (not significant) |
> | 7×7 ($n = 147$) | **$+11.36\%$** | **9/10** | $3.37$ | **$0.0083$** |
>
> Null on the small grid, significant on the large one — the predicted signature. Script:
> `scratchpad/probe_pm2.m`.

This is the item that answers the Associate Editor's second point: the analysis *produced* a
parameter rule that was then confirmed independently, rather than describing the algorithm
after the fact.

---

## 8. What is retired, and why

| Submitted result | Status | Reason |
|---|---|---|
| $\Phi(h)$, $h^\star$, Lemma 3 | **removed** | certifies 0/150 generations; $h_t \equiv 1 < h^\star$ |
| Definition 2(c) ($h_t$) | **removed** | hop distance is a proxy for a magnitude-ordered prune; only correlated |
| Lemma 4 (one-step LQR increase) | **removed** | needs $K_{t-1}$ entrywise $= K_{\mathrm d}$ on its support; also uses block-level $h$ to bound scalar entries (P-IV-5) |
| Prop. 1 (net improvement) | **removed** | consequence of Lemmas 3, 4 |
| Prop. 2 (improvement probability) | **replaced** by Prop. A.1 | crossover is uniform, not single-point (P-III-5), so the stated path is wrong; the $1-p_c$ bypass gives a cleaner one |
| Theorem 1 (drift + pruning depth) | **replaced** by Theorems A, C | see §0(c) |
| $\mathbf{a} = \mathbf{s} = \mathbf{1}$ | **dropped** | assumes away the coordinate that does all the work |
| $\Upsilon$, $\rho$, $L_J$, $\sigma_{\mathrm{crit}}$, $\mathcal{S}$ | **unused in Section IV** | $\Upsilon,\rho$ remain in Section V; $\mathcal{S}$'s undetermined $c$ (P-IV-2) stops being a blocker |
| main.tex `:212` ("the bound is quite close to the true convergence rate") | **deleted** | measured $2\times10^3$–$5\times10^3\times$ at best; falsifiable as written |

---

## 9. Symbols (repo rule: define before use)

| Symbol | Meaning | Define at |
|---|---|---|
| $\mathcal{R}$ | actuator index set $\{1,\dots,N_u\}$ | §1, start of Sec. IV |
| $S$, $S^\star$, $S_{\mathrm{opt}}$ | retained actuator set; local optimum; global optimum | §1 |
| $K(S)$ | $\Pi_\ell(K_{\mathrm d})$ masked to rows $S$, columns $\mathbf{s}$ | §1 |
| $f(S)$ | $-J_{\mathrm{LQR}}(K(S))/J_{\mathrm{LQR}}(K_{\mathrm d})$ | §1 |
| $m_u$ | links contributed by actuator $u$ | §1 |
| $M(S)$ | $\sum_{u\in S}(w_a\mathbf{1}\{m_u>0\} + w_c m_u)$, modular | Lemma 1 |
| $U(S)$ | $f(S) - M(S)$ | (IV.1) |
| $\gamma_{\mathrm f}$ | submodularity ratio of $f$ | Assumption A1 |
| $\mathcal{B}$ | architecture sizes the search occupies | Assumption A1 |
| $\xi_u$, $\Xi(T)$ | local-optimality slack | (IV.4) |
| $q$ | per-generation elite-child hitting probability | Prop. A.1 |
| $\mathcal{H}_{t-1}$ | history $\sigma$-algebra | Prop. A.1 |
| $H(R)$ | removal-form objective (optional) | Lemma 3 |

Note: $\mathcal{R}$ replaces the $\Omega$ used in earlier drafts, because $(\Omega,\beta)$ is
already the stability pair in Section V. Likewise $\gamma_{\mathrm f}$, not $\gamma$, because
Assumption 1(b) already uses $\gamma$ for $R \succeq \gamma I$.

Retired, so no longer needing definitions in Section IV: $\Phi$, $h^\star$, $h_t$, $h(\ell)$,
$X_t$, $p_{\mathrm{imp}}$, $P_{\mathrm{imp}}$, $\Delta_{\min}$, $p_{\mathrm{el}}$,
$\ell(h^\star)$, $\ell_{\mathrm{stab}}$ (the last stays in Section V).

---

## 10. Open items

1. ~~Direct measurement of the A2 monotonicity violation rate.~~ **Done** — 0/600 (5×5) and
   2/600 (7×7), violations of relative size $\le 0.76\%$, no destabilising additions.
2. ~~Sharpen Theorem C using the removal half of local optimality.~~ **Attempted; gives an
   incomparable form, not an improvement to the computable bound.** The removal half plus the
   submodular inequality $\sum_{u\in T}[f(S) - f(S\setminus u)] \le f(S) - f(S \setminus T)$
   for $T \subseteq S$ yields

   $$
   M(S^\star \setminus S_{\mathrm{opt}}) \;\le\; f(S^\star) - f(S^\star \cap S_{\mathrm{opt}})
   $$

   so (IV.2) also holds with $f(S^\star) - f(S^\star \cap S_{\mathrm{opt}})$ in place of
   $M(S^\star \setminus S_{\mathrm{opt}})$, and the smaller of the two applies. But that
   quantity still depends on $S_{\mathrm{opt}}$, and relaxing it to $f(S^\star) - f(\emptyset)
   = \infty$ is useless, so **Corollary C.1 is unchanged**. Worth one line as a remark.
3. ~~$\gamma_{\mathrm f}$ is a percentile, not the Das–Kempe minimum.~~ **Resolved** — the
   estimator was also the wrong quantity and sampled the wrong band; all three defects fixed.
4. ~~Sensor coordinate.~~ **Decided: Theorem C is stated for actuators only.** Theorem A still
   covers both coordinates. The restriction is forced: on the sensor coordinate $f$ is
   $-\infty$ on most subsets (only 10–24 of 250 sampled triples usable), so no submodularity
   ratio exists there. Figure legend and theorem statement must both say
   $\min_{\mathbf a} J_{\mathrm{EA}}$.
5. ~~Interaction with $\ell$.~~ **Decided: state the measurement and condition on $\ell$.**
   Over 150 generations $\ell$ moves $276 \to 259$ (5×5) and $734 \to 734$ (7×7), using 17 and
   0 of the 750 steps available, while $N_c$ collapses $258 \to 12$ and $692 \to 16$.
6. **Repair interaction — one or two sentences, deferred.** With Algorithm 2 active $K$ is no
   longer $K(S)$, so $f$ is not the value function being optimised. Section VI already scopes
   its claims to unrepaired elites.
7. ~~10-seed paired test of the 7×7 horizon gain.~~ **Done** — $+5.15\%$, 7/10 seeds improved,
   $p = 0.044$.
8. **Fix the $\xi$ test to charge the induced $w_s$.** `gap_predictor`'s slack test omits the
   extra $w_s$ that adding an actuator can trigger by activating an empty sensor column,
   which is why greedy shows $\Xi > 0$ on 2 of 6 runs where its own oracle correctly finds no
   improvement. The omission makes the bound conservative, so it is safe, but the two
   should agree.
