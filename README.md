# EA co-design of actuators, sensors and communication

Code for *An Evolutionary Algorithm for Actuator-Sensor-Communication Co-Design
in Distributed Control* (Wu and Li, CDC 2026).

Given a plant `(A, B)` and weights `Q, R`, the algorithm synthesises the dense
LQR gain `K_d` once and then searches for the best **pruning** of it, trading LQ
performance against the number of actuators, sensors and communication links:

```
J_EA(theta) = J_LQR(K_s(theta)) / J_LQR(K_d)
            + w_a N_a(K) + w_s N_s(K) + w_c N_c(K)
```

The search variable is `theta = [ell, a, s]`: a link count `ell` selecting the
`ell` largest entries of `K_d`, an actuator mask `a`, and a sensor mask `s`.

## Requirements

MATLAB with the Control System Toolbox (`dlqr`, `dlyap`).

**External dependency.** The grid plant generators `generate_grid_topology.m`
and `generate_grid_plant.m` are not in this repository; they come from the SLS
toolbox, `sls-code/matlab/shared_tools/plant_generators`. Put that directory on
the MATLAB path before running anything under `results/`.

## Entry points

| Script | Produces |
|---|---|
| `analysis_perf_bounds.m` | Fig. 1 (three plants, EA vs greedy vs baselines) and Fig. 3 (scaling), plus `plot_only/fig_data.mat` |
| `run_multi_seed_unstable.m` | Fig. 2 (Gershgorin repair on an open-loop unstable plant) |
| `plot_only/plot_figs.m` | redraws Fig. 1 from the saved `.mat` without re-running the search |
| `plot_only/plot_fig2_pdf.m` | redraws Fig. 2 likewise |

## Core functions (`EA functions/`)

| File | Role |
|---|---|
| `ea_lqr_codesign_gershgorin.m` | Algorithm 1, the EA itself |
| `gersgorin_stabilize_K.m` | Algorithm 2, the repair for unstable genes |
| `greedy_prune.m` | Algorithm 3, the monotone greedy baseline |
| `gap_predictor.m` | the certified optimality gap of Theorem 4 |
| `get_lqr_cost.m` | `trace(P Sigma)` via one discrete Lyapunov solve; `Inf` if unstable |
| `truncate.m` | the kappa-hop truncation baseline |
| `build_ieee13bus_system.m` | the IEEE 13-bus plant |

`cost_funcs/cost_EA.m` assembles `J_EA` from the LQ term and the structural
penalties.

## Options that matter

Set on the `options` struct passed to `ea_lqr_codesign_gershgorin`:

- `linkDecode` — `'paper'` decodes as Section III does; `'canonical'` (default in
  the figure scripts) additionally resets the gene's `ell` to the deepest rank
  still present, which costs nothing and restores a search gradient on `ell`.
- `cacheElites` — reuse the elites' costs instead of recomputing them. Exact
  when the repair is off, and roughly halves the runtime. Defaults to
  `~useGersRepair`, because a successful repair rewrites the elite's gene.
- `initMasks` — `'near-dense'` (Bernoulli(0.99), the historical behaviour),
  `'uniform-half'`, or `'half-split'`.
- `useGersRepair` — enable Algorithm 2.

## Reproducing the figures

```matlab
addpath(genpath(pwd));
addpath('<path to>/sls-code/matlab/shared_tools/plant_generators');
analysis_perf_bounds     % Fig. 1 and Fig. 3
run_multi_seed_unstable  % Fig. 2
```

Outputs land in `results/`; the paper sources are under `paper/`.
