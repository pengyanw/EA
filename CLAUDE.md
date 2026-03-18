# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a **Sparse LQR Co-Design Framework** using Evolutionary Algorithms (EA). The goal is to find sparse feedback controllers K that stabilize networked dynamical systems while minimizing both LQR performance cost and hardware complexity (number of communication links, active actuators/sensors).

All code is **MATLAB** (`.m` files). There is no build system — scripts are run directly in MATLAB.

## Running Experiments

Primary experiment scripts (run from MATLAB):

```matlab
% Multi-seed comparison: Baseline EA vs Gershgorin-repair EA
run_multi_seed_experiment

% Test on open-loop unstable systems (rho(A) > 1)
run_multi_seed_unstable

% Performance bounds: Dense LQR vs Diagonal LQR vs EA-LQR
analysis_perf_bounds

% Power grid application (IEEE 13-bus)
run_ieee13bus_ea   % in ieee13bus/ folder
```

Test/validation:
```matlab
test_scripts/test_performance_bounds
```

## External Dependencies

The repo calls external system-generation functions not defined here:
- `generate_grid_topology(gridSize, connectThresh, seed)` — generates adjacency matrix
- `generate_grid_plant(actuatedNodes, adjMtx, ...)` — returns discrete-time `(A, B, C, D)` state-space
- `plot_graph(adjMtx, nodeCoords, color)` — network visualization

These live in a shared utilities folder outside this repo. The repo itself only contains the EA co-design layer.

## Architecture

### Gene Encoding

Each individual in the population is a chromosome:
```
[diag(Q), diag(R), n_links, b_act, b_sens]
```
- `diag(Q)`, `diag(R)` — continuous weights for LQR design (base variant)
- `n_links` — target number of non-zeros in K
- `b_act`, `b_sens` — binary masks for active actuators/sensors

The Gershgorin variant fixes Q=I, R=I and only evolves `[n_links, b_act, b_sens]`.

### Fitness / Cost Function

```
J_total = alpha * (J_LQR / J_benchmark)
        + (1-alpha) * (0.05*nnz(K) + 0.4*n_act + 0.2*n_sens)
        + beta * L_gershgorin
```

- `alpha` ∈ [0,1]: tradeoff between performance and sparsity
- Implemented in `cost_funcs/cost_EA.m`
- LQR cost computed via `dlyap()` in `EA functions/get_lqr_cost.m`
- Fast surrogate (Hutchinson trace) in `cost_funcs/get_lqr_cost_new.m`

### Three EA Variants

| File | Description |
|------|-------------|
| `EA functions/ea_lqr_codesign.m` | Base: co-evolves Q, R, structure |
| `EA functions/ea_lqr_codesign_gershgorin.m` | **Primary**: Fixed Q=I R=I, Gershgorin repair on unstable individuals |
| `EA functions/ea_lqr_khop_mask.m` | Constrains K to k-hop graph neighborhoods |

### Gershgorin Stability Repair (`EA functions/gersgorin_stabilize_K.m`)

When `rho(A+BK) >= 1`, applies a 3-phase repair:
1. Scale rows/columns to satisfy Gershgorin disk conditions
2. Project back to feasible sparsity pattern
3. Verify stability; if still unstable, assign large penalty cost

This avoids eigenvalue computation during evolution by using norm-based bounds.

### Genetic Operators (`sub_gene_next_hop/`)

- **Selection**: Softmax `P(i) = exp(-J_i / tau) / Z` (temperature `tau` controls selection pressure)
- **Crossover**: Fitness-weighted BLX for continuous genes; per-bit uniform for binary masks
- **Mutation**: Gaussian for continuous, bit-flip for binary

### Key Hyperparameters

```matlab
popSize  = 20     % Population size
maxGen   = 150    % Generations
pMutate  = 0.05   % Per-gene mutation probability
pCross   = 0.8    % Crossover probability
nTop     = 10     % Elites preserved each generation
alpha    = 0..1   % Performance vs. hardware tradeoff
```

## Output & Results

- Figures saved to `figures/` as `.png`
- Workspace data saved to `results/` as `.mat`
- Cached models in `cache/` (RF models, population buffers)
- Each experiment saves convergence curves: `bestCost`, `avgCost`, `unstableCount` per generation

## Code Organization Notes

- `EA functions/` — core EA implementations (start here for algorithm changes)
- `cost_funcs/` — cost function variants (edit here to change objective)
- `sub_gene_next_hop/` — crossover/mutation operators
- `ieee13bus/` — self-contained power grid application
- `test_scripts/` — validation against known bounds
- `main script/` — legacy scripts, kept for reference only
