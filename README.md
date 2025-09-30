# README

## 📌 Project Overview

This project implements an **Evolutionary Algorithm (EA)** for **LQR controller co-design**.  
The objective is to balance **control performance (LQR cost)** and **hardware efficiency (number of communication links)**.

The algorithm automatically tunes:

- The diagonal entries of the state-weighting matrix \(Q\)
- The diagonal entries of the control-weighting matrix \(R\)
- The number of active communication links in the sparse controller \(K\)

The result is a **sparse LQR controller** that stabilizes the system while minimizing cost.

---

## ⚙️ Inputs

### 1. System Parameters

- `gridSize`: Defines an \(N \times N\) grid of nodes
- `connectThresh`: Connectivity threshold for graph topology
- `Ts`: Sampling time
- `actDensity`: Ratio of actuated nodes
- `seed`: Random seed for reproducibility

### 2. EA Parameters

- `popSize`: Population size
- `maxGen`: Maximum generations
- `pMutate`: Mutation probability per gene
- `pCross`: Crossover probability
- `nTop`: Number of elites preserved each generation
- `alpha`: Weight between performance cost and hardware cost
- `max_Q_val`, `min_Q_val`: Bounds for diagonal entries of \(Q\)
- `max_R_val`, `min_R_val`: Bounds for diagonal entries of \(R\)
- `min_links`, `max_links`: Bounds for number of nonzero entries in \(K\)

---

## 🛠️ Key Functions

- `generate_grid_topology(gridSize, connectThresh, seed)`  
  Generates graph topology, node coordinates, and system parameters.

- `plot_graph(adjMtx, nodeCoords, color)`  
  Visualizes the network topology.

- `generate_grid_plant(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts)`  
  Builds the discrete-time state-space model.

- `dlqr(A, B, Q, R)`  
  Solves for the dense LQR controller.

- `get_lqr_cost(A, B, Q, R, K)`  
  Evaluates the closed-loop LQR cost for a given controller.

---

## 🔄 Algorithm Workflow

The optimization follows an **Evolutionary Algorithm (EA)** with an **alternating method** idea (fix \(P\), optimize \(K\); then fix \(K\), optimize \(P\)).

1. **Initialization**

   - Each chromosome encodes:  
     \([ \text{diag}(Q), \text{diag}(R), \text{num\_links} ]\)

2. **Fitness Evaluation**

   - Decode chromosome into \(Q, R, \text{num_links}\)
   - Compute dense controller \(K\_{dense}\) using LQR
   - Sparsify \(K\_{dense}\) to keep `num_links` largest entries
   - Check stability (\(\rho(A+B K\_{sparse}) < 1\))
   - If stable: compute cost  
     \[
     J = \alpha \frac{J*{\text{LQR}}}{J*{\text{BM}}} \;+\; (1-\alpha)\frac{\text{nnz}(K)}{\text{max_links}}
     \]
   - If unstable: assign a large penalty

3. **Selection and Reproduction**

   - Sort population by fitness
   - Keep `nTop` elites
   - Generate new individuals via crossover and mutation

4. **Iteration**

   - Track:
     - Best cost per generation
     - Average stable cost
     - Number of unstable individuals
     - \(\|K*{new}-K*{old}\| / \|K\_{old}\|\)

5. **Results**
   - Report best solution, generation index, and improvement over benchmark
   - Save evolution plots

---

## 📤 Outputs

### Numerical Results

- `Best cost found`: Minimum cost achieved
- `Best generation`: Generation where best solution occurred
- `Advantage over benchmark`: Relative improvement vs. dense LQR baseline

### Figures

All figures saved under `figures/`:

- `evo_bestcost_grid{gridSize}seed{seed}.png`  
  Best cost evolution
- `evo_avgcost_grid{gridSize}seed{seed}.png`  
  Average stable cost evolution
- `unstable_count_grid{gridSize}seed{seed}.png`  
  Number of unstable individuals per generation

---

## 📊 Example Run

```matlab
clear; clc; close all;
addpath(genpath(pwd));

% Configure system and EA parameters
gridSize = 5; connectThresh = 0.5; Ts = 0.2; actDensity = 1; seed = 17;
popSize = 20; maxGen = 150; alpha = 0.5;

% Run optimization
main_script;  % your provided script

% Example output:
% Benchmark LQR cost (dense controller): 123.456
% Gen 100: Best Cost=0.85, Avg Cost=1.20, Unstable=3/20
% Best cost found: 0.80 at generation 120
% Best advantage over benchmark: 35%
```
