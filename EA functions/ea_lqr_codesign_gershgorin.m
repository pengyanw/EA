function [result, history] = ea_lqr_codesign_gershgorin(A, B, ea_params, options)
% EA_LQR_CODESIGN_GERSHGORIN  Evolutionary sparse LQR co-design (Algorithm 1)
%
% Design points:
%   1) Q, R fixed to benchmark values (dlqr computed once, not per individual)
%   2) Gene = [n_links, b_act, b_sens] -- no Q/R in chromosome
%   3) Softmax parent selection: P(i) = exp(-J_i / tau) / Z
%   4) Segment-aware crossover:
%        Integer  (n_links) : fitness-weighted interpolation
%        Binary   (b_act, b_sens) : fitness-weighted per-bit uniform
%   5) 3-phase Gershgorin hard repair during evolution
%      (gersgorin_stabilize_K applied to every unstable individual)
%
% Syntax:
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params)
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params, options)
%
% Inputs:
%   A          - State matrix (Nx x Nx)
%   B          - Input matrix (Nx x Nu)
%   ea_params  - Struct with EA parameters:
%                .popSize     - Population size (default: 20)
%                .maxGen      - Maximum generations (default: 150)
%                .pMutate     - Mutation probability (default: 0.05)
%                .pCross      - Crossover probability (default: 0.8)
%                .nTop        - Number of elite individuals (default: 10)
%                .alpha       - Performance vs hardware tradeoff (default: 0)
%   options    - Struct with optional parameters:
%                .verbose       - Display progress (default: true)
%                .Q_bm          - Fixed Q for LQR (default: eye(Nx))
%                .R_bm          - Fixed R for LQR (default: eye(Nu))
%                .gers_weight   - Gershgorin soft penalty weight in cost_EA
%                                 (default: 0). Set >0 to add soft penalty.
%                .useGersRepair - 3-phase hard repair on every unstable
%                                 individual during evolution (default: true)
%
% Outputs:
%   result     - Struct containing:
%                .K_sparse         - Best sparse controller
%                .Q_best           - Q matrix used (= Q_bm)
%                .R_best           - R matrix used (= R_bm)
%                .cost             - Best cost value
%                .gene             - Best gene [n_links, b_act, b_sens]
%                .active_actuators - Active actuator indices
%                .active_sensors   - Active sensor indices
%                .num_links        - Number of non-zero elements in K
%                .costBM           - Benchmark cost
%                .improvement      - Relative improvement over benchmark
%   history    - Struct containing:
%                .bestCost      - Best cost per generation
%                .avgCost       - Average cost per generation
%                .unstableCount - Unstable individuals per generation
%                .deltaK_norm   - K perturbation norm per generation
%                .generations   - Generation indices
%
% Example:
%   ea_params.popSize = 20;
%   ea_params.maxGen  = 100;
%   [result, history] = ea_lqr_codesign_new(A, B, ea_params);
%
% Author: Pengyang Wu
% Date: 2026-02-10

%% ==================== Parse inputs ====================
[Nx, Nu] = size(B);

% EA parameters
if ~isfield(ea_params, 'popSize'),  ea_params.popSize = 20;    end
if ~isfield(ea_params, 'maxGen'),   ea_params.maxGen  = 150;   end
if ~isfield(ea_params, 'pMutate'),  ea_params.pMutate = 0.05;  end
if ~isfield(ea_params, 'pCross'),   ea_params.pCross  = 0.8;   end
if ~isfield(ea_params, 'nTop'),     ea_params.nTop    = 10;    end
if ~isfield(ea_params, 'alpha'),    ea_params.alpha   = 0;     end

% Options
if nargin < 4, options = struct(); end
if ~isfield(options, 'verbose'),       options.verbose       = true;  end
if ~isfield(options, 'Q_bm'),          options.Q_bm          = eye(Nx); end
if ~isfield(options, 'R_bm'),          options.R_bm          = eye(Nu); end
if ~isfield(options, 'gers_weight'),   options.gers_weight   = 0;     end
if ~isfield(options, 'useGersRepair'), options.useGersRepair = true;   end

% Gene structure: [n_links(1), b_act(Nu), b_sens(Nx)]
idx_nL   = 1;
idx_act  = 1 + (1:Nu);
idx_sens = 1 + Nu + (1:Nx);
geneLength = 1 + Nu + Nx;

%% ==================== Benchmark (computed once) ====================
Q_fixed = options.Q_bm;
R_fixed = options.R_bm;

% K_d is the dense LQR gain with |K_ij| < denseTol treated as zero. Thresholding
% here (rather than only when computing the benchmark cost) makes K_d a single
% matrix throughout: the truncation family Pi_ell(K_d), the link budget, the
% ||K_d||_0 initialization, and the plotted dense baseline all refer to the same
% support. Previously the EA truncated the raw dlqr gain (N_u N_x links) while
% the baseline used the thresholded one, so the cost curve started a factor
% ~2.4x above the dense line instead of at it.
if ~isfield(options, 'denseTol'), options.denseTol = 1e-3; end

K_dense = -dlqr(A, B, Q_fixed, R_fixed);
K_dense(abs(K_dense) <= options.denseTol) = 0;

costBM  = get_lqr_cost(A, B, Q_fixed, R_fixed, K_dense);
K_bm    = K_dense;      % retained for the result struct / summary below

% Link bounds: at most ||K_d||_0 links can be retained, since Pi_ell(K_d) draws
% from the support of K_d only.
min_links = 1;
max_links = nnz(K_dense);

% Pre-sort K_dense magnitude (shared across all individuals)
[~, K_sorted_idx] = sort(abs(K_dense(:)), 'descend');
K_rank = zeros(numel(K_dense), 1);
K_rank(K_sorted_idx) = 1:numel(K_dense);   % rank of each entry in that order

% --- How ell interacts with the masks (see Section VII discussion) ---
% 'paper'     : K <- Pi_ell(K_d), then zero the masked rows/columns. This is the
%               decode of Section III. Because the masks remove entries that
%               Pi_ell had already selected, N_c << ell in general, and ell sits
%               on a plateau of width 85-254 on which J_EA is bit-for-bit
%               constant -- selection gives ell no gradient at all.
% 'canonical' : same decode, then rewrite the gene's ell to the deepest value
%               that leaves K unchanged (the largest rank still present). The
%               phenotype, the cost, and the image of the decode map are all
%               untouched; only the redundant slack in the genotype is removed,
%               so a subsequent -delta mutation deletes real links. Modularity
%               of M at fixed (ell, s) is preserved, hence Theorem 4 still holds.
% 'masked'    : apply the masks to K_d FIRST, then keep the top-ell survivors.
%               Makes N_c = min(ell, nnz(K_d o mask)) exactly, so ell is fully
%               live. WARNING: at fixed ell the entries actuator u retains then
%               depend on which OTHER actuators are on (they compete for the
%               same budget), so m_u is no longer a function of u alone and
%               Lemma 1 (M exactly modular) fails, taking Theorem 4 with it.
if ~isfield(options, 'linkDecode'), options.linkDecode = 'paper'; end
if ~ismember(options.linkDecode, {'paper','canonical','masked'})
    error('ea_lqr_codesign_gershgorin:linkDecode', ...
          'options.linkDecode must be ''paper'', ''canonical'' or ''masked''.');
end

if options.verbose
    fprintf('========== EA-LQR Co-Design (New) ==========\n');
    fprintf('System: Nx=%d, Nu=%d\n', Nx, Nu);
    fprintf('Q, R fixed to benchmark (dlqr computed once)\n');
    fprintf('Benchmark LQR cost: %.6f\n', costBM);
    fprintf('Gene: [n_links, b_act(%d), b_sens(%d)], length=%d\n', Nu, Nx, geneLength);
    fprintf('EA: popSize=%d, maxGen=%d, pMut=%.3f, pCross=%.2f\n', ...
        ea_params.popSize, ea_params.maxGen, ea_params.pMutate, ea_params.pCross);
    fprintf('Gershgorin hard repair: %s (soft penalty weight: %.1f)\n', ...
        mat2str(options.useGersRepair), options.gers_weight);
    fprintf('Crossover: softmax selection + segment-aware\n');
    fprintf('============================================\n\n');
end

%% ==================== Initialize population ====================
if options.verbose, fprintf('Initializing population...\n'); end

% 'paper'  : ell = ||K_d||_0 for every individual, masks Bernoulli (Section III).
% 'random' : ell ~ U{1, Nu*Nx} independently per individual (the previous
%            behaviour). This matters for Section IV: the certified region of
%            Theorem 1 is {h_t > h*}, which only holds while ell is large, so a
%            random ell_0 starts the population outside the region the theory
%            describes.
% Mask initialisation. 'near-dense' is what this file has always done: every
% individual gets Bernoulli(0.99) masks, i.e. ~1% of bits off, so the whole
% population starts at essentially the all-ones architecture and there is no
% mask diversity at generation 0. 'half-split' gives the first half of the
% population the exact all-ones masks and draws the second half at
% Bernoulli(1/2), so roughly half the actuators and sensors start switched off.
% 'uniform-half' draws every individual at Bernoulli(1/2), which is the only one
% of the three that matches the sentence in Section III ("for each individual
% ... draw from a Bernoulli distribution") while still starting with real mask
% diversity.
if ~isfield(options, 'initMasks'),        options.initMasks        = 'near-dense'; end
if ~ismember(options.initMasks, {'near-dense','half-split','uniform-half'})
    error('ea_lqr_codesign_gershgorin:initMasks', ...
          ['options.initMasks must be ''near-dense'', ''half-split'' ' ...
           'or ''uniform-half''.']);
end
if ~isfield(options, 'initLinks'),        options.initLinks        = 'paper'; end
if ~isfield(options, 'gateLinkMutation'), options.gateLinkMutation = false;   end
if isfield(ea_params, 'mutRange'), mutRange = ea_params.mutRange; else, mutRange = 5; end

pop = cell(ea_params.popSize, 1);
ell0 = min(max(nnz(K_dense), min_links), max_links);
for i = 1:ea_params.popSize
    switch options.initLinks
        case 'paper'
            num_links = ell0;
        case 'random'
            num_links = randi([min_links, max_links]);
        otherwise
            error('ea_lqr_codesign_gershgorin:initLinks', ...
                  'options.initLinks must be ''paper'' or ''random''.');
    end
    switch options.initMasks
        case 'uniform-half'
            pOff = 0.5;
        case 'half-split'
            pOff = 0.5 * (i > ceil(ea_params.popSize/2));   % 0 => exact all-ones
        otherwise
            pOff = 0.01;
    end
    act_bin   = double(rand(1, Nu) > pOff);
    sens_bin  = double(rand(1, Nx) > pOff);
    pop{i} = [num_links, act_bin, sens_bin];
end

% History tracking
historyBestCost = zeros(ea_params.maxGen, 1);
historyAvgCost  = zeros(ea_params.maxGen, 1);
errbuffer       = zeros(ea_params.maxGen, 1);
repairBuffer    = zeros(ea_params.maxGen, 1);
delta_K_norm    = zeros(ea_params.maxGen, 1);

% Per-generation trajectories of the best individual's structural coordinates.
% bestLinks is the gene coordinate ell*_t; bestNnz/bestNa/bestNs are measured on
% the decoded (unrepaired) controller, so bestNnz <= bestLinks in general because
% the actuator/sensor masks are applied after the top-ell truncation.
bestLinks       = zeros(ea_params.maxGen, 1);
bestNnz         = zeros(ea_params.maxGen, 1);
bestNa          = zeros(ea_params.maxGen, 1);
bestNs          = zeros(ea_params.maxGen, 1);
bestGeneHist    = zeros(ea_params.maxGen, geneLength);

% Repair outcome histogram: column iGen counts how many repair calls that
% generation exited in phase 0/1/2/3 or failed (-1). Phases 2 and 3 are fallbacks
% that do not appear in Algorithm 2 of the paper; phase 3 additionally leaves
% K_S(theta) because it re-sparsifies against the dense reference.
repairPhase     = zeros(5, ea_params.maxGen);   % rows: phase -1, 0, 1, 2, 3
repairIters     = zeros(ea_params.maxGen, 1);   % total subgradient iterations

% Reuse the elites' costs instead of recomputing them. The n_e elites are carried
% into the next generation with their genes unchanged, and the cost is a
% deterministic function of the gene (the repair is deterministic and the
% canonicalisation is idempotent), so the recomputation returns the same number.
% It is therefore exact, not an approximation: the search trajectory is
% bit-identical and only the evaluation count and the runtime change. From
% generation 2 onwards it removes n_e of the N_p solves per generation.
% NOTE for the greedy comparison of Section VII: the "no memoisation" rule
% imposed there is justified by the EA paying for these re-evaluations, so
% enabling this changes what a like-for-like budget axis means.
% Default: on exactly when the repair is off. With the repair on, a successful
% repair overwrites the elite's gene (lines 342-346), so the elite is NOT carried
% over unchanged, re-evaluation is not idempotent, and caching would change the
% trajectory. Measured on the unstable plant: caching shifts the final cost by
% 5.8 and the unstable count by 47. The same write-back also costs elitism its
% monotonicity there --- the best cost rises on 5 of 149 generations --- which is
% a defect of the write-back, not of the caching.
if ~isfield(options, 'cacheElites')
    options.cacheElites = ~options.useGersRepair;
end
if options.cacheElites && options.useGersRepair
    warning('ea_lqr_codesign_gershgorin:cacheElites', ...
        ['cacheElites with useGersRepair changes the search: the repair ' ...
         'rewrites elite genes, so they are not carried over unchanged.']);
end
eliteCost     = [];    % costs of pop{1..nTop}, valid from generation 2 on
eliteUnstable = [];    % their J >= 1e9 flags, so the diagnostics are unchanged
eliteRepaired = [];    % and their repair-success flags
nSkipped      = 0;     % Lyapunov solves avoided, reported when verbose

% Gershgorin hard repair options
gers_opts.maxIter   = 80;
gers_opts.targetRho = 0.95;
gers_opts.verbose   = false;
% Inf = the pure Polyak step of Eq. (37), which is what Proposition 5 proves.
% Set options.gersEtaCap = 0.5 for the previously implemented clamped step.
if isfield(options, 'gersEtaCap')
    gers_opts.etaCap = options.gersEtaCap;
else
    gers_opts.etaCap = Inf;
end

%% ==================== Main evolutionary loop ====================
if options.verbose, fprintf('Starting evolution...\n'); end

for iGen = 1:ea_params.maxGen
    costs = zeros(ea_params.popSize, 1);
    unstable_count = 0;
    repair_count   = 0;
    
    % --- Fitness Evaluation ---
    % Positions 1..nTop hold the previous generation's elites, unchanged, so
    % their costs are already known when caching is on.
    unstableFlag = false(ea_params.popSize, 1);
    repairedFlag = false(ea_params.popSize, 1);
    nCached = 0;
    if options.cacheElites && ~isempty(eliteCost)
        nCached = numel(eliteCost);
        costs(1:nCached)         = eliteCost;
        unstableFlag(1:nCached)  = eliteUnstable;
        repairedFlag(1:nCached)  = eliteRepaired;
        unstable_count = unstable_count + sum(eliteUnstable);
        repair_count   = repair_count   + sum(eliteRepaired);
        nSkipped       = nSkipped + nCached;
    end

    for i = (nCached + 1):ea_params.popSize
        gene = pop{i};

        num_links = round(gene(idx_nL));
        act_bin   = logical(gene(idx_act));
        sens_bin  = logical(gene(idx_sens));
        
        switch options.linkDecode
        case 'masked'
            % Masks first, then keep the top-ell of what survives, so every
            % unit of ell is a link that actually exists.
            K_masked = K_dense;
            K_masked(~act_bin, :)  = 0;
            K_masked(:, ~sens_bin) = 0;
            surv = find(K_masked);
            [~, ord] = sort(abs(K_masked(surv)), 'descend');
            keep_idx = surv(ord(1:min(num_links, numel(ord))));
            K_sparse = zeros(Nu, Nx);
            K_sparse(keep_idx) = K_dense(keep_idx);
        otherwise
            % Sparsify the fixed K_dense by strongest elements
            K_sparse = zeros(Nu, Nx);
            keep_idx = K_sorted_idx(1:min(num_links, end));
            K_sparse(keep_idx) = K_dense(keep_idx);

            % Apply actuator/sensor masks
            K_sparse(~act_bin, :) = 0;
            K_sparse(:, ~sens_bin) = 0;

            if strcmp(options.linkDecode, 'canonical')
                % Collapse the redundant slack in ell: move it down to the
                % deepest value that leaves K_sparse identical. Phenotype and
                % cost are unchanged; only the genotype is canonicalised, so
                % the next -delta mutation removes real links.
                %
                % Canonicalising the MASKS the same way was tried and rejected:
                % it is equally phenotype-preserving, but the dormant bits it
                % clears (18-48 sensor bits at the returned elite) act as a
                % diversity reservoir, and removing them cost 5.8% mean final
                % cost over sixteen instances (9 worse, 3 better, 4 tied).
                nzr = find(K_sparse);
                if ~isempty(nzr)
                    num_links   = max(K_rank(nzr));
                    gene(idx_nL) = num_links;
                    pop{i}       = gene;
                end
            end
        end

        % --- 3-phase hard repair for unstable individuals ---
        if options.useGersRepair
            rho_pre = max(abs(eig(A + B * K_sparse)));
            if rho_pre >= 1.0
                [K_sparse, gers_info] = gersgorin_stabilize_K( ...
                    A, B, K_sparse, K_dense, gers_opts);
                repairPhase(gers_info.phase + 2, iGen) = ...
                    repairPhase(gers_info.phase + 2, iGen) + 1;
                repairIters(iGen) = repairIters(iGen) + gers_info.iters;
                if gers_info.success
                    repair_count = repair_count + 1;
                    repairedFlag(i) = true;
                    % Sync genotype to match repaired phenotype
                    gene(idx_nL)   = nnz(K_sparse);
                    gene(idx_act)  = double(any(K_sparse, 2))';
                    gene(idx_sens) = double(any(K_sparse, 1));
                    pop{i} = gene;
                end
            end
        end
        
        % Cost evaluation (returns 1e9 if still unstable after repair)
        J = cost_EA(A, B, Q_fixed, R_fixed, K_sparse, costBM, ...
            ea_params.alpha, max_links, options.gers_weight);
        costs(i) = J;
        
        if J >= 1e9
            unstable_count = unstable_count + 1;
            unstableFlag(i) = true;
        end
    end

    errbuffer(iGen)    = unstable_count;
    repairBuffer(iGen) = repair_count;
    
    % --- Selection: sort by cost ---
    [sortedCosts, sortedIdx] = sort(costs);
    pop   = pop(sortedIdx);
    costs = costs(sortedIdx);
    unstableFlag = unstableFlag(sortedIdx);
    repairedFlag = repairedFlag(sortedIdx);

    % The first nTop entries become the next generation's elites, carried over
    % with their genes unchanged, so their evaluations can be reused.
    if options.cacheElites
        eliteCost     = costs(1:ea_params.nTop);
        eliteUnstable = unstableFlag(1:ea_params.nTop);
        eliteRepaired = repairedFlag(1:ea_params.nTop);
    end

    bestCost = sortedCosts(1);
    historyBestCost(iGen) = bestCost;
    historyAvgCost(iGen)  = mean(costs(costs < 1e9));
    
    if options.verbose
        fprintf('Gen %3d: Best=%.4f, Avg=%.4f, Unstable=%d/%d, Repaired=%d\n', ...
            iGen, bestCost, historyAvgCost(iGen), unstable_count, ea_params.popSize, repair_count);
    end
    
    % --- Delta K norm for best individual ---
    % The history block must decode exactly as the fitness loop did, or
    % bestNnz/bestNa/bestNs describe a controller that was never evaluated.
    bestGene  = pop{1};
    bg_act    = logical(bestGene(idx_act));
    bg_sens   = logical(bestGene(idx_sens));
    K_sp_best = zeros(Nu, Nx);
    if strcmp(options.linkDecode, 'masked')
        K_masked_b = K_dense;
        K_masked_b(~bg_act, :)  = 0;
        K_masked_b(:, ~bg_sens) = 0;
        surv_b = find(K_masked_b);
        [~, ord_b] = sort(abs(K_masked_b(surv_b)), 'descend');
        keep_best  = surv_b(ord_b(1:min(round(bestGene(idx_nL)), numel(ord_b))));
        K_sp_best(keep_best) = K_dense(keep_best);
    else
        keep_best = K_sorted_idx(1:min(round(bestGene(idx_nL)), end));
        K_sp_best(keep_best) = K_dense(keep_best);
        K_sp_best(~bg_act, :)  = 0;
        K_sp_best(:, ~bg_sens) = 0;
    end
    delta_K_norm(iGen) = norm(K_sp_best - K_dense, 2) / (norm(K_dense, 2) + eps);

    % --- Structural coordinates of the best individual (for the Section IV
    %     predictions). bestGeneHist keeps the full gene so the actuator/sensor
    %     removal certificate, which needs the actual masks and hence the actual
    %     row/column norms of K_s(theta), can be evaluated per generation.
    bestLinks(iGen) = round(bestGene(idx_nL));
    bestNnz(iGen)   = nnz(K_sp_best);
    bestNa(iGen)    = nnz(any(K_sp_best, 2));
    bestNs(iGen)    = nnz(any(K_sp_best, 1));
    bestGeneHist(iGen, :) = bestGene;
    
    % =========================================================
    %  Create next generation: softmax + segment-aware crossover
    % =========================================================
    newPop = cell(ea_params.popSize, 1);
    newPop(1:ea_params.nTop) = pop(1:ea_params.nTop);   % Elitism
    
    % Adaptive softmax temperature
    validCosts = costs(costs < 1e9);
    if numel(validCosts) >= 2
        tau = max(std(validCosts), 1e-6);
    else
        tau = 1.0;
    end
    
    % Softmax selection probabilities (computed once per generation)
    logits = -costs / tau;
    logits(costs >= 1e9) = -inf;
    logits = logits - max(logits);
    weights = exp(logits);
    sumW = sum(weights);
    if sumW < eps
        probs = ones(ea_params.popSize, 1) / ea_params.popSize;
    else
        probs = weights / sumW;
    end
    cumP = cumsum(probs);
    
    for i = (ea_params.nTop + 1):ea_params.popSize
        
        % --- Softmax parent selection ---
        pidx1 = find(cumP >= rand(), 1, 'first');
        pidx2 = find(cumP >= rand(), 1, 'first');
        p1 = pop{pidx1};   J1 = costs(pidx1);
        p2 = pop{pidx2};   J2 = costs(pidx2);
        
        % Fitness advantage weight: w1 -> 1 when parent1 much fitter
        if J1 < 1e9 && J2 < 1e9
            w1 = exp(-J1/tau) / (exp(-J1/tau) + exp(-J2/tau));
        else
            w1 = 0.5;
        end
        
        % --- Segment-aware crossover ---
        child = zeros(1, geneLength);
        
        if rand() >= ea_params.pCross
            if J1 <= J2, child = p1; else, child = p2; end
        else
            % Integer (n_links): fitness-weighted interpolation
            child(idx_nL) = round(w1 * p1(idx_nL) + (1 - w1) * p2(idx_nL));
            
            % Binary (b_act, b_sens): per-bit, P(parent1) = w1
            for j = [idx_act, idx_sens]
                if rand() < w1
                    child(j) = p1(j);
                else
                    child(j) = p2(j);
                end
            end
        end
        
        % --- Mutation ---
        % Section III defines the mutated gene as
        %   theta_m = [sat(ell_c + delta), a XOR m_a, s XOR m_s],
        % with delta ~ Unif{-d,...,d} applied UNCONDITIONALLY and only the mask
        % flips governed by p_m. Gating delta on p_m (options.gateLinkMutation =
        % true, the previous behaviour) leaves ell nearly immobile: at p_m = 0.05
        % only ~5% of children get any ell perturbation, and delta has zero mean,
        % so ell moves only by whatever elitism ratchets.
        if ~options.gateLinkMutation || rand < ea_params.pMutate
            child(idx_nL) = child(idx_nL) + randi([-mutRange, mutRange]);
        end
        for j = [idx_act, idx_sens]
            if rand < ea_params.pMutate
                child(j) = 1 - child(j);
            end
        end
        
        % --- Bound enforcement ---
        child(idx_nL) = min(max(round(child(idx_nL)), min_links), max_links);
        
        newPop{i} = child;
    end
    
    pop = newPop;
end

%% ==================== Extract results ====================
bestGene       = pop{1};
num_links_best = round(bestGene(idx_nL));
act_bin_best   = logical(bestGene(idx_act));
sens_bin_best  = logical(bestGene(idx_sens));

keep_idx = K_sorted_idx(1:min(num_links_best, end));
K_sparse_best = zeros(Nu, Nx);
K_sparse_best(keep_idx) = K_dense(keep_idx);
K_sparse_best(~act_bin_best, :)  = 0;
K_sparse_best(:, ~sens_bin_best) = 0;

% Final safety net: hard repair on extracted result
if options.useGersRepair
    rho_final = max(abs(eig(A + B * K_sparse_best)));
    if rho_final >= 1.0
        [K_sparse_best, ~] = gersgorin_stabilize_K( ...
            A, B, K_sparse_best, K_dense, gers_opts);
    end
end

% Pack results
result.K_sparse         = K_sparse_best;
result.Q_best           = Q_fixed;
result.R_best           = R_fixed;
result.cost             = historyBestCost(end);
result.gene             = bestGene;
result.active_actuators = find(act_bin_best);
result.active_sensors   = find(sens_bin_best);
result.num_links        = nnz(K_sparse_best);
cost_BM = get_lqr_cost(A, B, Q_fixed, R_fixed, K_bm);
result.costBM       = cost_EA(A, B, Q_fixed, R_fixed, K_bm, cost_BM, ea_params.alpha, max_links);
result.improvement  = (result.costBM - result.cost) / result.costBM;

history.bestCost      = historyBestCost;
history.avgCost       = historyAvgCost;
history.unstableCount = errbuffer;
history.repairCount   = repairBuffer;
history.deltaK_norm   = delta_K_norm;
history.generations   = 1:ea_params.maxGen;
history.bestLinks     = bestLinks;
history.bestNnz       = bestNnz;
history.bestNa        = bestNa;
history.bestNs        = bestNs;
history.bestGene      = bestGeneHist;   % maxGen x geneLength
history.repairPhase   = repairPhase;    % rows: phase -1, 0, 1, 2, 3
history.repairIters   = repairIters;

%% ==================== Summary ====================
if options.verbose
    fprintf('\n========== EA Summary ==========\n');
    fprintf('Best cost: %.6f\n', result.cost);
    fprintf('Benchmark cost: %.6f\n', result.costBM);
    fprintf('Improvement: %.2f%%\n', result.improvement * 100);
    fprintf('Active actuators: %d/%d\n', length(result.active_actuators), Nu);
    fprintf('Active sensors: %d/%d\n', length(result.active_sensors), Nx);
    fprintf('Controller links: %d/%d\n', result.num_links, Nu*Nx);
    fprintf('Max delta K/K: %.4f\n', max(delta_K_norm));
    fprintf('Best link count ell*: %d -> %d\n', bestLinks(1), bestLinks(end));
    fprintf('Best nnz(K):          %d -> %d\n', bestNnz(1), bestNnz(end));
    if options.useGersRepair
        pcount = sum(repairPhase, 2);
        ptot   = sum(pcount);
        if ptot > 0
            fprintf('Repair exits (%d calls): fail=%d (%.1f%%), p0=%d, p1=%d (%.1f%%), p2=%d (%.1f%%), p3=%d (%.1f%%)\n', ...
                ptot, pcount(1), 100*pcount(1)/ptot, pcount(2), ...
                pcount(3), 100*pcount(3)/ptot, pcount(4), 100*pcount(4)/ptot, ...
                pcount(5), 100*pcount(5)/ptot);
            fprintf('  (p1 = Algorithm 2 as written; p2/p3 are fallbacks absent from the paper, p3 leaves K_S(theta))\n');
        end
    end
    fprintf('================================\n');
end

end
