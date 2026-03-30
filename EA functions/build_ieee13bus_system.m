function [sys, busInfo] = build_ieee13bus_system(cfg)
% BUILD_IEEE13BUS_SYSTEM  Discrete-time oscillator model on IEEE 13-bus feeder.
%
%   Two model modes controlled by cfg.model:
%
%   'swing' (default) — Full swing-equation with heterogeneous inertia
%       and damping, discretized via ZOH:
%       theta_dot = omega
%       omega_dot = M^{-1}( -D*omega - L*theta + B_in*u )
%
%   'shin' — Shin et al. (2023) compatible model with uniform mass,
%       no damping, Euler-forward discretization:
%       A = [I, dt*I; -dt*L, I],   B = [0; dt*eta*B_in]
%       Q = eta^2*[I 0; 0 0],      R = I
%
% Syntax:
%   [sys, busInfo] = build_ieee13bus_system()
%   [sys, busInfo] = build_ieee13bus_system(cfg)
%
% Optional cfg fields:
%   .model        - 'swing' (default) or 'shin'
%   .Ts           - Sampling time (default: 0.2 for swing, 0.005 for shin)
%   .actuatedBuses- Vector of actuated bus indices 1..13 (default: all)
%   .inertiaBase  - Base inertia constant (default: 1.0, ignored in shin)
%   .dampingBase  - Base damping constant (default: 0.5, ignored in shin)
%   .rhoTarget    - Scale A so rho(A)=rhoTarget; 0=no scaling (default: 0)
%   .eta          - Input/output scaling (shin mode only, default: 0.5)
%
% Outputs:
%   sys      - Struct with fields: A, B2, Nx, Nu, Q, R
%   busInfo  - Struct with fields: nBus, busNames, adjMtx, Laplacian,
%              actuatedBuses, nodeCoords

if nargin < 1, cfg = struct(); end
model       = getf(cfg, 'model',        'swing');
actBuses    = getf(cfg, 'actuatedBuses', 1:13);
rhoTarget   = getf(cfg, 'rhoTarget',   0);

if strcmp(model, 'shin')
    Ts       = getf(cfg, 'Ts',  0.005);
    eta      = getf(cfg, 'eta', 0.5);
else
    Ts        = getf(cfg, 'Ts',          0.2);
    inertBase = getf(cfg, 'inertiaBase', 1.0);
    dampBase  = getf(cfg, 'dampingBase', 0.5);
end

%% ==================== IEEE 13-Bus Topology ====================
% Bus index mapping (1-based):
%   1:650  2:632  3:633  4:634  5:645  6:646
%   7:671  8:680  9:684  10:611 11:652 12:692 13:675
nBus = 13;
busNames = {'650','632','633','634','645','646', ...
            '671','680','684','611','652','692','675'};

% Edges: [from, to, susceptance]
%   susceptance ~ 1/impedance, larger = stronger coupling
%   Values derived from IEEE 13-bus line configurations and lengths
edgeData = [
    1,  2,   5.0;   % 650-632,  2000ft cfg601 (3-phase, thick)
    2,  3,   8.0;   % 632-633,  500ft  cfg602
    3,  4,  12.0;   % 633-634,  transformer (strong coupling)
    2,  5,   6.0;   % 632-645,  500ft  cfg603
    5,  6,   7.0;   % 645-646,  300ft  cfg603
    2,  7,   5.0;   % 632-671,  2000ft cfg601
    7,  8,   4.0;   % 671-680,  1000ft cfg601
    7,  9,   9.0;   % 671-684,  300ft  cfg604
    9, 10,   8.0;   % 684-611,  300ft  cfg605
    9, 11,   5.5;   % 684-652,  800ft  cfg607
    7, 12,  20.0;   % 671-692,  switch (very low impedance)
   12, 13,   7.0;   % 692-675,  500ft  cfg606
];
nEdge = size(edgeData, 1);

% Build adjacency and susceptance-weighted Laplacian
adjMtx = zeros(nBus);
W      = zeros(nBus);
for e = 1:nEdge
    i = edgeData(e, 1);
    j = edgeData(e, 2);
    b = edgeData(e, 3);
    adjMtx(i,j) = 1;   adjMtx(j,i) = 1;
    W(i,j) = b;         W(j,i) = b;
end
L = diag(sum(W, 2)) - W;

%% ==================== Actuator selection ====================
Nu = numel(actBuses);
B_in = zeros(nBus, Nu);
for k = 1:Nu
    B_in(actBuses(k), k) = 1;
end

Nx = 2 * nBus;   % x = [theta(13); omega(13)]
I_n = eye(nBus);

if strcmp(model, 'shin')
    %% ============ Shin-compatible: no damping, uniform mass, Euler ============
    inertias = ones(nBus, 1);
    dampings = zeros(nBus, 1);

    A  = [I_n,     Ts * I_n; ...
          -Ts * L, I_n];
    B2 = [zeros(nBus, Nu); Ts * eta * B_in];

    C_out = [eta * I_n, zeros(nBus)];
    Q = C_out' * C_out;            % eta^2 * [I 0; 0 0]
    R = eye(Nu);
else
    %% ============ Swing-equation: heterogeneous M, D, ZOH ============
    inertias = inertBase * ones(nBus, 1);
    dampings = dampBase  * ones(nBus, 1);

    inertias(1)  = 2.0  * inertBase;   % 650: substation
    inertias(4)  = 0.8  * inertBase;   % 634: after transformer
    inertias(8)  = 0.6  * inertBase;   % 680: feeder end
    inertias(10) = 0.7  * inertBase;   % 611: feeder end
    inertias(11) = 0.7  * inertBase;   % 652: feeder end
    inertias(13) = 0.8  * inertBase;   % 675: feeder end

    dampings(1)  = 0.8  * dampBase;    % 650: heavily damped source
    dampings(8)  = 0.3  * dampBase;    % 680: lightly damped
    dampings(10) = 0.3  * dampBase;    % 611: lightly damped

    M_inv = diag(1 ./ inertias);
    D_mat = diag(dampings);

    Ac = [zeros(nBus), I_n; ...
          -M_inv * L,  -M_inv * D_mat];
    Bc = [zeros(nBus, Nu); M_inv * B_in];

    sysc = ss(Ac, Bc, eye(Nx), zeros(Nx, Nu));
    sysd = c2d(sysc, Ts, 'zoh');
    A  = sysd.A;
    B2 = sysd.B;

    Q = eye(Nx);
    R = eye(Nu);
end

%% ==================== Optional: scale spectral radius ====================
if rhoTarget > 0
    rho_orig = max(abs(eig(A)));
    A = A * (rhoTarget / rho_orig);
end

%% ==================== Node coordinates for plotting ====================
% Approximate geographic layout of IEEE 13-bus
nodeCoords = [
    0.0,  5.0;   % 1: 650
    2.0,  5.0;   % 2: 632
    3.0,  6.0;   % 3: 633
    4.0,  6.5;   % 4: 634
    3.0,  4.0;   % 5: 645
    4.0,  3.5;   % 6: 646
    5.0,  5.0;   % 7: 671
    7.0,  5.0;   % 8: 680
    5.5,  3.5;   % 9: 684
    6.5,  2.5;   % 10: 611
    5.5,  2.0;   % 11: 652
    6.0,  6.0;   % 12: 692
    7.0,  6.5;   % 13: 675
];

%% ==================== Pack output ====================
sys.A  = A;
sys.B2 = B2;
sys.Nx = Nx;
sys.Nu = Nu;
sys.Q  = Q;
sys.R  = R;

busInfo.nBus          = nBus;
busInfo.busNames      = busNames;
busInfo.adjMtx        = adjMtx;
busInfo.Laplacian      = L;
busInfo.actuatedBuses = actBuses;
busInfo.nodeCoords    = nodeCoords;
busInfo.inertias      = inertias;
busInfo.dampings      = dampings;
busInfo.edgeData      = edgeData;

end

function v = getf(s, f, d)
    if isfield(s, f), v = s.(f); else, v = d; end
end
