function result = check_stability_margin(A, B, K_sparse, K_ref)
% CHECK_STABILITY_MARGIN  Check stability margin condition for perturbed controller.
%
% Usage:
%   result = check_stability_margin(A, B, K_sparse, K_ref)
%
% Inputs:
%   A, B        - System matrices
%   K_sparse    - Current (possibly perturbed / sparse) feedback gain
%   K_ref       - Reference (nominal) stabilizing feedback gain
%
% Outputs:
%   result      - Struct with margin evaluation fields:
%       .rho_Acl        spectral radius of nominal Acl
%       .rho_pert       spectral radius of perturbed system
%       .lhs            = ||ΔK||₂ / ||K_ref||₂
%       .rhs            = (1 - ρ(Acl)) / (κ(V) σₘₐₓ(B) ||K_ref||₂)
%       .margin         = rhs - lhs
%       .satisfied      = (margin > 0)
%
% Example:
%   res = check_stability_margin(A, B, K_sparse, K_ref);
%   if res.satisfied
%       disp('Perturbation within stable margin.');
%   else
%       disp('Potential instability detected.');
%   end
%
% Author: Pengyang Wu
% ===============================================================

    % ΔK perturbation
    dK = K_sparse - K_ref;

    % Closed-loop and perturbed systems
    Acl = A + B * K_ref;
    Acl_pert = A + B * K_sparse;

    % Spectral radii
    rho_Acl = max(abs(eig(Acl)));
    rho_pert = max(abs(eig(Acl_pert)));

    % Spectral condition number of Acl eigenbasis
    [~, V] = eig(Acl);
    kappa_V = cond(V);

    % Largest singular value of B
    sigma_B = svds(B, 1);

    % Norms
    K_norm = norm(K_ref, 2);
    dK_norm = norm(dK, 2);

    % Margin evaluation
    lhs = dK_norm / K_norm;
    rhs = (1 - rho_Acl) / (kappa_V * sigma_B * K_norm);
    margin = -rhs + lhs;
    satisfied = (margin > 0);

    % Package results
    result = struct('rho_Acl', rho_Acl, ...
                    'rho_pert', rho_pert, ...
                    'lhs', lhs, ...
                    'rhs', rhs, ...
                    'margin', margin, ...
                    'satisfied', satisfied);

    % Optional print summary
    % fprintf('[check_stability_margin] ρ(Acl)=%.4f, ρ(pert)=%.4f, margin=%.4e (%s)\n', ...
    %     rho_Acl, rho_pert, margin, string(satisfied));
end
