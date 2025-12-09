clear; clc;

% --- 1. Toy system ---
A = [0.2 -0.1  0.0;
     0.1  0.25 -0.05;
     0.0  0.2   0.15];

K = [ 0    0.4   0;
      0    0.1  -0.3;
      0.2  0     0];

Adj = [1 1 0;
       1 1 1;
       0 1 1];

K_dense = [ 0   -0.5   0;
            0   -0.3  -0.4;
            0.2   0   -0.2];

I = eye(3);

% compute Gersgorin radii function
gers = @(Acl) abs(diag(Acl)) + sum(abs(Acl),2) - abs(diag(Acl));

% --- 2. Iterate Gersgorin row-transfer ---
for iter = 1:20
    Acl = A + K; 
    R = gers(Acl);

    fprintf("\nIter %d: Gersgorin radii = [%f %f %f]\n", iter, R);

    % find most dangerous row
    [~, i_star] = max(R);

    fprintf("  Most dangerous row = %d\n", i_star);

    % find offdiag columns in that row
    J = find(K(i_star,:) ~= 0 & (1:3) ~= i_star);

    if isempty(J)
        fprintf("  No offdiag to transfer. Stopping.\n");
        break;
    end
    
    best_improvement = 0;
    best_K = K;

    % try transferring each offdiag element
    for j = J
        % candidate receiving rows
        L = find(Adj(:,j));  % rows that are allowed to hold K(:,j)
        
        for ell = L'
            
            if ell == i_star, continue; end  % skip same row
            
            % Build candidate K'
            K2 = K;
            val = K2(i_star,j);
            K2(i_star,j) = 0;

            % fill using K_dense direction
            K2(ell,j) = K_dense(ell,j);

            % evaluate new Gersgorin radius
            Acl2 = A + K2;
            R2 = gers(Acl2);
            old_max = max(R);
            new_max = max(R2);

            fprintf("    Try moving (%d,%d)->(%d,%d): new max R = %.4f\n", ...
                     i_star,j, ell,j, new_max);

            if new_max < old_max
                improvement = old_max - new_max;
                if improvement > best_improvement
                    best_improvement = improvement;
                    best_K = K2;
                end
            end
        end
    end

    if best_improvement > 0
        fprintf("  >>> Applying best transfer, improvement = %.4f\n", best_improvement);
        K = best_K;
    else
        fprintf("  No improving move. Algorithm stops.\n");
        break;
    end
end

fprintf("\nFinal K = \n");
disp(K);
