function [K_new, updatedBin, updatedCts] = pole_placement_correct_K(K, A, B, C, currPopBin, currPopCts, alpha)
    % 默认 C 为单位阵
    if nargin < 4 || isempty(C)
        C = eye(size(A));
    end

    Nx = size(A, 1);
    Nu = size(B, 2);
    A_cl = A + B * K * C;
    eig_vals = eig(A_cl);
    max_eig = max(abs(eig_vals));

    % 默认输出为原 K 的结构
    K_new = K;
    updatedBin = logical(K_new ~= 0);
    updatedCts = [sum(any(updatedBin, 2)), sum(any(updatedBin, 1))];

    if max_eig <= alpha
        return;
    end

    % ========== 是否完全可控 ==========
    if rank(ctrb(A, B)) == Nx
        eig_A = eig(A);
        eig_new = eig_A;
        unstable_idx = find(abs(eig_A) > alpha);
        eig_new(unstable_idx) = alpha * exp(1i * angle(eig_A(unstable_idx)));

        try
            K_tmp = place_CCF(A, B, eig_new);
            K_new = K_tmp;
        catch
            warning('Pole placement failed (full system). Returning original K.');
            return;
        end
    else
        % ========== 可控性分解 ==========
        [P, Abar, Bbar, ~, k] = ctrbf(A, B, C);
        Ac  = Abar(1:k, 1:k);
        Bc  = Bbar(1:k, :);

        eig_Ac = eig(Ac);
        eig_new = eig_Ac;
        if checkctrb(A,B) == 1
            eig_old =  [];
        else
        eig_old = eig(Abar((k+1):end,(k+1):end));
        end
        unstable_idx = find(real(eig_Ac) > alpha);
        eig_old 
        eig_new(unstable_idx) = real(eig_Ac)/alpha * (eig_Ac(unstable_idx));

        try
            Kc = place_CCF(Ac, Bc, eig_new); %pole placement for Ac
            Kbar = zeros(size(Bbar')); 
            Kbar(:, 1:k) = Kc;
            K_tmp = Kbar * inv(P);
            

            K_new = zeros(size(K));
            K_new(1:size(K_tmp,1), 1:size(K_tmp,2)) = K_tmp; %turn back to form of normal K 
        catch
            warning('Pole placement failed (controllable part). Returning original K.');
            return;
        end
    end
   K_new(abs(K_new)<1e-6) = 0;
    % ========== update Bin and Cts ==========
    updatedBin = logical(abs(K_new) > 1e-6);    
    updatedCts = [sum(any(updatedBin, 2)), sum(any(updatedBin, 1))];  % [#Sensors, #Actuators]
    updatedBin = reshape(updatedBin,1,[]);
end
