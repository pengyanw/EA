function K_new = apply_action_to_K(K, z, editableIdx)
% APPLY mask-based structural action (no value modification)
%
% Inputs:
%   K            : current sparse feedback gain
%   z            : actor latent action (dAct × 1)
%   editableIdx  : M×2 list of [i,j] entries that are allowed to flip
%
% Output:
%   K_new        : K after applying mask flips (value unchanged)
 K_new = K;
    M = size(editableIdx,1);   % number of editable links

    % sort z
    [~, ord] = sort(z, 'descend');

    % ensure indices never exceed 1..M
    ord = ord(ord <= M);       % *** <- 关键：过滤无效索引 ***

    if isempty(ord)
        return;
    end

    % choose top-k among valid indices
    k = min( min(length(z), M), 5 );   % you can tune 5
    chosen = ord(1:k);

    for idx = chosen'
        i = editableIdx(idx,1);
        j = editableIdx(idx,2);

        if abs(K_new(i,j)) > 0
            K_new(i,j) = 0;  
        else
            % you may want to use K_dense(i,j)
            % K_new(i,j) = K_dense(i,j);
            K_new(i,j) = K_new(i,j);
        end
    end
end
