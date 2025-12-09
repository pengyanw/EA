function L_gers = gersgorin_penalty(A, B, K)
    Acl = A + B*K;
    Nx = size(Acl,1);

    radii = zeros(Nx,1);
    for i = 1:Nx
        diag_i = abs(Acl(i,i));
        offsum = sum(abs(Acl(i,:))) - abs(Acl(i,i));
        radii(i) = diag_i + offsum;
    end

    rmax = max(radii);
    L_gers = max(rmax - 1, 0);   % hinge penalty
end
