function s = embed_state(A,B,adj,K,phi_sys,K_dense)

    Acl = A + B*K;
    rhoCL = max(abs(eig(Acl)));

    sparsity = nnz(K)/numel(K);
    relNorm  = norm(K,'fro')/(norm(K_dense,'fro')+1e-9);

    rowNZ = sum(abs(K)>0,2);
    colNZ = sum(abs(K)>0,1);

    stats = [ ...
        rhoCL; sparsity; relNorm; ...
        mean(rowNZ); std(rowNZ); ...
        mean(colNZ); std(colNZ); ...
    ];

    s = [phi_sys; stats];
end
