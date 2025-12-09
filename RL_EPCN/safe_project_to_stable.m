function K_final = safe_project_to_stable(A,B,K_dense,K_prop)

    eta = 1.0;
    K_try = K_prop;

    while eta>1e-3
        Acl = A + B*K_try;
        if max(abs(eig(Acl))) < 1.1
            K_final = K_try;
            return;
        end
        eta = eta/2;
        K_try = (1-eta)*K_dense + eta*K_prop;
    end

    K_final = K_dense;
end
