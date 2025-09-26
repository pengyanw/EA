function K = place_ccf(A, B, poles)
%PLACE_CCF Pole placement via Controllable Canonical Form for single-input systems
%   K = PLACE_CCF(A,B,poles) returns the state feedback gain K 
%   such that eigenvalues of (A+B*K) match poles.

    n = size(A,1);
    assert(all(size(B) == [n,1]), 'B must be an n×1 single-input matrix');
    assert(numel(poles) == n, 'The number of poles must match the system dimension');

    % Check controllability
    Wc = ctrb(A,B);
    assert(rank(Wc)==n, 'System is uncontrollable, cannot place poles.');

    % Build P^{-1} from last row of Wc^{-1}
    ri = inv(Wc);
    r = ri(end,:);
    Pinv = zeros(n);
    for k = 1:n
        Pinv(k,:) = r * (A^(k-1));
    end
    P = inv(Pinv);
    disp(Pinv)
    % Transform A into canonical form
    A_c = Pinv * A * P;
    disp(A_c)
    % Characteristic polynomials
    a = charpoly(A_c);
    d = poly(poles(:).');

    % Gains in canonical form
    K_tilde = d(2:end) - a(2:end);
  
    % Map back to original coordinates
    K = -flip(K_tilde) * Pinv;
end
