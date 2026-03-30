function K = truncate(K, G, Is, Js, kappa, problem)
% Truncate feedback matrix K based on locality radius kappa
% K(rows,cols) = 0 if node distance(i,j) > kappa
% Is{i} — column indices of K (state indices) related to node i
% Js{j} — row indices of K (actuator indices) related to node j

    if nargin < 5, kappa = 0; end
    if nargin < 6, problem = struct('name', 'grid'); end
    if strcmp(problem.name,'2d-mesh')
       c = 0.05;
       D = round(distances(G) / c);
   else
       c = 1;
       D = distances(G, 'Method', 'unweighted');
   end
   
   % All-pairs shortest path distance
    n = numnodes(G);

    Is = reshape(Is, 1, []);
    Js = reshape(Js, 1, []);

    for i = 1:n
        d = D(i, :);  % distances from node i to others
        for j = 1:n
            if d(j) > kappa
                for ii = 1:length(Is{i})
                    for jj = 1:length(Js{j})
                         cols = Is{i}(ii);
                         rows = Js{j}(jj);
                         K(rows, cols) = 0;
                    end
                end
            end
        end
    end
end
