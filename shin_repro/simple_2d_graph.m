function G = simple_2d_graph(n)
%SIMPLE_2D_GRAPH Build an n×n 2D grid graph with edge weight 0.05.
% Returns a graph with n^2 nodes. Node (i,j) has index n*(i-1)+j.

    N = n * n;
    A = zeros(N);
    weight = 0.05;

    for i = 1:n
        for j = 1:n
            node = n*(i-1) + j;
            % horizontal edge to (i, j+1)
            if j < n
                neighbor = n*(i-1) + j + 1;
                A(node, neighbor) = weight;
                A(neighbor, node) = weight;
            end
            % vertical edge to (i+1, j)
            if i < n
                neighbor = n*i + j;
                A(node, neighbor) = weight;
                A(neighbor, node) = weight;
            end
        end
    end

    G = graph(A);
end
