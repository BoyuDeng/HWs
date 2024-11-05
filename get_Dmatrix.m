function [D,w] = get_Dmatrix(N,L)

[x,w] = lglnodes(N,0,L);



% Compute the differentiation matrix
    D = zeros(N, N);
    c = [2; ones(N-1, 1); 2] .* (-1).^(0:N)';  % c_j values for differentiation matrix

    for i = 1:N
        for j = 1:N
            if i ~= j
                D(i,j) = c(j) / c(i) / (x(i) - x(j));
            end
        end
    end

    % Diagonal elements
    for i = 1:N
        D(i,i) = -sum(D(i, :));
    end
end
