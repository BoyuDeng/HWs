function x_mapped = map_gll(n, L, m)
    % Generate GLL nodes on the reference interval [-1, 1]
    n = n+1;
    xi = lglnodes(n); % Assuming lglnodes returns exactly n nodes in [-1, 1]
    
    % Initialize mapped nodes vector
    x_mapped = zeros(n * m, 1);
    
    % Loop over each element
    for e = 1:m
        % Define left and right boundaries of the element in [0, L]
        x_left = (e - 1) * (L / m);
        x_right = e * (L / m);
        
        % Map nodes to the interval [x_left, x_right]
        mapped_nodes = (x_right - x_left) / 2 * (xi + 1) + x_left;
        mapped_nodes = flip(mapped_nodes);
        % Store mapped nodes in the correct location within x_mapped
        start_idx = (e - 1) * n + 1;
        end_idx = e * n;
        x_mapped(start_idx:end_idx) = mapped_nodes;
    end
end
