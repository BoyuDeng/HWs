function x_mapped = map_gllse(n, L, m)
    % Generate GLL nodes on the reference interval [-1, 1]
    n = n + 1;
    xi = lglnodes(n); % Assuming lglnodes returns exactly n nodes in [-1, 1]
    
    % Initialize mapped nodes vector with appropriate length
    x_mapped = zeros(n * m - (m - 1), 1); % Account for deleting (m - 1) entries
    
    % Loop over each element
    idx = 1;
    for e = 1:m
        % Define left and right boundaries of the element in [0, L]
        x_left = (e - 1) * (L / m);
        x_right = e * (L / m);
        
        % Map nodes to the interval [x_left, x_right]
        mapped_nodes = (x_right - x_left) / 2 * (xi + 1) + x_left;
        mapped_nodes = flip(mapped_nodes);
        
        % Store mapped nodes, skipping the first entry of each element after the first
        if e == 1
            x_mapped(idx:idx + n - 1) = mapped_nodes;
            idx = idx + n;
        else
            x_mapped(idx:idx + n - 2) = mapped_nodes(2:end);
            idx = idx + n - 1;
        end
    end
end
