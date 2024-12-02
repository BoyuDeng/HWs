function x_mapped = map_gllAD(n, L, m)
    % Generate GLL nodes on the reference interval [-1, 1]
    n = n + 1; % Number of GLL nodes per element
    xi = lglnodes(n); % Assuming lglnodes returns exactly n nodes in [-1, 1]
    
    % Determine the split between 90% and 10% of the domain
    m_90 = floor(2/3 * m); % Number of elements in the 90% region
    m_10 = m - m_90; % Number of elements in the 10% region
    
    % Lengths of the subdomains
    L_90 = 0.9 * L; % 90% of the domain
    L_10 = 0.1 * L; % 10% of the domain
    
    % Initialize mapped nodes vector with appropriate length
    x_mapped = zeros(n * m - (m - 1), 1); % Account for deleting (m - 1) entries
    
    % Loop over each element
    idx = 1;
    for e = 1:m
        % Determine which region the element belongs to
        if e <= m_90
            % Element in the 90% region
            x_left = (e - 1) * (L_90 / m_90);
            x_right = e * (L_90 / m_90);
        else
            % Element in the 10% region
            e_10 = e - m_90; % Local element index in the 10% region
            x_left = L_90 + (e_10 - 1) * (L_10 / m_10);
            x_right = L_90 + e_10 * (L_10 / m_10);
        end
        
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
    x_mapped(:)=x_mapped(:)-L/2;
end
