function global_D = assemble_global_matrix(local_D, num_elements)
    % Inputs:
    % local_D: A 2D matrix (size NxN) for one element (local differentiation matrix).
    % num_elements: The number of elements in the domain.
    
    % Get the size of the local differentiation matrix (N+1, typically the number of GLL points)
    Nplus1 = size(local_D, 1);
    N = Nplus1 - 1;  % The number of internal points per element is N
    
    % Total number of global degrees of freedom (DoFs) in the global system
    num_global_points = num_elements * N + 1;
    
    % Initialize the global differentiation matrix with zeros
    global_D = zeros(num_global_points, num_global_points);
    
    % Loop over each element
    for e = 1:num_elements
        % Compute the global degrees of freedom for this element
        % Elements share boundaries, so for element e:
        % - Element 1 maps to global DoFs 1, 2, 3 (N+1 points),
        % - Element 2 maps to DoFs 3, 4, 5, etc.
        global_indices = (e-1)*N + (1:N+1);
        
        % Add the local matrix into the global matrix
        global_D(global_indices, global_indices) = global_D(global_indices, global_indices) + local_D;
    end
end
