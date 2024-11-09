function L2_error = compute_L2(u_num, u_ref, x)
    % Compute the L2 norm of the error between numerical and reference solutions
    % while handling duplicate x locations and variable grid spacing.
    % Inputs:
    % u_num - Numerical solution vector
    % u_ref - Reference or exact solution vector
    % x     - Spatial grid vector
    %
    % Output:
    % L2_error - L2 norm of the error

    % Ensure inputs are column vectors
    if isrow(u_num), u_num = u_num'; end
    if isrow(u_ref), u_ref = u_ref'; end
    if isrow(x), x = x'; end
    
    % Remove duplicate x values by averaging corresponding u_num and u_ref
    [x_unique, ~, idx_unique] = unique(x);
    u_num_unique = accumarray(idx_unique, u_num, [], @mean);
    u_ref_unique = accumarray(idx_unique, u_ref, [], @mean);

    % Calculate the spacing between each unique x point
    dx = diff(x_unique);

    % Calculate squared error between numerical and reference solutions at each point
    error_squared = (u_num_unique(1:end-1) - u_ref_unique(1:end-1)).^2;

    % Compute L2 norm using variable spacing
    L2_error = sqrt(sum(error_squared .* dx));
end

