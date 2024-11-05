function [D_scaled, D2_scaled, D3_scaled] = scale_derivative_matrices(D, D2, D3, a, b)
    % Input:
    %   D - First derivative matrix on [-1, 1]
    %   D2 - Second derivative matrix on [-1, 1]
    %   D3 - Third derivative matrix on [-1, 1]
    %   a, b - The target domain [a, b]
    %
    % Output:
    %   D_scaled - First derivative matrix on [a, b]
    %   D2_scaled - Second derivative matrix on [a, b]
    %   D3_scaled - Third derivative matrix on [a, b]
    
    % Compute the scaling factor for the new domain
    scaling_factor = 2 / (b - a);

    % Scale the derivative matrices
    D_scaled = D * scaling_factor;           % First derivative
    D2_scaled = D2 * scaling_factor;       % Second derivative
    D3_scaled = D3 * scaling_factor;       % Third derivative
end
