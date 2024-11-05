function B = DGmass(v, m, L)
    % createDiagMatrixFromVector - Creates a large diagonal matrix
    % by repeating a vector along the main diagonal
    % v - input vector
    % m - number of times to repeat vector as a diagonal block
    % B - resulting diagonal matrix with m repetitions of v

    % Ensure v is a column vector
    v = v(:);
    
    % Repeat the vector m times by creating a larger vector
    repeated_v = repmat(v, m, 1);
    
    % Create a diagonal matrix from the repeated vector
    B = diag(repeated_v);
    B = B.*L/(2*m);
end
