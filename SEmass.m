function B = SEmass(v, m, L)
    % createDiagMatrixFromVector - Creates a large diagonal matrix
    % by repeating a vector along the main diagonal with adjustments
    % v - input vector
    % m - number of times to repeat vector as a diagonal block
    % B - resulting diagonal matrix with m repetitions of v
    
    % Ensure v is a column vector
    v = v(:);
    n = length(v);  % Length of original vector v
    
    
    % Initialize the adjusted vector
    adjusted_v = [];
    
    % Construct adjusted vector with combined entries
    for i = 1:m
        if i == 1
            % For the first block, just add v as it is
            adjusted_v = [adjusted_v; v];
        else
            % For subsequent blocks, combine the last entry of previous block with the first entry of current
            combined_entry = adjusted_v(end) + v(1);
            adjusted_v(end) = combined_entry;  % Replace the last entry with the combined entry
            adjusted_v = [adjusted_v; v(2:end)];  % Append the rest of the current block
        end
    end
    
    % Create a diagonal matrix from the adjusted vector
    B = diag(adjusted_v);
    B = B .* L / (2 * m);
end

