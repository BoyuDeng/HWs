function B = SEmassAD(v, m, L)
    % createDiagMatrixFromVector - Creates a large diagonal matrix
    % by repeating a vector along the main diagonal with adjustments
    % v - input vector
    % m - number of times to repeat vector as a diagonal block
    % B - resulting diagonal matrix with m repetitions of v
    
    % Ensure v is a column vector
    v = v(:);
    n = length(v);  % Length of original vector v
    
        v1 = v * (L-0.2) / (2 * m *2/3);
        v2 = v * (0.2) / (2 * m *1/3);
% Initialize adjusted vector
adjusted_v = [];

for i = 1:m
    if i == 1
        % For the first block, just add v as it is
        adjusted_v = [adjusted_v; v1];
    else
        % Determine the vector to use (v1 or v2)
        if i <= m * 2 / 3
            vec_to_append = v1;  % Use v1
        else
            vec_to_append = v2;  % Use v2
        end

        % Combine the last entry of the previous block with the first entry of the current vector
        combined_entry = adjusted_v(end) + vec_to_append(1);
        adjusted_v(end) = combined_entry;  % Replace the last entry with the combined entry
        adjusted_v = [adjusted_v; vec_to_append(2:end)];  % Append the rest of the current vector
    end
end

    
    % Create a diagonal matrix from the adjusted vector
    B = diag(adjusted_v);

end

