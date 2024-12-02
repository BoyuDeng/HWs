function B = makeGsym(A, m)
    % createShiftedBlockDiag - Creates a shifted block diagonal matrix with overlapping
    % A - n*n matrix (input small block)
    % m - number of blocks
    % w - weight to multiply A by
    % B - resulting block matrix with A repeated m times, each shifted up by one row and
    %     with the last entry of each block overlapping with the first entry of the next block.
    
    % Get the size of the input matrix A
    [n, n_col] = size(A);

    % Check if A is square
    if n ~= n_col
        error('Input matrix A must be square');
    end
    
    % Calculate the size of the full matrix, including zero rows
    total_size = n * m + (m - 1);  % Account for row shifts
    
    % Initialize the matrix with zeros
    B = zeros(total_size, n * m);
    
    % Loop to place each block in the matrix with a shift
    for i = 1:m
        row_start = (i - 1) * n - (i - 1) + 1;  % Adjust row start for the upward shift
        col_start = (i - 1) * n - (i - 1) + 1;  % Adjust column start for the leftward shift
        
        % Place the block A
        B(row_start:row_start + n - 1, col_start:col_start + n - 1) = ...
            B(row_start:row_start + n - 1, col_start:col_start + n - 1) + A;
        
        % If this is not the last block, combine the overlapping entry
        % of this block with the next one
        % if i < m
        %     B(row_start + n - 1, col_start + n) = ...
        %         B(row_start + n - 1, col_start + n) + B(row_start + n - 1, col_start + n - 1);
        % end
    end
    
    % Remove zero columns at the right by finding the last column with data
    last_col = find(any(B, 1), 1, 'last');  % Find the last non-zero column
    B = B(:, 1:last_col);  % Trim to include only columns up to the last non-zero column
    
    % Remove zero rows at the bottom by finding the last row with data
    last_row = find(any(B, 2), 1, 'last');  % Find the last non-zero row
    B = B(1:last_row, :);  % Trim to include only rows up to the last non-zero row


    %for periodic boundary condition
    % B(1,1) = B(1,1)+B(end,end);
    % B(end,end) = B(1,1);
    % 
    % % Copy the first half of the first row to the last row (same columns)
    % half_point = floor(size(B, 2) / 2);
    % B(end, 1:half_point) = B(1, 1:half_point);
    % 
    % % Copy the second half of the last row to the first row (same columns)
    % B(1, half_point+1:end) = B(end, half_point+1:end);


    %for fix boundary condition
    B(1,:) = 0;
    B(1,1) = 1;
    B(end,:) = 0;
    B(end,end) = 1;

end 