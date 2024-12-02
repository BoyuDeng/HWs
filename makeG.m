function B = makeG(A, m)
    % makeG - Creates a shifted block diagonal matrix with overlapping
    % A - n*n matrix (input small block)
    % m - number of blocks
    % B - resulting block matrix with A repeated m times, each shifted up by one row and
    %     with the last entry of each block overlapping with the first entry of the next block.

    % Get the size of the input matrix A
    [n, n_col] = size(A);

    % Check if A is square
    if n ~= n_col
        error('Input matrix A must be square');
    end

    % Define A1 and A2
    A1 = A * (2 * m * 2 / 3) / (2 - 0.2);
    A2 = A * (2 * m * 1 / 3) / 0.2;

    % Calculate the size of the full matrix, including zero rows
    total_size = n * m + (m - 1); % Account for row shifts

    % Initialize the matrix with zeros
    B = zeros(total_size, n * m);

    % Loop to place each block in the matrix with a shift
    for i = 1:m
        % Choose A1 or A2 based on the index i relative to 2/3*m
        if i <= floor((2 / 3) * m)
            currentA = A1;
        else
            currentA = A2;
        end

        % Calculate starting positions for rows and columns
        row_start = (i - 1) * n - (i - 1) + 1; % Adjust row start for the upward shift
        col_start = (i - 1) * n - (i - 1) + 1; % Adjust column start for the leftward shift

        % Place the block
        B(row_start:row_start + n - 1, col_start:col_start + n - 1) = ...
            B(row_start:row_start + n - 1, col_start:col_start + n - 1) + currentA;
    end

    % Remove zero columns at the right by finding the last column with data
    last_col = find(any(B, 1), 1, 'last'); % Find the last non-zero column
    B = B(:, 1:last_col); % Trim to include only columns up to the last non-zero column

    % Remove zero rows at the bottom by finding the last row with data
    last_row = find(any(B, 2), 1, 'last'); % Find the last non-zero row
    B = B(1:last_row, :); % Trim to include only rows up to the last non-zero row

    % For fixed boundary conditions
    B(1, :) = 0;
    B(1, 1) = 1;
    B(end, :) = 0;
    B(end, end) = 1;
end
