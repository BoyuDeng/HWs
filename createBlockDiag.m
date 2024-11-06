function B = createBlockDiag(A, m, w, L)
    % createBlockDiag - Creates a block diagonal matrix
    % A - n*n matrix (input small block)
    % m - number of blocks
    % B - resulting block diagonal matrix with A repeated m times
    
    % Get the size of the input matrix A
    [n, n_col] = size(A);
    
    % Check if A is square
    if n ~= n_col
        error('Input matrix A must be square');
    end
    
    % Initialize an empty block diagonal matrix
    B = [];
    A = (2*m/L).*A;
    A = A.*w;
    % Loop to create the block diagonal matrix
    for i = 1:m
        B = blkdiag(B, A);
        B(n*i,n*i) = B(n*i,n*i) + 1;
        
    end
    for i = 1:m-1
        B(n*i+1, n*i) = B(n*i+1, n*i)-1;
    end
    B(1,end) = B(1,end) -1;
end
