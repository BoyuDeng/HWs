function [D, D2, D3] = derv(nterm, x)
    % Input:
    %   nterm - The order of polynomial approximation N (total points = N + 1)
    %   x - Vector of GLL grid points on the interval [-1, 1]
    %
    % Output:
    %   D - First derivative matrix
    %   D2 - Second derivative matrix
    %   D3 - Third derivative matrix
    % Initialize matrices
    D = zeros(nterm+1, nterm+1);
    D2 = zeros(nterm+1, nterm+1);
    D3 = zeros(nterm+1, nterm+1);
    
    % Compute coefficients C_i
    c = zeros(1, nterm+1);
    for k = 0:nterm
        prod = 1;
        xk = x(k+1); % MATLAB is 1-based indexing
        for l = 0:nterm
            xl = x(l+1);
            if l ~= k
                prod = prod * (xk - xl);
            end
        end
        c(k+1) = prod;
    end

    % FIRST DERIVATIVE MATRIX (D)
    % Off-diagonal elements
    for k = 0:nterm
        xk = x(k+1);
        for j = 0:nterm
            xj = x(j+1);
            if k ~= j
                D(k+1,j+1) = c(k+1) / (c(j+1) * (xk - xj));
            else
                D(k+1,j+1) = 0;
            end
        end
    end
    
    % Diagonal elements
    for k = 0:nterm
        sum = 0;
        for j = 0:nterm
            if k ~= j
                sum = sum + D(k+1,j+1);
            end
        end
        D(k+1,k+1) = -sum;
    end

    % SECOND DERIVATIVE MATRIX (D2)
    m = 2;
    for k = 0:nterm
        xk = x(k+1);
        for j = 0:nterm
            xj = x(j+1);
            if k ~= j
                D2(k+1,j+1) = m * (D(k+1,k+1) * D(k+1,j+1) - D(k+1,j+1) / (xk - xj));
            else
                D2(k+1,j+1) = 0;
            end
        end
    end
    
    % Diagonal elements for D2
    for k = 0:nterm
        sum = 0;
        for j = 0:nterm
            if k ~= j
                sum = sum + D2(k+1,j+1);
            end
        end
        D2(k+1,k+1) = -sum;
    end

    % THIRD DERIVATIVE MATRIX (D3)
    m = 3;
    for k = 0:nterm
        xk = x(k+1);
        for j = 0:nterm
            xj = x(j+1);
            if k ~= j
                D3(k+1,j+1) = m * (D2(k+1,k+1) * D(k+1,j+1) - D2(k+1,j+1) / (xk - xj));
            else
                D3(k+1,j+1) = 0;
            end
        end
    end
    
    % Diagonal elements for D3
    for k = 0:nterm
        sum = 0;
        for j = 0:nterm
            if k ~= j
                sum = sum + D3(k+1,j+1);
            end
        end
        D3(k+1,k+1) = -sum;
    end
end
