function V = legendre_vandermonde(x, n)
    % legendre_vandermonde: Constructs a Legendre Vandermonde matrix.
    % Inputs:
    %   x: Vector of points at which the Legendre polynomials are evaluated.
    %   n: Number of columns (degree of the polynomial + 1).
    % Output:
    %   V: Legendre Vandermonde matrix.

    % Number of points
    m = length(x);
    
    % Initialize the Vandermonde matrix
    V = zeros(m, n);
    
    % Fill the Vandermonde matrix
    for k = 1:n
        % Evaluate the (k-1)-th Legendre polynomial at points x
        V(:, k) = legendre_poly(x, k-1);
    end
end