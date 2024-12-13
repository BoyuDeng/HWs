function P = legendre_poly(x, degree)
    % legendre_poly: Evaluates the Legendre polynomial of given degree at x.
    % Inputs:
    %   x: Vector of points.
    %   degree: Degree of the Legendre polynomial.
    % Output:
    %   P: Values of the Legendre polynomial at x.

    if degree == 0
        P = ones(size(x)); % P0(x) = 1
    elseif degree == 1
        P = x; % P1(x) = x
    else
        % Recurrence relation for Legendre polynomials
        Pn_2 = ones(size(x)); % P0(x)
        Pn_1 = x;             % P1(x)
        for n = 2:degree
            P = ((2*n - 1) .* x .* Pn_1 - (n - 1) .* Pn_2) / n;
            Pn_2 = Pn_1;
            Pn_1 = P;
        end
    end
end