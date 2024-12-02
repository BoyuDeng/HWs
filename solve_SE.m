
function U_final = solve_SE(M, N, Ini)
    % Solve the advection problem using the Spectral Element (S.E.) method with BDF/AB3
    % Inputs:
    % M   - Number of elements
    % N   - Number of nodes per element
    % Ini - Initial condition selection (1 for Gaussian, 2 for Cone)
    
    % Parameters
    L = 1000;      
    lambda0 = 30;  
    s = 1/25;      
    x0 = 250;   
    t = 0;         
    t_final = 2000; 
    h = 0.05;       
    
    % Generate the nodes and weights for each element
    [xi, w] = lglnodes(N + 1); 
    xi = flipud(xi);         

    % Differentiation matrix and stiffness/mass matrices
    D = derv(N, xi);
    Ase = StiffSE(D, M, w);
    SMass = SEmass(w, M, L);
    xse = map_gllse(N, L, M);  % Map nodes to the spatial domain

    % Set initial condition based on Ini value
    if Ini == 1
        % Gaussian initial condition
        u0 = exp(-0.5 * ((xse - x0).^2) / lambda0^2);
    elseif Ini == 2
        % Cone initial condition
        u0 = max(0, 1 - s * abs(xse - x0));
    else
        error('Invalid value for Ini. Use 1 for Gaussian or 2 for Cone.');
    end

    % Pre-allocate solution matrix for time-stepping
    num_steps = floor(t_final / h);   % Total number of time steps
    U = zeros(length(xse), num_steps + 1);  % Solution matrix
    U(:, 1) = u0;  % Set initial condition as first column in U

    % Precompute the matrix product for efficiency
    H = inv(SMass) * Ase;

    % BDF3 and AB3 coefficients
    a0 = 11/6; a1 = -3; a2 = 3/2; a3 = -1/3;
    b = [3, -3, 1];

    % Set initial conditions for the first few time steps
    U(:, 2) = u0; 
    U(:, 3) = u0;  

    % Time-stepping loop
    n = 3;  % Start at the third time step
    while t < t_final
        RHS = 0;
        for j = 1:3
            RHS = RHS + b(j) * H * U(:, n - j + 1);
        end
 
        U(:, n + 1) = (h * RHS - a1 * U(:, n) - a2 * U(:, n - 1) - a3 * U(:, n - 2)) / a0;


        t = t + h;
        n = n + 1;
    end

    U_final = U(:, end);

end
