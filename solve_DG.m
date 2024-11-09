function U_final = solve_DG(M,N, Ini)

% Parameters
L = 1000;      % Domain length (meters)
lambda0 = 30;  % Lengthscale (meters)
s = 1/25;      % Scaling factor
C = 1;         % Velocity (m/s)
x0 = 250;      % Initial center (meters)
t_final = 2000;              % Final time for simulation
t = 0;             % Initial time
h = 0.05;           % Time step size


% Parameters        % Length of the domain
M = M;         % Number of elements
N = N;         % Number of nodes per element




% Generate the nodes and weights for each element
[xi, w] = lglnodes(N+1); % Gauss-Legendre nodes and weights on [-1, 1]

% Map nodes from [-1, 1] to each element in [0, L]
xi = flipud(xi);


% Differentiation matrix for Legendre-Gauss-Lobatto points
D = derv(N,xi);

Adg = StiffDG(D, M, w);

DMass = DGmass(w,M,L);

x = map_gll(N,L,M);

    % Set initial condition based on Ini value
    if Ini == 1
        % Gaussian initial condition
        u0 = exp(-0.5 * ((x - x0).^2) / lambda0^2);
    elseif Ini == 2
        % Cone initial condition
        u0 = max(0, 1 - s * abs(x - x0));
    else
        error('Invalid value for Ini. Use 1 for Gaussian or 2 for Cone.');
    end

% Define initial conditions 
U = u0;   % Initial state vector



% Time-stepping loop with RK4 and animation
while t < t_final
    % Calculate k1
    k1 = h * (DMass \ (Adg * U));

    % Calculate k2
    k2 = h * (DMass \ (Adg * (U + 0.5 * k1)));

    % Calculate k3
    k3 = h * (DMass \ (Adg * (U + 0.5 * k2)));

    % Calculate k4
    k4 = h * (DMass \ (Adg * (U + k3)));

    % Update U
    U = U + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4);              
         

    % Update time
    t = t + h;
end

U_final = U;