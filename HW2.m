close all; clear; clc;

% Parameters
L = 1000;      % Domain length (meters)
lambda0 = 30;  % Lengthscale (meters)
s = 1/25;      % Scaling factor
C = 1;         % Velocity (m/s)
x0 = 250;      % Initial center (meters)
T = L / C;     % Time for one cycle




% Parameters        % Length of the domain
M = 10;         % Number of elements
N = 16;         % Number of nodes per element
C = 1;         % Constant advection velocity



% Generate the nodes and weights for each element
[xi, w] = lglnodes(N+1); % Gauss-Legendre nodes and weights on [-1, 1]

% Map nodes from [-1, 1] to each element in [0, L]
xi = flipud(xi);
%w = w*(L/M)/2;
% Initial Condition (Gaussian bump)


% Differentiation matrix for Legendre-Gauss-Lobatto points
D = derv(N,xi);

Adg = StiffDG(D, M, w);
Ase = StiffSE(D,M, w);
DMass = DGmass(w,M,L);
SMass = SEmass(w,M,L);


x = map_gll(N,L,M);
xse = map_gllse(N,L,M);

%%
% Initial Condition: Gaussian bump
u0_gaussian = exp(-0.5 * ((x - x0).^2) / lambda0^2);

% Initial Condition: Cone
u0_cone = max(0, 1 - s * abs(x - x0));

% Define initial conditions 
U = u0_gaussian;   % Initial state vector
t = 0;             % Initial time
t_final = 1200;      % Final time for simulation
h = 0.1;           % Time step size


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

hold on
plot(x, u0_gaussian);
plot(x, U);
hold off

%%
% Define the spatial grid and initial Gaussian condition
u0_gaussian = exp(-0.5 * ((xse - x0).^2) / lambda0^2);

% Initialize other parameters
t = 0;                       % Initial time
t_final = 1200;              % Final time for simulation
h = 0.1;                    % Time step size
num_steps = floor(t_final / h);  % Total number of time steps

% Pre-allocate U for all time steps
U = zeros(length(xse), num_steps + 1);  % Pre-allocate with size for efficiency
U(:,1) = u0_gaussian;  % Initial state vector at t=0

% Precompute matrices
I = size(Ase);
H = inv(SMass)*Ase;

% BDF3 coefficients (for implicit step)
a0 = 11/6;
a1 = -3;
a2 = 3/2;
a3 = -1/3;

% AB3 coefficients (for explicit step)
b(1) = 3;
b(2) = -3;
b(3)= 1;

% Set initial conditions for the first few time steps
U(:,2) = u0_gaussian;  % Placeholder for u^{n-1}
U(:,3) = u0_gaussian;  % Placeholder for u^{n-2}

% Time-stepping loop
n = 3;  % Start at third time step (for AB3/BDF3)
while t < t_final
RHS = 0;
    for j = 1:3
        RHS = RHS + b(j)*H*U(:,n-j+1);
    end

    U(:,n+1) = (h*RHS - a1*U(:,n) - a2*U(:,n-1) - a3*U(:,n-2))/a0;
    
    % Increment time and shift indices
    t = t + h;
    n = n + 1;
end

% Plot initial and final solutions
hold on
plot(xse, u0_gaussian, 'DisplayName', 'Initial Condition');
plot(xse, U(:,end), 'DisplayName', 'Final Solution');
hold off
legend;
title('BDF/AB3 Time Stepping Solution');
xlabel('x');
ylabel('u(x,t)');







%%

% Define the spatial grid and initial Gaussian condition
u0_gaussian = exp(-0.5 * ((xse - x0).^2) / lambda0^2);

% Initialize other parameters
U = u0_gaussian;              % Initial state vector
t = 0;                        % Initial time
t_final = 1200;                % Final time for simulation
h = 0.01;                     % Time step size

% Precompute the implicit update matrix

% For Backward Euler: (I - h * SMass \ Ase)
%implicitMatrix = (SMass - h * Ase);


% Time-stepping loop with animation
I = size(Ase);
while t < t_final
    % Solve for U_new using the implicit time-stepping formula
    % (SMass - h * Ase) * U_new = SMass * U
    %U = implicitMatrix \ (SMass * U);
    U = (h*inv(SMass)*Ase+eye(I(1)))*U;
    t = t + h;
end

hold on
plot(xse, u0_gaussian);
plot(xse, U);
hold off



%%
% Define the grid for the complex plane
real_z = linspace(-3, 3, 500); % Real part range
imag_z = linspace(-3, 3, 500); % Imaginary part range
[Re, Im] = meshgrid(real_z, imag_z);
z = Re + 1i * Im;

% Calculate the amplification factor for RK4
R = 1 + z + (z.^2) / factorial(2) + (z.^3) / factorial(3) + (z.^4) / factorial(4);

% Calculate the magnitude of R
R_magnitude = abs(R);

% Plot the stability region (where |R| <= 1)
% contourf(Re, Im, R_magnitude, [0, 1], 'LineColor', 'none');
% hold on;
% colormap([1 1 1; 0.6 0.8 1]) % White outside the stability region, blue inside
% colorbar;
% title('Stability Region of RK4 Method');
% xlabel('Real part of \lambda \Delta t');
% ylabel('Imaginary part of \lambda \Delta t');
% axis equal;
% hold off;
I = size(Adg);
adMatrix = inv(DMass)*Adg;
% Calculate eigenvalues of the D.G. advection operator
eigenvalues = eig(adMatrix);  % where A_DG is your discretized advection matrix

% Scale eigenvalues by the timestep
dt = 0.3;
ze = eigenvalues * dt;

% Plot stability region and scaled eigenvalues
figure;
hold on;
% Stability region plot
contourf(Re, Im, R_magnitude, [0, 1], 'LineColor', 'none');
colormap([1 1 1; 0.6 0.8 1])
colorbar;


% Eigenvalues plot
plot(real(ze), imag(ze), 'ro', 'MarkerSize', 5, 'LineWidth', 1.5);
title('RK4 Stability Region with Eigenvalues of Advection Operator, M =20, N =32');
xlabel('Real part of \lambda \Delta t');
ylabel('Imaginary part of \lambda \Delta t');
axis equal;
hold off;



































