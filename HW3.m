close all; clear; clc;

% Parameters
C = 1;                  % Advection velocity
kappa = 0.01;           % Diffusivity
f = 1;                  % Forcing term
L = 2;                  % Domain size (-1 to 1)
N = 12;
M = 12;
t_final = 3;
h = 0.001;


% Generate the nodes and weights for each element
    [xi, w] = lglnodes(N + 1); 
    xi = flipud(xi);         
    t=0;
    % Differentiation matrix and stiffness/mass matrices
    D = derv(N, xi);
    Ase = StiffSE(D, M, w);
    SMass = SEmass(w, M, L);
    x = map_gllse(N, L, M)-1;  % Map nodes to the spatial domain

    % Set initial condition based on Ini value
    u0 = sin(2 * pi * (x + 1));

    % Pre-allocate solution matrix for time-stepping
    num_steps = floor(t_final / h);   % Total number of time steps
    U = zeros(length(x), num_steps + 1);  % Solution matrix
    U(:, 1) = u0;  % Set initial condition as first column in U

    % Precompute the matrix product for efficiency
    H = inv(SMass) * Ase;
    H(1,1) = 1;
    H(end,end) = 1;
    G = D'*diag(w)*D.* (2 * M)/L;
    Gglobalsym = makeGsym(G,M);
    Aad = inv(SMass) * Gglobalsym;

    A = h*Aad*kappa+eye(length(u0));
    A(1,1) = 1;
    A(end,end) = A(1,1);

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
        RHS = RHS + f * ones(size(U(:, 1)));
 
        Ustar = (h * RHS - a1 * U(:, n) - a2 * U(:, n - 1) - a3 * U(:, n - 2)) / a0;

        U(:, n + 1) = A\(Ustar);
        %U(2:(end-1), n + 1) = U(2:(end-1), n + 1) + h*f; %add forcing term with boundary condtion
        U(1, n+1) = 0;
        U(end, n+1) = 0;


        t = t + h;
        n = n + 1;

    end

    U_final = U(:, end);
% Compute the steady-state solution
u_steady = (x + 1) - 2 * (exp((x - 1) / kappa) - exp(-2 / kappa)) ./ (1 - exp(-2 / kappa));

% Plot the numerical and steady-state solutions
figure;
hold on;
plot(x, U(:, end), 'b-', 'LineWidth', 2); % Numerical solution
plot(x, u_steady, 'r--', 'LineWidth', 2); % Steady-state solution
xlabel('x');
ylabel('u(x, t \rightarrow \infty)');
title('Comparison of Numerical and Steady-State Solutions');
legend('Numerical Solution', 'Steady-State Solution', 'Location', 'Best'); % Add legend
grid on;
hold off;

mag=sum(U(:,end));


%%


%Generate the nodes and weights for each element
[xi, w] = lglnodes(N + 1); 
xi = flipud(xi);         

%Differentiation matrix and stiffness/mass matrices
D = derv(N, xi);
Ase = StiffSE(D, M, w);
SMass = SEmassAD(w, M, L);
x = map_gllAD(N, L, M);  % Map nodes to the spatial domain





% Initial condition
u0 = sin(2 * pi * (x + 1));

h=0.001;
t=0;
t_final= 5;
    num_steps = floor(t_final / h);   % Total number of time steps
    U = zeros(length(x), num_steps + 1);  % Solution matrix
    U(:, 1) = u0;  % Set initial condition as first column in U

    % Precompute the matrix product for efficiency
    B = inv(SMass) * Ase;
    B(1,1) = 1;
    B(end,end) = B(1,1);


    G = D'*diag(w)*D;
    Gglobal = makeG(G,M);
    Aad = inv(SMass) * Gglobal;

    A = h*Aad*kappa+eye(length(u0));
    A(1,1) = 1;
    A(end,end) = A(1,1);

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
            RHS = RHS + b(j) * B * U(:, n - j + 1);
        end
        RHS = RHS + f * ones(size(U(:, 1)));
 
        Ustar = (h * RHS - a1 * U(:, n) - a2 * U(:, n - 1) - a3 * U(:, n - 2)) / a0;
        U(:, n + 1) = A\(Ustar);
        U(1, n+1) = 0;
        U(end, n+1) = 0;
        t = t + h;
        n = n + 1;

    end

    U_final = U(:, end);




mag=sum(abs(U));
mag=mag-mag(1);
mag=sum(abs(mag));









% Steady-state analytical solution
u_steady = (x + 1) - 2 * (exp((x - 1) / kappa) - exp(-2 / kappa)) ./ (1 - exp(-2 / kappa));

% Plot the numerical and steady-state solutions
figure;
hold on;
plot(x, U(:, end), 'b-', 'LineWidth', 2); % Numerical solution
plot(x, u_steady, 'r--', 'LineWidth', 2); % Steady-state solution
xlabel('x');
ylabel('u(x, t \rightarrow \infty)');
title('Comparison of Numerical and Steady-State Solutions');
legend('Numerical Solution', 'Steady-State Solution', 'Location', 'Best'); % Add legend
grid on;
hold off;

%%

real_z = linspace(-3, 3, 500);
imag_z = linspace(-3, 3, 500); 
[Re, Im] = meshgrid(real_z, imag_z);
z = Re + 1i * Im;

% Calculate the amplification factor for RK4
R = 1 + z + (z.^2) / factorial(2) + (z.^3) / factorial(3) + (z.^4) / factorial(4);

% Calculate the magnitude of R
R_magnitude = abs(R);


I = size(H);
adMatrix = H;
eigenvalues = eig(adMatrix);  

% Scale eigenvalues by the timestep
dt = 0.005;
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












