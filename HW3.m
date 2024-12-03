close all; clear; clc;

% Parameters
C = 1;                  % Advection velocity
kappa = 0.0001;           % Diffusivity
f = 1;                  % Forcing term
L = 2;                  % Domain size (-1 to 1)
N = 12;
M = 24;
t_final = 2.5;
h = 0.0001;





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

t=0;
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


% Steady-state analytical solution
u_steady = (x + 1) - 2 * (exp((x - 1) / kappa) - exp(-2 / kappa)) ./ (1 - exp(-2 / kappa));

% Compute the L_inf norm (maximum absolute difference)
L_inf_norm = max(abs(U(:, end) - u_steady));

% Display the result
fprintf('L_inf norm between U(:,end) and U_steady: %.6f\n', L_inf_norm);


% Plot the numerical and steady-state solutions
figure;
hold on;
plot(x, U(:, end), 'b-', 'LineWidth', 2); % Numerical solution
plot(x, u_steady, 'r--', 'LineWidth', 2); % Steady-state solution
xlabel('x');
ylabel('u(x, t \rightarrow 2.5)');
title('Comparison of Numerical and Steady-State Solutions');
legend('Numerical Solution', 'Steady-State Solution', 'Location', 'Best'); % Add legend
grid on;
hold off;

% Select five equally spaced time steps between 1 and 3000
time_steps_to_plot = floor(linspace(1, 2500, 5));

figure;
hold on;
for i = 1:length(time_steps_to_plot)
    step = time_steps_to_plot(i);
    plot(x, U(:, step), 'LineWidth', 2, 'DisplayName', sprintf('t = %.3f', (step - 1) * h));
end
xlabel('x');
ylabel('u(x, t)');
title('Evolution of Numerical Solution');
legend('Location', 'Best'); % Add legend for all selected timesteps
grid on;
hold off;
















