close all;clear;clc;


[x8,w8] = lglnodes(9);

[x16,w16] = lglnodes(17);

[x32,w32] = lglnodes(33);


% Plot the nodes for each case
figure;
hold on;
plot(x8, zeros(size(x8)), 'ro', 'MarkerSize', 8, 'DisplayName', 'N=8');   % Red circles for N=8
plot(x16, zeros(size(x16)) + 0.1, 'go', 'MarkerSize', 8, 'DisplayName', 'N=16');  % Green circles for N=16
plot(x32, zeros(size(x32)) + 0.2, 'bo', 'MarkerSize', 8, 'DisplayName', 'N=32');  % Blue circles for N=32

% Set plot properties
title('Gauss-Lobatto-Legendre Nodes for N=8, N=16, and N=32');
xlabel('x');
ylabel('y (artificial separation)');
legend('show');
ylim([-0.1, 0.3]);  % Add some separation for clarity
grid on;
hold off;

[D,D2,D3] = derv(8,x8);


up1 = (2 * pi) * cos(2 * pi * x8);
ug8 = D*x8;

a = 0;

% Parameters
N = 8;  % Degree of the polynomial (total points = N + 1)
L = 1;  % Length of the domain
b = L;
% Generate GLL nodes (this assumes you have a function to get GLL nodes)
[x, ~] = lglnodes(9); % Gauss-Lobatto-Legendre nodes
x = flipud(x);
% Compute the first, second, and third derivative matrices
[D, D2, D3] = derv(N, x);
[D, D2, D3] = scale_derivative_matrices(D,D2,D3,a,b);
x = (x+1)/2*(b-a)+a;

% Function u(x) = sin(2 * pi * x / L) evaluated at GLL points
u = sin(2 * pi * x / L);

% First derivative of u(x) using the spectral differentiation matrix
u_prime_numerical = D * u;

% Exact derivative for comparison
u_prime_exact = (2 * pi / L) * cos(2 * pi * x / L);

% Display results
disp('Numerical Derivative (u_prime_numerical):');
disp(u_prime_numerical);

disp('Exact Derivative (u_prime_exact):');
disp(u_prime_exact);

% Plot to compare numerical and exact results
figure;
plot(x, u_prime_exact, 'r-', 'LineWidth', 2); hold on;
plot(x, u_prime_numerical, 'bo', 'MarkerSize', 5);
legend('Exact Derivative', 'Numerical Derivative (Spectral)');
title('First Derivative of u(x) = sin(2\pi x / L)');
xlabel('x (GLL points)');
ylabel('Derivative');
grid on;












% Rescale to given interval



%%


% Parameters
alpha = 1.0;    % Diffusivity
tf = 0.25;      % Final time
Nt = 1000;      % Number of time steps
N = 50;        % Number of spatial nodes (LGL nodes)
dt = tf / Nt;   % Time step size
a = 0;
b = 1;
% Initial and Boundary Conditions
initial_condition = @(x) 1 - x - (1 / pi) * sin(2 * pi * x);
time_indices = round([1, Nt/3, 2*Nt/3, Nt+1]);
time_values = time_indices * dt;

% Create the spatial grid using LGL nodes
[x,w] = lglnodes(N+1);
x = (flipud(x)+1)/2*(b-a)+a;
w = w/2*(b-a);

% Initial condition
C = initial_condition(x);
C_at_times = zeros(N+1, 4);  % 4 columns for 4 time points
C_at_times(:, 1) = C;       % Store the initial condition

% Second derivative matrix
[~,D2,~] = derv(N, x);

% Construct the coefficient matrix for the implicit method
A = eye(N+1) - alpha * dt * D2;

% Time-stepping loop (fully implicit scheme)
for n = 1:Nt
    % Apply boundary conditions
    C(1) = 1;
    C(end) = 0;
    
    % Solve linear system directly: A * C_new = C_old
    C = A \ C;
    
    % Store the solution at the selected time steps
    if ismember(n+1, time_indices)
        index = find(time_indices == n+1);
        C_at_times(:, index) = C;
    end
end

% Analytic solution function
analytic_solution = @(x, t) 1 - x - (1 / pi) * exp(-4 * pi^2 * t) * sin(2 * pi * x);

% Compute the analytic solution at the selected time points
C_analytic = zeros(N+1, 4);
for k = 1:4
    C_analytic(:, k) = analytic_solution(x, time_values(k));
end

% Plot the numerical and analytic solutions at the selected time points
figure;
hold on;
colors = lines(4);  % Use different colors for the plots

for i = 1:4
    % Plot numerical solution
    plot(x, C_at_times(:, i), '-o', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Numerical, t = ', num2str(time_values(i))]);

    % Plot analytic solution
    plot(x, C_analytic(:, i), '--', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Analytic, t = ', num2str(time_values(i))]);
end

hold off;
xlabel('x');
ylabel('C(x,t)');
title('Comparison of Numerical and Analytic Solutions at Selected Time Points');
legend show;
grid on;






%%


% Parameters
alpha = 1.0;    % Diffusivity
tf = 0.25;      % Final time
Nt = 1000;      % Number of time steps
N = 10;         % Number of Legendre polynomials (degree)
dt = tf / Nt;   % Time step size
a = 0;
b = 1;
% Create the spatial grid using LGL nodes
[x,w] = lglnodes(N+1);
x = (flipud(x)+1)/2*(b-a)+a;
w = w/2*(b-a);

% Initial condition
C = initial_condition(x);
C_at_times = zeros(N+1, 4);  % 4 columns for 4 time points
C_at_times(:, 1) = C;       % Store the initial condition

[D,~,~] = derv(N, x);

% Construct the coefficient matrix for the implicit method
G = diag(w)^-1*D'*diag(w)*D;

A = dt*G*alpha+eye(N+1);

% Time-stepping loop (fully implicit scheme)
for n = 1:Nt
    % Apply boundary conditions
    C(1) = 1;
    C(end) = 0;
    
    % Solve linear system directly: A * C_new = C_old
    C = A \ C;
    
    % Store the solution at the selected time steps
    if ismember(n+1, time_indices)
        index = find(time_indices == n+1);
        C_at_times(:, index) = C;
    end
end

% Analytic solution function
analytic_solution = @(x, t) 1 - x - (1 / pi) * exp(-4 * pi^2 * t) * sin(2 * pi * x);

% Compute the analytic solution at the selected time points
C_analytic = zeros(Nx+1, 4);
for k = 1:4
    C_analytic(:, k) = analytic_solution(x, time_values(k));
end

% Plot the numerical and analytic solutions at the selected time points
figure;
hold on;
colors = lines(4);  % Use different colors for the plots

for i = 1:4
    % Plot numerical solution
    plot(x, C_at_times(:, i), '-o', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Numerical, t = ', num2str(time_values(i))]);

    % Plot analytic solution
    plot(x, C_analytic(:, i), '--', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Analytic, t = ', num2str(time_values(i))]);
end

hold off;
xlabel('x');
ylabel('C(x,t)');
title('Comparison of Numerical and Analytic Solutions at Selected Time Points');
legend show;
grid on;








%%
% Parameters
alpha = 1.0;    % Diffusivity
tf = 0.25;      % Final time
Nt = 1000;      % Number of time steps
N =8;        % Number of spatial nodes (LGL nodes)
dt = tf / Nt;   % Time step size
a = -1;
b = 1;
m= 2;
npt = N+1+(N)*(m-1);
% Initial and Boundary Conditions
initial_condition = @(x) 1 - x - (1 / pi) * sin(2 * pi * x);
time_indices = round([1, Nt/3, 2*Nt/3, Nt+1]);
time_values = time_indices * dt;

% % Create the spatial grid using LGL nodes
[xlocal,w] = lglnodes(N+1);
xlocal = (flipud(xlocal)+1)/2*(b-a)+a;
w = w/2*(b-a);
x = linspace(0,1,npt);
x = x';
% Initial condition
C = initial_condition(x);
C_at_times = zeros(npt, 4);  % 4 columns for 4 time points
C_at_times(:, 1) = C;       % Store the initial condition

[D,~,~] = derv(N, xlocal);

% Construct the coefficient matrix for the implicit method
G = (2/((1)/m))*diag(w)^-1*D'*diag(w)*D;
BG = assemble_global_matrix(G,m);
% BG(1,:) = 0;
% BG(end,:) = 0;
% BG(1,1) = 1;
% BG(9,:) = BG(9,:)./2;
% BG(end,end) = 1;
A = dt*D+eye(npt);


% Time-stepping loop (fully implicit scheme)
for n = 1:Nt
    % Apply boundary conditions
    C(1) = 1;
    C(end) = 0;
    
    % Solve linear system directly: A * C_new = C_old
    C = A \ C;
    
    % Store the solution at the selected time steps
    if ismember(n+1, time_indices)
        index = find(time_indices == n+1);
        C_at_times(:, index) = C;
    end
end

% Analytic solution function
analytic_solution = @(x, t) 1 - x - (1 / pi) * exp(-4 * pi^2 * t) * sin(2 * pi * x);

% Compute the analytic solution at the selected time points
C_analytic = zeros(npt, 4);
for k = 1:4
    C_analytic(:, k) = analytic_solution(x, time_values(k));
end

% Plot the numerical and analytic solutions at the selected time points
figure;
hold on;
colors = lines(4);  % Use different colors for the plots

for i = 1:4
    % Plot numerical solution
    plot(x, C_at_times(:, i), '-o', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Numerical, t = ', num2str(time_values(i))]);

    % Plot analytic solution
    plot(x, C_analytic(:, i), '--', 'LineWidth', 2, 'Color', colors(i, :), ...
         'DisplayName', ['Analytic, t = ', num2str(time_values(i))]);
end

hold off;
xlabel('x');
ylabel('C(x,t)');
title('Comparison of Numerical and Analytic Solutions at Selected Time Points');
legend show;
grid on;






%%






