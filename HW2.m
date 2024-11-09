close all; clear; clc;

% Parameters
L = 1000;      % Domain length (meters)
lambda0 = 30;  % Lengthscale (meters)
s = 1/25;      % Scaling factor
C = 1;         % Velocity (m/s)
x0 = 250;      % Initial center (meters)



[xi, w] = lglnodes(N+1); % Gauss-Legendre nodes and weights on [-1, 1]

% Map nodes from [-1, 1] to each element in [0, L]
xi = flipud(xi);


D = derv(N,xi);

Adg = StiffDG(D, M, w);
Ase = StiffSE(D,M, w);
DMass = DGmass(w,M,L);
SMass = SEmass(w,M,L);


x = map_gll(N,L,M);
xse = map_gllse(N,L,M);



%%
%%time step
% Define the grid for the complex plane
real_z = linspace(-3, 3, 500);
imag_z = linspace(-3, 3, 500); 
[Re, Im] = meshgrid(real_z, imag_z);
z = Re + 1i * Im;

% Calculate the amplification factor for RK4
R = 1 + z + (z.^2) / factorial(2) + (z.^3) / factorial(3) + (z.^4) / factorial(4);

% Calculate the magnitude of R
R_magnitude = abs(R);


I = size(Adg);
adMatrix = inv(DMass)*Adg;
eigenvalues = eig(adMatrix);  

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

%%
% Define parameters
Ms = [5, 10, 20];  
Ns = [4, 8, 16, 32];  
L2_errors_SE = zeros(length(Ms), length(Ns));
L2_errors_DG = zeros(length(Ms), length(Ns));


% Loop over each M and N to compute solutions and errors
for i = 1:length(Ms)
    M = Ms(i);
    
    for j = 1:length(Ns)
        N = Ns(j);
        
        % Compute solutions for S.E. and D.G. based on M and N
        u_SE = solve_SE(M,N,1);
        u_DG = solve_DG(M,N,1);
        
        % Assuming HW2 has a function to compute the reference solution
        xse = map_gllse(N,L,M);
        xdg = map_gll(N,L,M);
        u_abs_SE = exp(-0.5 * ((xse - x0).^2) / lambda0^2);
        u_abs_DG = exp(-0.5 * ((xdg - x0).^2) / lambda0^2);
        
        % Calculate L2 norm of the error
        L2_errors_SE(i, j) = compute_L2(u_SE, u_abs_SE, xse);
        L2_errors_DG(i, j) = compute_L2(u_DG, u_abs_DG, xdg);
    end
end

% Plotting the L2 norm of the error as a function of N for each M
for i = 1:length(Ms)
    M = Ms(i);
    figure;
    plot(Ns, L2_errors_SE(i, :), '-o', 'DisplayName', 'Spectral Element (S.E.)');
    hold on;
    plot(Ns, L2_errors_DG(i, :), '-s', 'DisplayName', 'Discontinuous Galerkin (D.G.)');
    xlabel('Polynomial Degree N');
    ylabel('L2 Norm of Error');
    title(['L2 Error Norm for M = ' num2str(M) ' Elements']);
    legend;
    hold off;
end

%%
M=10;
N=16;
% Generate spatial grids for SE and DG methods
x = map_gllse(N, L, M);    % Spectral Element (S.E.) grid
xdg = map_gll(N, L, M);     % Discontinuous Galerkin (D.G.) grid

% Compute solutions
Udg = solve_DG(M, N, 1);    % D.G. solution
Use = solve_SE(M, N, 1);    % S.E. solution
u_abs = exp(-0.5 * ((x - x0).^2) / lambda0^2); % Exact solution (Gaussian)

% Handle duplicates in xdg by averaging Udg values at each unique x
[x_unique, ~, idx_unique] = unique(xdg);
Udg = accumarray(idx_unique, Udg, [], @mean);

% Plot the results
figure;
hold on;
plot(x_unique, Udg, '-s', 'LineWidth', 1.5, 'DisplayName', 'Discontinuous Galerkin (D.G.)');
plot(x, Use, '-o', 'LineWidth', 1.5, 'DisplayName', 'Spectral Element (S.E.)');
plot(x, u_abs, '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2], 'DisplayName', 'Exact Solution');
hold off;

% Add plot enhancements
title('Comparison of Discontinuous Galerkin (D.G.) and Spectral Element (S.E.) Methods');
xlabel('Spatial Position (x)');
ylabel('Solution Value');
legend('Location', 'Best');
grid on;


%%

% Pre-allocate arrays to store differences
mass_diff_SE = zeros(length(Ms), length(Ns));
mass_diff_DG = zeros(length(Ms), length(Ns));
energy_diff_SE = zeros(length(Ms), length(Ns));
energy_diff_DG = zeros(length(Ms), length(Ns));

% Loop over each value of M and N
for i = 1:length(Ms)
    M = Ms(i);
    
    for j = 1:length(Ns)
        N = Ns(j);

        % S.E. Method
        x_se = map_gllse(N, L, M);
        u0_se = exp(-0.5 * ((x_se - x0).^2) / lambda0^2);  % Initial condition
        U_final_SE = solve_SE(M, N, 1);                    % Final solution using S.E.
        
        % Compute initial and final mass and energy for S.E.
        mass_initial_SE = trapz(x_se, u0_se);
        mass_final_SE = trapz(x_se, U_final_SE);
        energy_initial_SE = trapz(x_se, u0_se.^2);
        energy_final_SE = trapz(x_se, U_final_SE.^2);
        
        % Store mass and energy differences for S.E.
        mass_diff_SE(i, j) = abs(mass_final_SE - mass_initial_SE);
        energy_diff_SE(i, j) = abs(energy_final_SE - energy_initial_SE);

        % D.G. Method
        x_dg = map_gll(N, L, M);
        u0_dg = exp(-0.5 * ((x_dg - x0).^2) / lambda0^2);   % Initial condition
        U_final_DG = solve_DG(M, N, 1);                    % Final solution using D.G.
        
        % Handle duplicate points in D.G. grid (if applicable)
        [x_dg_unique, ~, idx_dg] = unique(x_dg);
        u0_dg = accumarray(idx_dg, u0_dg, [], @mean);
        U_final_DG = accumarray(idx_dg, U_final_DG, [], @mean);

        % Compute initial and final mass and energy for D.G.
        mass_initial_DG = trapz(x_dg_unique, u0_dg);
        mass_final_DG = trapz(x_dg_unique, U_final_DG);
        energy_initial_DG = trapz(x_dg_unique, u0_dg.^2);
        energy_final_DG = trapz(x_dg_unique, U_final_DG.^2);
        
        % Store mass and energy differences for D.G.
        mass_diff_DG(i, j) = abs(mass_final_DG - mass_initial_DG);
        energy_diff_DG(i, j) = abs(energy_final_DG - energy_initial_DG);
    end
end

% Plot results for each value of M
for i = 1:length(Ms)
    M = Ms(i);

    % Create a new figure for each value of M
    figure;

    % Plot mass difference for S.E. and D.G. as a function of N
    subplot(2, 1, 1);
    hold on;
    plot(Ns, mass_diff_SE(i, :), '-o', 'LineWidth', 1.5, 'DisplayName', ['S.E., M = ' num2str(M)]);
    plot(Ns, mass_diff_DG(i, :), '-s', 'LineWidth', 1.5, 'DisplayName', ['D.G., M = ' num2str(M)]);
    hold off;
    title(['Mass Conservation Error for M = ', num2str(M)]);
    xlabel('Polynomial Order (N)');
    ylabel('Mass Difference (|final - initial|)');
    legend('Location', 'Best');
    grid on;

    % Plot energy difference for S.E. and D.G. as a function of N
    subplot(2, 1, 2);
    hold on;
    plot(Ns, energy_diff_SE(i, :), '-o', 'LineWidth', 1.5, 'DisplayName', ['S.E., M = ' num2str(M)]);
    plot(Ns, energy_diff_DG(i, :), '-s', 'LineWidth', 1.5, 'DisplayName', ['D.G., M = ' num2str(M)]);
    hold off;
    title(['Energy Conservation Error for M = ', num2str(M)]);
    xlabel('Polynomial Order (N)');
    ylabel('Energy Difference (|final - initial|)');
    legend('Location', 'Best');
    grid on;
end


%%


M=10;
N=16;
% Generate spatial grids for SE and DG methods
x = map_gllse(N, L, M);    % Spectral Element (S.E.) grid
xdg = map_gll(N, L, M);     % Discontinuous Galerkin (D.G.) grid

% Compute solutions
Udg = solve_DG(M, N, 1);    % D.G. solution
Use = solve_SE(M, N, 1);    % S.E. solution
u_abs = exp(-0.5 * ((x - x0).^2) / lambda0^2); % Exact solution (Gaussian)

% Handle duplicates in xdg by averaging Udg values at each unique x
[x_unique, ~, idx_unique] = unique(xdg);
Udg = accumarray(idx_unique, Udg, [], @mean);

% Plot the results
figure;
hold on;
plot(x_unique, Udg, '-s', 'LineWidth', 1.5, 'DisplayName', 'Discontinuous Galerkin (D.G.)');
plot(x, Use, '-o', 'LineWidth', 1.5, 'DisplayName', 'Spectral Element (S.E.)');
plot(x, u_abs, '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2], 'DisplayName', 'Exact Solution');
hold off;

% Add plot enhancements
title('CFL=1.5');
xlabel('Spatial Position (x)');
ylabel('Solution Value');
legend('Location', 'Best');
grid on;

%%

M = 10;        % Number of elements
Ns = [4, 8, 16];  % Polynomial orders to consider

% Loop over each N value to compute and plot solutions
for j = 1:length(Ns)
    N = Ns(j);

    % Spatial grids for S.E. and D.G.
    x_se = map_gllse(N, L, M);   % S.E. grid
    x_dg = map_gll(N, L, M);     % D.G. grid

    % Cone initial condition for both methods
    u0_se = max(0, 1 - s * abs(x_se - x0));
    u0_dg = max(0, 1 - s * abs(x_dg - x0));

    % Compute final solutions using S.E. and D.G.
    U_final_SE = solve_SE(M, N, 2); 
    U_final_DG = solve_DG(M, N, 2);  
    
    % Handle duplicates in D.G. grid by averaging U_final_DG at each unique x
    [x_dg_unique, ~, idx_dg] = unique(x_dg);
    U_final_DG = accumarray(idx_dg, U_final_DG, [], @mean);

    % Plot the results for each N
    figure;
    hold on;
    plot(x_dg_unique, U_final_DG, '-s', 'LineWidth', 1.5, 'DisplayName', 'D.G. Solution');
    plot(x_se, U_final_SE, '-o', 'LineWidth', 1.5, 'DisplayName', 'S.E. Solution');
    plot(x_se, u0_se, '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2], 'DisplayName', 'Initial Cone Condition');
    hold off;
    
    % Add titles and labels
    title(['S.E. and D.G. Solutions with Cone Initial Condition, M = 10, N = ', num2str(N)]);
    xlabel('Spatial Position (x)');
    ylabel('Solution Value');
    legend('Location', 'Best');
    grid on;
end




























