close all; clear; clc;


% Parameters
N=24;
M=10;
Lx = 2; % Domain size in x-direction
Ly = 1; % Domain size in y-direction
Nx = 128; % Number of grid points in x-direction
x0 = -0.5; % Initial bump center in x
y0 = 0; % Initial bump center in y
lambda0 = 1/8; % Width of the Gaussian bump
kappa = 5e-4; % Diffusion coefficient

C=1;
dt=0.001;
L=1;

    [xi, w] = lglnodes(N + 1); 
    xi = flipud(xi);         
    % Differentiation matrix and stiffness/mass matrices
    D = derv(N, xi);
    y = map_gllse(N, 1, M)-0.5;  % Map nodes to the spatial domain
    SMass = SEmass(w, M, L);
    G = D'*diag(w)*D.* (2 * M)/L;
    Gglobalsym = makeGsym(G,M);
    Aad = inv(SMass) * Gglobalsym;
    Ny=length(Aad(:,1));
    A = dt*Aad*kappa+eye(Ny);
    A(1,1) = 1;
    A(end,end) = A(1,1);



% Grid
x = linspace(-Lx/2, Lx/2, Nx+1); % x-domain
[X, Y] = meshgrid(x(1:end-1), y); % 2D grid



%%
% Parameters
lambda0 = 0.1; % Width of the Gaussian
u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2)); % Gaussian initial condition
T = 2; % Time at which we want the solution

% Advection and diffusion parameters
c = 1; % Velocity in x-direction
kappa_x = 5e-4; % Diffusion coefficient in x
kappa_y = 5e-4; % Diffusion coefficient in y
dt = 0.001; % Time step
Nt = round(T / dt); % Number of time steps

% Wavenumbers for Fourier differentiation (1D in x)
kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1]; % Wavenumbers in x
Laplacian_x = -kx.^2; % Spectral operator for diffusion in x

% Transform initial condition to spectral space row by row
u_hat = zeros(Ny, Nx); % Initialize spectral solution
for j = 1:Ny
    u_hat(j, :) = fft(u0(j, :)); % Compute 1D FFT along each row
end

% Time stepping loop
for t = 1:Nt
    % Update advection and diffusion in x (spectral space)
    for j = 1:Ny
        u_hat(j, :) = u_hat(j, :) .* exp(-1i * c * kx * dt);         % Advection in x
        u_hat(j, :) = u_hat(j, :) .* exp(kappa_x * Laplacian_x * dt); % Diffusion in x
    end

    % Transform back to physical space for diffusion in y
    u_real = zeros(Ny, Nx); % Initialize real-space solution
    for j = 1:Ny
        u_real(j, :) = real(ifft(u_hat(j, :))); % Compute 1D IFFT along each row
    end

    % Update diffusion in y (real space)
    for i = 1:Nx
        u_real(:, i) = A \ u_real(:, i); % Solve the linear system for each column
    end

    % Transform updated solution back to spectral space
    for j = 1:Ny
        u_hat(j, :) = fft(u_real(j, :)); % Compute 1D FFT along each row
    end
end

% Transform final solution back to physical space
u_final = zeros(Ny, Nx);
for j = 1:Ny
    u_final(j, :) = real(ifft(u_hat(j, :))); % Compute 1D IFFT along each row
end

% Plot results
figure;

% Plot initial condition
subplot(1, 2, 1);
contourf(X, Y, u0, 20, 'LineStyle', 'none'); % 20 contour levels, no lines
title('Initial Condition');
xlabel('x');
ylabel('y');
colorbar;
axis equal tight;

% Plot final solution
subplot(1, 2, 2);
contourf(X, Y, u_final, 20, 'LineStyle', 'none'); % 20 contour levels, no lines
title('Solution with Diffusion in y and Advection-Diffusion in x T=2');
xlabel('x');
ylabel('y');
colorbar;
axis equal tight;

% Set the same color scale for both plots
u_min = min([u0(:); u_final(:)]); % Minimum value across both plots
u_max = max([u0(:); u_final(:)]); % Maximum value across both plots
subplot(1, 2, 1); caxis([u_min, u_max]); % Apply color axis to first plot
subplot(1, 2, 2); caxis([u_min, u_max]); % Apply color axis to second plot;



%%
t=T;

% Adjust X to enforce periodicity in the range [-1, 1]
X_periodic = mod(X - C * t - x0 + 1, 2) - 1; 
% Explanation:
% (X - C * t + 1): Shift range to [0, 2].
% mod(..., 2): Apply periodicity within [0, 2].
% -1: Shift back to [-1, 1].

% Exact solution at t = 1
factor = lambda0^2 / (lambda0^2 + 2 * kappa * t);
u_exact_t1 = factor * ...
             exp(-((X_periodic).^2 + (Y).^2) / (2 * (lambda0^2 + 2 * kappa * t)));



% Plot
figure;
contourf(X, Y, u_exact_t1-u_final, 50, 'LineColor', 'none'); % Filled contour plot
colorbar; % Add color bar
title(['Error at t = ', num2str(t), ': Gaussian Bump']);
xlabel('x');
ylabel('y');
grid on;
axis equal tight;

% Compute the L-infinity norm
L_inf_error = max(max(abs(u_final - u_exact_t1)));

% Display results
disp(['L-inf norm of the error at T = ', num2str(2), ' is: ', num2str(L_inf_error)]);


