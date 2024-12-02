close all; clear; clc;


% Parameters
N=12;
M=10;
Lx = 2; % Domain size in x-direction
Ly = 1; % Domain size in y-direction
Nx = 128; % Number of grid points in x-direction
Ny = 120; % Number of grid points in y-direction
x0 = -0.5; % Initial bump center in x
y0 = 0; % Initial bump center in y
lambda0 = 1/8; % Width of the Gaussian bump
kappa = 5e-4; % Diffusion coefficient
t = 0; % Time at which we want the solution
C=1;
dt=0.001;
L=1;

    [xi, w] = lglnodes(N + 1); 
    xi = flipud(xi);         
    % Differentiation matrix and stiffness/mass matrices
    D = derv(N, xi);
    y = map_gllse(N, 1, M)-0.5;  % Map nodes to the spatial domain
    SMass = SEmass(w, M, L);


% Grid
x = linspace(-Lx/2, Lx/2, Nx+1); % x-domain
[X, Y] = meshgrid(x(1:end-1), y); % 2D grid

% Initial condition
u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2));

% Plot
figure;
contourf(X, Y, u0, 50, 'LineColor', 'none'); % Filled contour plot
colorbar; % Add color bar
title('Initial Condition: Gaussian Bump');
xlabel('x');
ylabel('y');
grid on;

%%


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
contourf(X, Y, u_exact_t1, 50, 'LineColor', 'none'); % Filled contour plot
colorbar; % Add color bar
title(['Exact Solution at t = ', num2str(t), ': Gaussian Bump']);
xlabel('x');
ylabel('y');
grid on;

%%

 % Parameters
% Lx = 2; % Domain length
% Nx = 128; % Number of grid points
% x = linspace(0, Lx, Nx+1); % Grid points (include endpoint for periodicity)
% x = x(1:end-1); % Remove duplicate point for periodic domain
% dx = Lx / Nx; % Grid spacing
% Q_values = [1, 4, 16]; % Wavelength counts
% 
% % Initialize error storage
% L_inf_errors = zeros(length(Q_values), 3); % Store errors for FFT, 1st, and 2nd derivatives
% 
% % Loop over Q values
% for q_idx = 1:length(Q_values)
%     Q = Q_values(q_idx);
% 
%     % Generate sine function for given Q
%     u = sin(2 * pi * Q * x / Lx);
% 
%     % Analytical derivatives
%     u_exact_1st = 2 * pi * Q/Lx * cos(2 * pi * Q * x / Lx); % 1st derivative
%     u_exact_2nd = -4 * pi^2 * Q^2/Lx^2 * sin(2 * pi * Q * x / Lx); % 2nd derivative
% 
%     % FFT of u
%     u_hat = fft(u);
% 
%     % Wavenumbers
%     kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1];
% 
%     % First derivative in Fourier space
%     u_hat_1st = 1i * kx .* u_hat;
%     u_num_1st = real(ifft(u_hat_1st)); % Back to physical space
% 
%     % Second derivative in Fourier space
%     u_hat_2nd = -(kx.^2) .* u_hat;
%     u_num_2nd = real(ifft(u_hat_2nd)); % Back to physical space
% 
%     % Inverse FFT to reconstruct the original function
%     u_reconstructed = ifft(u_hat);
% 
%     % Compute errors
%     L_inf_errors(q_idx, 1) = max(abs(u - u_reconstructed)); % Error in original reconstruction
%     L_inf_errors(q_idx, 2) = max(abs(u_num_1st - u_exact_1st)); % Error in 1st derivative
%     L_inf_errors(q_idx, 3) = max(abs(u_num_2nd - u_exact_2nd)); % Error in 2nd derivative
% end
% 
% % Display errors
% disp('L_inf errors for Q = 1, 4, 16:');
% disp(array2table(L_inf_errors, 'VariableNames', {'FFT_Error', '1st_Deriv_Error', '2nd_Deriv_Error'}, ...
%     'RowNames', {'Q=1', 'Q=4', 'Q=16'}));


%%
t=0;
t_final=0.1;
    G = D'*diag(w)*D.* (2 * M)/L;
    Gglobalsym = makeGsym(G,M);
    Aad = inv(SMass) * Gglobalsym;
    A = dt*Aad*kappa+eye(121);
    A(1,1) = 1;
    A(end,end) = A(1,1);

uhat=fft(u0);
while t<t_final
    % Wavenumbers
    kx = [0:Nx/2-1, -Nx/2:-1];
    uhatadv = -(1i * kx) .* uhat.*dt + uhat;
    uhat=uhatadv;
    % uhatdiff = ((ones(1,length(kx))-dt.*kx.^2).^-1).*uhatadv;
    % ureal = real(uhatdiff);
    % uimg = imag(uhatdiff);
    % 
    % 
    % for n=1:128
    %     u0r=ureal(:,n);
    %     u0i=uimg(:,n);
    %     uhatr = A\u0r;
    %     uhati = A\u0i;
    %     uhat(:,n) = uhatr+1i*uhati;
    % 
    % end

t = t + dt;
end

utrue=real(ifft(uhatadv));


% Plot
figure;
contourf(X, Y, utrue, 50, 'LineColor', 'none'); % Filled contour plot
colorbar; % Add color bar
title('Initial Condition: Gaussian Bump');
xlabel('x');
ylabel('y');
grid on;











%%
% 1D Advection Fourier-Galerkin Method
close all; clear; clc;
% Parameters
L = 2;           % Length of the domain
N = 128;         % Number of Fourier modes (grid points)
c = 1.0;         % Advection speed
T = 1;           % Total simulation time
dt = 0.0001;       % Time step size

% Spatial discretization
x = linspace(0, L, N+1); % Spatial grid (periodic boundary)
x = x(1:end-1);          % Remove duplicate point
dx = L / N;              % Grid spacing
k = [0:N/2-1, -N/2:-1] * (2*pi/L); % Fourier wavenumbers

% Initial condition
u0 = exp(-((x-L/2).^2) / (2*0.1^2)); % Gaussian pulse

% Fourier transform of the initial condition
u_hat = fft(u0);

% Time-stepping loop
u_hat_new = u_hat; % Initialize
time = 0;
while time < T
    % Compute the Fourier derivative

    % Update in Fourier space using the explicit Euler method
    u_hat_new = u_hat - dt * c * 1i * k .* u_hat;

    % Update time and solution
    time = time + dt;
    u_hat = u_hat_new;
end

% Inverse Fourier transform to obtain the solution in physical space
u = ifft(u_hat);

% Visualization
figure;
plot(x, u0, 'k--', 'LineWidth', 1.5); hold on;
plot(x, real(u), 'b-', 'LineWidth', 1.5);
xlabel('x'); ylabel('u');
title(['1D Advection using Fourier-Galerkin, Time = ', num2str(T)]);
legend('Initial Condition', 'Numerical Solution');
grid on;

%%
% Parameters
Lx = 2; % Domain length
Nx = 128; % Number of grid points
x = linspace(0, Lx, Nx+1); % Grid points (include endpoint for periodicity)
x = x(1:end-1); % Remove duplicate point for periodic domain
dx = Lx / Nx; % Grid spacing
sigma = 0.1; % Width of the Gaussian
u0 = exp(-((x - Lx/2).^2) / (2 * sigma^2)); % Initial condition

% Advection parameters
c = 1; % Velocity
dt = 0.0001; % Time step
T = 2; % Total simulation time
Nt = round(T / dt); % Number of time steps

% Wavenumbers for Fourier differentiation
kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1];

% Initial condition in spectral space
u_hat = fft(u0);

% Time stepping loop
for t = 1:Nt
    % Update in spectral space using advection equation
    u_hat = u_hat .* (1-1i * c * kx * dt);  
end

% Transform back to physical space
u_final = real(ifft(u_hat));

% Plot results
figure;
plot(x, u0, 'b', 'LineWidth', 1.5); hold on;
plot(x, u_final, 'r--', 'LineWidth', 1.5);
title('1D Advection Equation Solution (Spectral Method)');
legend('Initial Condition', 'Final Solution');
xlabel('x');
ylabel('u(x)');
grid on;

%%
% Parameters
% Lx = 2; % Domain length
% Nx = 128; % Number of grid points
% x = linspace(0, Lx, Nx+1); % Grid points (include endpoint for periodicity)
% x = x(1:end-1); % Remove duplicate point for periodic domain
% 
% % Initial condition
% sigma = 0.1; % Width of the Gaussian
% u0 = exp(-((x - Lx/2).^2) / (2 * sigma^2)); % Gaussian function
% 
% % Analytical derivatives
% u_exact_1st = -((x - Lx/2) / sigma^2) .* u0; % 1st derivative
% u_exact_2nd = ((x - Lx/2).^2 / sigma^4 - 1 / sigma^2) .* u0; % 2nd derivative
% 
% % FFT of u0
% u_hat = fft(u0);
% 
% % Wavenumbers
% kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1];
% 
% % Numerical derivatives using FFT
% u_hat_1st = 1i * kx .* u_hat;
% u_num_1st = real(ifft(u_hat_1st)); % Back to physical space
% 
% u_hat_2nd = -(kx.^2) .* u_hat;
% u_num_2nd = real(ifft(u_hat_2nd)); % Back to physical space
% 
% % Plot comparisons
% figure('Name', 'Gaussian Initial Condition Derivatives', 'NumberTitle', 'off');
% 
% % Original function comparison
% subplot(3, 1, 1);
% plot(x, u0, 'b', 'LineWidth', 1.5); hold on;
% plot(x, real(ifft(u_hat)), 'r--', 'LineWidth', 1.5);
% title('Original Function Comparison');
% legend('Analytical', 'Numerical (FFT Reconstruction)');
% xlabel('x');
% ylabel('u(x)');
% grid on;
% 
% % First derivative comparison
% subplot(3, 1, 2);
% plot(x, u_exact_1st, 'b', 'LineWidth', 1.5); hold on;
% plot(x, u_num_1st, 'r--', 'LineWidth', 1.5);
% title('1st Derivative Comparison');
% legend('Analytical', 'Numerical');
% xlabel('x');
% ylabel("u'(x)");
% grid on;
% 
% % Second derivative comparison
% subplot(3, 1, 3);
% plot(x, u_exact_2nd, 'b', 'LineWidth', 1.5); hold on;
% plot(x, u_num_2nd, 'r--', 'LineWidth', 1.5);
% title('2nd Derivative Comparison');
% legend('Analytical', 'Numerical');
% xlabel('x');
% ylabel("u''(x)");
% grid on;

% Adjust spacin



%%
% % Parameters
% Lx = 2; % Domain size in x-direction
% Ly = 1; % Domain size in y-direction
% Nx = 128; % Number of grid points in x-direction
% Ny = 120; % Number of grid points in y-direction
% x = linspace(0, Lx, Nx+1); x = x(1:end-1); % x-grid (periodic domain)
% y = linspace(0, Ly, Ny+1); y = y(1:end-1); % y-grid (periodic domain)
% dx = Lx / Nx; % Grid spacing in x
% dy = Ly / Ny; % Grid spacing in y
% 
% [X, Y] = meshgrid(x, y); % 2D grid
% 
% % Initial condition
% x0 = Lx / 2; % Center in x
% y0 = Ly / 2; % Center in y
% lambda0 = 0.1; % Width of the Gaussian
% u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2)); % Gaussian initial condition
% 
% % Advection parameters
% c = 1; % Velocity in x-direction
% dt = 0.001; % Time step
% T = 1; % Total simulation time
% Nt = round(T / dt); % Number of time steps
% 
% % Wavenumbers for Fourier differentiation
% kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1]; % Wavenumbers in x
% ky = (2 * pi / Ly) * [0:Ny/2-1, -Ny/2:-1]; % Wavenumbers in y
% [KX, ~] = meshgrid(kx, ky); % 2D wavenumber grid (only KX matters for advection in x)
% 
% % Initial condition in spectral space
% u_hat = fft2(u0); % 2D FFT of initial condition
% 
% % Time stepping loop
% for t = 1:Nt
%     % Update in spectral space using the advection equation in x-direction
%     u_hat = u_hat .* exp(-1i * c * KX * dt);
% end
% 
% % Transform back to physical space
% u_final = real(ifft2(u_hat));
% 
% % Plot results
% figure;
% subplot(1, 2, 1);
% imagesc(x, y, u0);
% title('Initial Condition');
% xlabel('x');
% ylabel('y');
% colorbar;
% axis equal tight;
% 
% subplot(1, 2, 2);
% imagesc(x, y, u_final);
% title('Final Solution');
% xlabel('x');
% ylabel('y');
% colorbar;
% axis equal tight;

%%
% Parameters
Lx = 2; % Domain size in x-direction
Ly = 1; % Domain size in y-direction
Nx = 128; % Number of grid points in x-direction
Ny = 120; % Number of grid points in y-direction
x = linspace(0, Lx, Nx+1); x = x(1:end-1); % x-grid (periodic domain)
y = linspace(0, Ly, Ny+1); y = y(1:end-1); % y-grid (periodic domain)
dx = Lx / Nx; % Grid spacing in x
dy = Ly / Ny; % Grid spacing in y

[X, Y] = meshgrid(x, y); % 2D grid

% Initial condition
x0 = Lx / 2; % Center in x
y0 = Ly / 2; % Center in y
lambda0 = 0.1; % Width of the Gaussian
u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2)); % Gaussian initial condition

% Advection and diffusion parameters
c = 1; % Velocity in x-direction
kappa = 5e-4; % Diffusion coefficient in x
dt = 0.001; % Time step
T = 2; % Total simulation time
Nt = round(T / dt); % Number of time steps

% Wavenumbers for Fourier differentiation
kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1]; % Wavenumbers in x
ky = (2 * pi / Ly) * [0:Ny/2-1, -Ny/2:-1]; % Wavenumbers in y
[KX, KY] = meshgrid(kx, ky); % 2D wavenumber grid

% Precompute diffusion operator in x-direction
Laplacian_x = -KX.^2; % Spectral operator for diffusion in x

% Initial condition in spectral space
u_hat = fft2(u0); % 2D FFT of initial condition

% Time stepping loop
for t = 1:Nt
    % Update in spectral space using advection and diffusion in x
    u_hat = u_hat .* (1-1i * c * KX * dt);         % Advection term
    u_hat = u_hat .* (1+kappa * Laplacian_x * dt);  % Diffusion in x

end

% Transform back to physical space
u_final = real(ifft2(u_hat));

% Plot results
figure;
subplot(1, 2, 1);
imagesc(x, y, u0);
title('Initial Condition');
xlabel('x');
ylabel('y');
colorbar;
axis equal tight;

subplot(1, 2, 2);
imagesc(x, y, u_final);
title('Final Solution with Diffusion in x');
xlabel('x');
ylabel('y');
colorbar;
axis equal tight;






































