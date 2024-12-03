%%

% % Parameters
% % Lx = 2; % Domain size in x-direction
% % Ly = 1; % Domain size in y-direction
% % Nx = 128; % Number of grid points in x-direction
% % Ny = 120; % Number of grid points in y-direction
% % x = linspace(0, Lx, Nx+1); x = x(1:end-1); % x-grid (periodic domain)
% % y = linspace(0, Ly, Ny+1); y = y(1:end-1); % y-grid (periodic domain)
% % dx = Lx / Nx; % Grid spacing in x
% % dy = Ly / Ny; % Grid spacing in y
% % 
% % [X, Y] = meshgrid(x, y); % 2D grid
% % 
% % % Initial condition
% % x0 = Lx / 2; % Center in x
% % y0 = Ly / 2; % Center in y
% % lambda0 = 0.1; % Width of the Gaussian
% % u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2)); % Gaussian initial condition
% % 
% % % Advection and diffusion parameters
% % c = 1; % Velocity in x-direction
% kappa = 5e-3; % Diffusion coefficient in x
% dt = 0.001; % Time step
% T = 1; % Total simulation time
% Nt = round(T / dt); % Number of time steps
% % 
% % % Wavenumbers for Fourier differentiation
% kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1]; % Wavenumbers in x
% 
% % Precompute diffusion operator in x-direction
% Laplacian_x = -kx.^2; % Spectral operator for diffusion in x
% 
% % Initial condition in spectral space
% for n =1:128
%     u_hat(:,n) = fft(u0(:,n));
% end
% %u_hat = fft2(u0); % 2D FFT of initial condition
% 
% % Time stepping loop
% for t = 1:Nt
%     % Update in spectral space using advection and diffusion in x
%     u_hat = u_hat .* (1-1i * C * kx* dt);         % Advection term
%     u_hat = u_hat .* (1+kappa * Laplacian_x * dt);  % Diffusion in x
%     ureal = real(u_hat);
%     uimg = imag(u_hat);
%         % for n=1:128
%         %     u0r=ureal(:,n);
%         %     u0i=uimg(:,n);
%         %     uhatr = A\u0r;
%         %     uhati = A\u0i;
%         %     u_hat(:,n) = uhatr+1i*uhati;
%         % end
% 
% end
% 
% % Transform back to physical space
% %u_final = real(ifft2(u_hat));
% for n =1:128
%     u_final(:,n) = real(fft(u_hat(:,n)));
% end
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
% title('Final Solution with Diffusion in x');
% xlabel('x');
% ylabel('y');
% colorbar;
% axis equal tight;

%%
% Parameters
% u0 = exp(-((X - x0).^2 + (Y - y0).^2) / (2 * lambda0^2)); % Gaussian initial condition
% 
% % Advection and diffusion parameters
% c = 1; % Velocity in x-direction
% kappa = 5e-4; % Diffusion coefficient in x
% dt = 0.001; % Time step
% T = 0.01; % Total simulation time
% Nt = round(T / dt); % Number of time steps
% 
% % Wavenumbers for Fourier differentiation (1D)
% kx = (2 * pi / Lx) * [0:Nx/2-1, -Nx/2:-1]; % Wavenumbers in x
% Laplacian_x = -kx.^2; % Spectral operator for diffusion in x
% 
% % Transform initial condition to spectral space row by row
% u_hat = zeros(Ny, Nx); % Initialize spectral solution
% for j = 1:Ny
%     u_hat(j, :) = fft(u0(j, :)); % Compute 1D FFT along each row
% end
% 
% % Time stepping loop in spectral space
% for t = 1:Nt
%     % Update each row in spectral space
%     for j = 1:Ny
%         u_hat(1,:) = 0;
%         u_hat(end,:) = 0;
%         u_hat(j, :) = u_hat(j, :) .* exp(-1i * c * kx * dt);         % Advection term
%         u_hat(j, :) = u_hat(j, :) .* exp(kappa * Laplacian_x * dt);  % Diffusion in x
%             ureal = real(u_hat);
%             uimg = imag(u_hat);
%         for n=1:128
%             u0r=ureal(:,n);
%             u0i=uimg(:,n);
%             uhatr = A\u0r;
%             uhati = A\u0i;
%             u_hat(:,n) = uhatr+1i*uhati;
%         end
%     end
% end
% 
% % Transform back to physical space row by row
% u_final = zeros(Ny, Nx); % Initialize physical space solution
% for j = 1:Ny
%     u_final(j, :) = real(ifft(u_hat(j, :))); % Compute 1D IFFT along each row
% end
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
% title('Final Solution with Diffusion in x');
% xlabel('x');
% ylabel('y');
% colorbar;
% axis equal tight;

%%

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