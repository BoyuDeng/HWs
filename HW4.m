close all; clc; clear;

% Number of modal coefficients
N = 16;

% Initialize modal coefficients
a = ones(N+1, 1); % a_i = 1 for i = 0 to N

% Generate GLL points and weights (using MATLAB or external library)
[xi, w] = lglnodes(N + 1); 
x = flipud(xi); 

% Compute Legendre polynomials at GLL points
P = zeros(N+1, length(x));
for i = 0:N
    P(i+1, :) = legendreP(i, x);
end

% Evaluate u(x) at GLL points
u = zeros(size(x));
for i = 0:N
    u(i+1) = sum(P(:,i+1))*sum(P(:,i+1));
end

% % Plot the result
% figure;
% plot(x, u, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% title('Modal Representation of u(x) on GLL Grid');
% xlabel('x (GLL Points)');
% ylabel('u(x)');
% grid on;

% Generate the Legendre Vandermonde matrix
V = legendre_vandermonde(x, N+1);

a = V\u;






