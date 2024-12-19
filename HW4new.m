close all; clc; clear;

% Number of modal coefficients
N = 16;
M=64;

% Initialize modal coefficients
a = ones(N+1, 1); % a_i = 1 for i = 0 to N

% Generate GLL points and weights (using MATLAB or external library)
[xi, w] = lglnodes(N + 1); 
x = flipud(xi); 

[xiM, wM] = lglnodes(M + 1); 
xM = flipud(xiM); 

% Compute Legendre polynomials at GLL points
P = zeros(N+1, length(x));
for i = 0:N
    P(i+1, :) = legendreP(i, x);
end

% Evaluate u(x) at GLL points
u = zeros(size(x));
for i = 0:N
    u(i+1) = sum(P(:,i+1));
end


u2 = zeros(size(x));
for i = 0:N
    u2(i+1) = sum(P(:,i+1))*sum(P(:,i+1));
end

% Plot the result
figure;
plot(x, u, '-o', 'LineWidth', 2, 'MarkerSize', 8);
title('Modal Representation of u(x) on GLL Grid N=16');
xlabel('x (GLL Points)');
ylabel('u(x)');
grid on;

% Generate the Legendre Vandermonde matrix
V = legendre_vandermonde(x, N+1);
BM = legendre_vandermonde(xM, M+1);

a = inv(V)*u2;
mod = linspace(0,16,17);
figure;
plot(mod, abs(a), '-o', 'LineWidth', 2, 'MarkerSize', 8);
title('the modal coefficients u(x)^2 N=16');
xlabel('Pn');
ylabel('coe');
grid on;

B=V;
% W = diag(w);
% A = inv(B.'*W*B);
% a = A*u;

%%
%N to M

coeN = inv(B)*u;
coeM = zeros(M+1,1);
for i = 1:length(u)
    coeM(i) = coeN(i);
end

uM = BM*coeM;

figure; % Create a new figure
hold on; % Hold on to overlay plots on the same graph
plot(x, u, '-o', 'LineWidth', 2, 'MarkerSize', 8); % Plot the first dataset
plot(xM, uM, '-r', 'LineWidth', 2, 'MarkerSize', 8); % Plot the second dataset
title('Modal Representation of u(x) on GLL Grid N=16 and M=32'); % Add title
xlabel('x (GLL Points)'); % Label for x-axis
ylabel('u(x)'); % Label for y-axis
grid on; % Turn on the grid
hold off; % Release hold

S = uM.^2;

coeSM = inv(BM)*S;
coeSN=zeros(size(x));
for i = 1:length(u)
    coeSN(i) = coeSM(i);
end

SN = B*coeSN;

figure;
plot(mod, abs(coeSN), '-o', 'LineWidth', 2, 'MarkerSize', 8);
title('the modal coefficients u(x)^2 N=16 after filter');
xlabel('Pn');
ylabel('coe');
grid on;

figure; % Create a new figure
hold on; % Hold on to overlay plots on the same graph
plot(x, u2, '-o', 'LineWidth', 2, 'MarkerSize', 8); % Plot the first dataset
plot(x, SN, '-r', 'LineWidth', 2, 'MarkerSize', 8); % Plot the second dataset
title('Modal Representation of S on GLL Grid'); % Add title
xlabel('x (GLL Points)'); % Label for x-axis
ylabel('u(x)'); % Label for y-axis
legend('S', 'After filter', 'Location', 'Best'); % Add legend
grid on; % Turn on the grid
hold off; % Release hold









