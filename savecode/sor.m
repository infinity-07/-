function [x, iter] = sor(A, b, omega, x0, tol, max_iter)
% Input parameters:
% A - Coefficient matrix
% b - Right-hand side vector
% omega - Relaxation factor
% x0 - Initial guess
% tol - Tolerance for the stopping criterion
% max_iter - Maximum number of iterations

% Output parameters:
% x - Solution vector
% iter - Number of iterations performed

if nargin < 5
    tol = 1e-8;
end

if nargin < 4
    x0 = ones(length(b), 1);
end

if nargin < 3
    omega = 1.25;  % Default relaxation factor
end

% Initialize variables
D = diag(diag(A));       % Diagonal matrix
U = -triu(A, 1);         % Upper triangular matrix
L = -tril(A, -1);        % Lower triangular matrix

% Compute iteration matrix and constant vector for SOR
L_w = (D - omega * L) \ ((1 - omega) * D + omega * U);
f = omega * ((D - omega * L) \ b);

x = x0;
iter = 0;

for k = 1:max_iter
    x = L_w * x + f;
    iter = iter + 1;

    % Check for convergence
    if norm(A * x - b) < tol
        break;
    end
end

if iter == max_iter
    warning('Maximum number of iterations reached without convergence.');
end

% Return the final solution and the number of iterations
x = double(x);
end