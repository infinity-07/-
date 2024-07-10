function [x, iter] = jacobi(A, b, x0, tol, max_iter)
% Input parameters:
% A - Coefficient matrix
% b - Right-hand side vector
% x0 - Initial guess
% tol - Tolerance for the stopping criterion
% max_iter - Maximum number of iterations

% Output parameters:
% x - Solution vector
% iter - Number of iterations performed

if nargin < 4
    tol = 1e-8;
end

if nargin < 3
    x0 = ones(length(b), 1);
end

% Initialize variables
D = diag(diag(A));     % Diagonal matrix
U = -triu(A, 1);       % Upper triangular matrix
L = -tril(A, -1);      % Lower triangular matrix

B = D \ (L + U);
f = D \ b;

x = x0;
iter = 0;

for k = 1:max_iter
    x = B * x + f;
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