function [x, iter] = preconditioned_conjugate_gradient(A, b, M, x0, tol, max_iter)
    % Input parameters:
    % A - Coefficient matrix
    % b - Right-hand side vector
    % M - Preconditioner matrix
    % x0 - Initial guess
    % tol - Tolerance for the stopping criterion
    % max_iter - Maximum number of iterations
    
    % Output parameters:
    % x - Solution vector
    % iter - Number of iterations performed
    
    % Initialization
    x = x0;
    r = b - A * x0;
    z = M \ r;
    d = z;
    
    for k = 0:max_iter-1
        if norm(r) < tol
            break;
        end
        
        alpha = (r' * z) / (d' * (A * d));
        x = x + alpha * d;
        r_new = r - alpha * A * d;
        
        if norm(r_new) < tol
            break;
        end
        
        z_new = M \ r_new;
        beta = (r_new' * z_new) / (r' * z);
        d = z_new + beta * d;
        
        % Update for the next iteration
        r = r_new;
        z = z_new;
    end
    
    iter = k + 1;
end