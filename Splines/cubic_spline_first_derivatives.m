function s_prime = cubic_spline_first_derivatives(global_mesh, fvals, boundary_conditions, bc_type)
% CUBIC_SPLINE_FIRST_DERIVATIVES computes first derivative parameters for a cubic spline.
%
% Inputs:
%   global_mesh - (1 x (n+1)) array of mesh points.
%   f - Function handle for the given function.
%   boundary_conditions - (1x2) array specifying either:
%       * [s_0', s_n'] for 'first' derivative boundary conditions
%       * [f_0'', f_n''] for 'second' derivative boundary conditions
%   bc_type - String ('first' or 'second') indicating which boundary condition to use.
%
% Output:
%   s_prime - (1 x (n+1)) array of computed first derivatives for the cubic spline.

    % Number of points
    np1 = length(global_mesh); % np1 = n+1
    n = np1 - 1; % Number of intervals
    
    % Compute h array
    h = diff(global_mesh);

    % Compute lambda and mu
    lambda = h(2:end) ./ (h(1:end-1) + h(2:end));
    mu = h(1:end-1) ./ (h(1:end-1) + h(2:end));

    % Compute first-order divided differences
    div_diff = diff(fvals) ./ h; % f[x_i, x_{i+1}]

    % Compute g array
    g = 3 * (lambda .* div_diff(1:end-1) + mu .* div_diff(2:end));

    % System Setup Based on Boundary Condition Type
    if strcmp(bc_type, 'first')
        % **First Derivative Boundary Conditions**
        s_0_prime = boundary_conditions(1);
        s_n_prime = boundary_conditions(2);

        % Modify g for first derivative case
        g(1) = g(1) - lambda(1) * s_0_prime;
        g(end) = g(end) - mu(end) * s_n_prime;

        % Define system diagonals
        main_diag = 2 * ones(n-1, 1);
        sub_diag = lambda(2:end); % Lower diagonal (size n-2)
        super_diag = mu(1:end-1); % Upper diagonal (size n-2)

        % Solve the system using the Thomas algorithm
        s_inner = thomas_solver(sub_diag, main_diag, super_diag, g); % Size (n-1)

        % Append boundary conditions
        s_prime = [s_0_prime; s_inner; s_n_prime]; % Final array of size (n+1)

    elseif strcmp(bc_type, 'second')
        % **Second Derivative Boundary Conditions**
        f_0_double_prime = boundary_conditions(1);
        f_n_double_prime = boundary_conditions(2);

        % Ensure g_0 and g_n are **scalars**
        g_0 = (-h(1) / 2) * f_0_double_prime + 3 * div_diff(1);
        g_n = (h(end) / 2) * f_n_double_prime + 3 * div_diff(end);
        
        % Ensure g is a **column vector**
        g = [g_0; g(:); g_n]; % Fix concatenation issue

        % Define system diagonals
        main_diag = 2 * ones(n+1, 1);
        sub_diag = [lambda(:); 1]; % Append 1 to the end
        super_diag = [1; mu(:)];   % Append 1 to the beginning

        % Solve the system using the Thomas algorithm
        s_prime = thomas_solver(sub_diag, main_diag, super_diag, g); % Final array of size (n+1)

    else
        error('Invalid boundary condition type. Use "first" or "second".');
    end

end
