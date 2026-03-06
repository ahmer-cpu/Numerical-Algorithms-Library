function [spline_deriv_vals] = cubic_spline_derivative(global_mesh, fvals, s_prime, xcheck)
% CUBIC_SPLINE_DERIVATIVE Evaluates the derivative of the cubic spline at given xcheck points.
%
% Inputs:
%   global_mesh - (1 x (n+1)) array of sorted mesh points.
%   fvals - (1 x (n+1)) array of function values at mesh points.
%   s_prime - (1 x (n+1)) array of computed first derivatives at mesh points.
%   xcheck - (1 x m) array of x-values to evaluate the spline derivative.
%
% Outputs:
%   spline_deriv_vals - (1 x m) array of spline derivative values at xcheck.

    % Compute h array (interval widths)
    h = diff(global_mesh);

    % Initialize output arrays
    spline_deriv_vals = zeros(size(xcheck));

    % Loop through each xcheck value
    for j = 1:length(xcheck)
        x = xcheck(j);

        % Find the relevant interval i where x belongs to [x_{i-1}, x_i]
        i = find(global_mesh(1:end-1) <= x & x <= global_mesh(2:end), 1, 'last');

        if isempty(i)
            error('xcheck value out of global mesh range.');
        end

        % Interval endpoints
        x_left = global_mesh(i);
        x_right = global_mesh(i+1);
        h_i = h(i);
        
        % Function values and derivatives
        f_left = fvals(i);
        f_right = fvals(i+1);
        s_left = s_prime(i);
        s_right = s_prime(i+1);

        % Compute derivatives of Hermite basis functions
        psi_L_deriv = (2 * (x - x_right) / h_i^2) .* (1 + 2 * (x - x_left) / h_i) + ((x - x_right).^2 / h_i^2) * (2 / h_i);
        psi_R_deriv = (2 * (x - x_left) / h_i^2) .* (1 - 2 * (x - x_right) / h_i) - ((x - x_left).^2 / h_i^2) * (2 / h_i);
        Psi_L_deriv = ((x - x_right).^2 / h_i^2) + (2 * (x - x_right) / h_i^2) .* (x - x_left);
        Psi_R_deriv = ((x - x_left).^2 / h_i^2) - (2 * (x - x_left) / h_i^2) .* (x - x_right);

        % Compute the spline derivative value using the derivative basis functions
        spline_deriv_vals(j) = psi_L_deriv * f_left + psi_R_deriv * f_right + Psi_L_deriv * s_left + Psi_R_deriv * s_right;
    end
end
