function [spline_vals] = cubic_spline(global_mesh, fvals, s_prime, xcheck)
% CUBIC_SPLINE Evaluates the cubic spline at given xcheck points.
%
% Inputs:
%   global_mesh - (1 x (n+1)) array of sorted mesh points.
%   f - Function handle for the original function.
%   s_prime - (1 x (n+1)) array of computed first derivatives at mesh points.
%   xcheck - (1 x m) array of x-values to evaluate the spline.
%   num_intervals - Number of intervals used in the spline (n).
%
% Outputs:
%   spline_vals - (1 x m) array of spline values at xcheck.

    % Compute h array (interval widths)
    h = diff(global_mesh);

    % Initialize output arrays
    spline_vals = zeros(size(xcheck));

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

        % Compute Hermite basis functions
        psi_L = ((x - x_right).^2 / h_i^2) .* (1 + 2 * (x - x_left) / h_i);
        psi_R = ((x - x_left).^2 / h_i^2) .* (1 - 2 * (x - x_right) / h_i);
        Psi_L = ((x - x_right).^2 / h_i^2) .* (x - x_left);
        Psi_R = ((x - x_left).^2 / h_i^2) .* (x - x_right);

        % Compute spline value using the Hermite basis
        spline_vals(j) = psi_L * f_left + psi_R * f_right + Psi_L * s_left + Psi_R * s_right;
    end
end
