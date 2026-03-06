function [spline_vals] = B_spline(xcheck, alpha, global_mesh)
% B_SPLINE Evaluates the cubic B-spline at given xcheck points.
%
% Inputs:
%   xcheck - (1 x m) array of x-values to evaluate the spline.
%   alpha - (1 x (n+3)) array of B-spline coefficients (includes alpha_{-1} and alpha_{n+1}).
%   global_mesh - (1 x (n+1)) array of sorted mesh points.
%   f - Function handle for the original function.
%
% Outputs:
%   spline_vals - (1 x m) array of spline values at xcheck.
%   true_vals - (1 x m) array of true function values f(xcheck).

    % Number of intervals
    np1 = length(global_mesh); % np1 = n+1
    n = np1 - 1; % Number of intervals
    h = global_mesh(2) - global_mesh(1); % Uniform spacing

    % Extend the global mesh to include x_{-1} and x_{n+1}
    extended_mesh = [global_mesh(1) - h, global_mesh, global_mesh(end) + h];

    % Initialize output arrays
    spline_vals = zeros(size(xcheck));
    % Loop through each xcheck value
    for j = 1:length(xcheck)
        x = xcheck(j);

        % Check if xcheck exactly matches a global mesh point
        idx_exact = find(global_mesh == x, 1);

        if ~isempty(idx_exact)
            % Special case: xcheck is exactly at a mesh point
            i = idx_exact + 1; % Adjust for extended mesh indexing

            % Compute three contributing B-splines
            B_left = bell_spline(x, extended_mesh(i-1), h); % Left B-spline
            B_center = bell_spline(x, extended_mesh(i), h); % Center B-spline
            B_right = bell_spline(x, extended_mesh(i+1), h); % Right B-spline

            % Compute the weighted sum of B-splines
            spline_vals(j) = alpha(i-1) * B_left + alpha(i) * B_center + alpha(i+1) * B_right;

        else
            % General case: xcheck lies in an interval between mesh points
            i = find(extended_mesh(1:end-1) <= x & x <= extended_mesh(2:end), 1, 'last');

            if isempty(i)
                error('xcheck value out of global mesh range.');
            end

            % Compute the four contributing B-splines
            B_left_2 = bell_spline(x, extended_mesh(i-1), h); % Second left
            B_left_1 = bell_spline(x, extended_mesh(i), h);   % First left
            B_right_1 = bell_spline(x, extended_mesh(i+1), h); % First right
            B_right_2 = bell_spline(x, extended_mesh(i+2), h); % Second right

            % Compute the weighted sum of B-splines using alpha coefficients
            spline_vals(j) = alpha(i-1) * B_left_2 + alpha(i) * B_left_1 + ...
                             alpha(i+1) * B_right_1 + alpha(i+2) * B_right_2;
        end
    end
end
