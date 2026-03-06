function [pvals] = piecewise_p(xcheck, global_mesh, divdiffs, local_meshes, num_intervals, hermite_flag)
%PIECEWISE_P Evaluates piecewise Newton or Hermite interpolant at specified xcheck points
%
%   [pvals, ftrue_vals] = piecewise_p(f, xcheck, global_mesh, divdiffs, local_meshes, num_intervals, hermite_flag)
%
%   Inputs:
%       f             - Function handle to compute true values
%       xcheck        - Column vector of evaluation points
%       global_mesh   - Vector of interval endpoints (length = num_intervals + 1)
%       divdiffs      - Matrix of divided differences (size depends on interpolation type)
%       local_meshes  - Matrix of local interpolation nodes (size depends on interpolation type)
%       num_intervals - Number of subintervals
%       hermite_flag  - Boolean: true = Hermite interpolation, false = Newton
%
%   Outputs:
%       pvals         - Column vector of interpolated values at each xcheck point
%       ftrue_vals    - True function values at each xcheck point

    
    pvals = zeros(size(xcheck));

    for k = 1:length(xcheck)
        x = xcheck(k);

        % Find which interval x belongs to
        if x == global_mesh(end)
            i = num_intervals;
        else
            i = find(global_mesh <= x, 1, 'last');
            if i >= length(global_mesh)
                i = num_intervals;
            end
        end

        % Select corresponding local mesh and divided differences
        local_x = local_meshes(i, :);
        local_coeffs = divdiffs(i, :);

        % If Hermite, duplicate local_x as [a, a, b, b]
        if hermite_flag
            local_x = reshape([local_x; local_x], 1, []);
        end

        % Evaluate using nddeval
        pvals(k) = nddeval(local_coeffs, local_x, x);
    end
end
