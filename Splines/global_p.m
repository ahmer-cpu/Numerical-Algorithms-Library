function [pvals, ftruevals, lebesgue_vals, kappa_numer] = global_p(f, xcheck, a, b, n, node_type)
%GLOBAL_P Barycentric global interpolation over [a,b]
%
%   [pvals, ftruevals, lebesgue_vals, kappa_numer] = global_p(f, xcheck, a, b, n, node_type)
%
%   Inputs:
%       f         - Function handle to interpolate
%       xcheck    - Column vector of points to evaluate the interpolation at
%       a, b      - Interval [a, b]
%       n         - Degree of interpolating polynomial
%       node_type - 'chebyshev1' or 'chebyshev2'
%
%   Outputs:
%       pvals         - Interpolated polynomial values at xcheck
%       ftruevals     - True function values at xcheck
%       lebesgue_vals - Lebesgue function values at xcheck
%       kappa_numer   - Barycentric numerator sum at xcheck (Lebesgue kernel)

    % Generate interpolation nodes
    switch lower(node_type)
        case 'chebyshev1'
            xvals = chebyshev1(a, b, n);
        case 'chebyshev2'
            xvals = chebyshev2(a, b, n);
        otherwise
            error('Invalid node type. Use "chebyshev1" or "chebyshev2".');
    end

    % Sort and evaluate function at interpolation nodes
    xvals = sort(xvals);
    yvals = f(xvals);

    % True function values at evaluation points
    ftruevals = f(xcheck);

    % Barycentric weights
    gammainvs = bary1coeffs(xvals, n);

    % Preallocate output arrays
    kappa_numer = zeros(size(xcheck));
    lebesgue_vals = zeros(size(xcheck));

    % Evaluate barycentric interpolant
    [pvals, kappa_numer, lebesgue_vals] = bary1eval( ...
        xvals, yvals, gammainvs, n, xcheck, length(xcheck), ...
        kappa_numer, lebesgue_vals ...
    );
end
