function [p,kappa, kappa_one] = bary1Eval(xvals,fvals,nodes,gammas)
%BARY1EVAL Evaluates the interpolating polynomial in
%Barycentric 1 form.
%   This routine takes in fvals, nodes, and gammas and outputs a vector of
%   evaluations at the array xvals.

% Compute kappa and kappa_one in double precision first
    diffMat_d = xvals.' - nodes; % Default is double precision
    termMat_abs_d = abs(fvals .* gammas ./ diffMat_d);
    termMat_one_d = abs(gammas ./ diffMat_d);

    kappa_one = sum(termMat_one_d, 2);
    kappa = sum(termMat_abs_d, 2);

    omega_abs_d = abs(polyProduct(1, nodes, xvals)); % Compute omega in double precision
    kappa = kappa' .* omega_abs_d;
    kappa_one = kappa_one' .* omega_abs_d;

    % Now convert parameters to single precision for all other computations
    xvals = single(xvals);
    fvals = single(fvals);
    nodes = single(nodes);
    gammas = single(gammas);

    omega = single(polyProduct(1, nodes, xvals));
    diffMat = single(xvals.' - nodes);
    numVec = single(fvals .* gammas);
    termMat = single(numVec ./ diffMat);

    p = sum(termMat, 2);
    p = p' .* omega;

    tol = single(1e-12);
    dists = abs(single(repmat(xvals, length(nodes), 1)') - single(repmat(nodes, length(xvals), 1)));

    [minDist, minDistIdx] = min(dists, [], 2);
    closeMask = (minDist < tol);
    p(closeMask) = fvals(minDistIdx(closeMask));
end