function [xcheck, ftrue, pvals_global, pvals_piecewise, pvals_hermite, spline_vals1, spline_vals2] = run_all_interpolants( ...
    f, f_prime, a, b, ...
    n_global, node_type_global, ...
    num_intervals_piecewise, s_piecewise, node_type_piecewise, ...
    bd_type, boundary_conditions, ...
    ncheck)

    % Evaluation points
    xcheck = linspace(a, b, ncheck)';
    ftrue = f(xcheck);

    %% ---- Global Polynomial Interpolation ----
    [pvals_global, ~, ~, ~] = global_p(f, xcheck, a, b, n_global, node_type_global);
    pvals_global = pvals_global(:);  % Force column vector

    %% ---- Piecewise Polynomial Interpolation ----
    global_mesh = linspace(a, b, num_intervals_piecewise + 1);

    % Piecewise interpolation with user-defined degree
    hermite_flag = false;
    [divdiffs, local_meshes] = local_divdiffs(f, global_mesh, s_piecewise, node_type_piecewise, hermite_flag, f_prime);
    pvals_piecewise = piecewise_p(xcheck, global_mesh, divdiffs, local_meshes, num_intervals_piecewise, hermite_flag);

    % Hermite piecewise (fixed at degree 3 for now)
    hermite_flag = true;
    [divhermite, localhermite] = local_divdiffs(f, global_mesh, 3, 'uniform', hermite_flag, f_prime);
    pvals_hermite = piecewise_p(xcheck, global_mesh, divhermite, localhermite, num_intervals_piecewise, hermite_flag);

    %% ---- Cubic Spline 1 (First Derivative or Natural) ----
    fvals = f(global_mesh);
    s_prime = cubic_spline_first_derivatives(global_mesh, fvals, boundary_conditions, bd_type);
    spline_vals1 = cubic_spline(global_mesh, fvals, s_prime, xcheck);

    %% ---- Cubic Spline 2 (B-spline) ----
    h = (b - a) / num_intervals_piecewise;
    alpha = B_spline_alpha(global_mesh, fvals, boundary_conditions, h, bd_type);
    spline_vals2 = B_spline(xcheck, alpha, global_mesh);
end
