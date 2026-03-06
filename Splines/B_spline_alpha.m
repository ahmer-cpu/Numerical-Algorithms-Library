function alpha = B_spline_alpha(global_mesh, fvals, boundary_conditions, h, bc_type)
% B_SPLINE_ALPHA Computes the B-spline alpha parameters by solving a linear system.
%
% Inputs:
%   global_mesh - (1 x (n+1)) array of sorted mesh points.
%   f - Function handle for function evaluation.
%   boundary_conditions - (1x2) array specifying:
%       * [f_0', f_n'] for 'first' derivative boundary conditions
%       * [f_0'', f_n''] for 'second' derivative boundary conditions
%   h - Uniform interval length.
%   bc_type - String ('first' or 'second') indicating boundary condition type.
%
% Output:
%   alpha - (1 x (n+3)) array of computed B-spline parameters.

    % Number of intervals
    np1 = length(global_mesh); % np1 = n+1
    n = np1 - 1; % Number of intervals
    N = n + 3; % Total system size

    % Construct b vector (same for both cases)
    b = [boundary_conditions(1); fvals(:); boundary_conditions(2)]; % Ensure column vector

    % Construct A matrix
    A = zeros(N, N);

    if strcmp(bc_type, 'first')
        % First row: [-3/h, 0, 3/h, 0, ..., 0]
        A(1,1) = -3/h;
        A(1,3) = 3/h;
        % Last row: [0, ..., -3/h, 0, 3/h]
        A(N,N-2) = -3/h;
        A(N,N)   = 3/h;
    elseif strcmp(bc_type, 'second')
        % First row: [6/h^2, -12/h^2, 6/h^2, 0, ..., 0]
        A(1,1) = 6/h^2;
        A(1,2) = -12/h^2;
        A(1,3) = 6/h^2;
        % Last row: [0, ..., 6/h^2, -12/h^2, 6/h^2]
        A(N,N-2) = 6/h^2;
        A(N,N-1) = -12/h^2;
        A(N,N)   = 6/h^2;

    else
        error('Invalid boundary condition type. Use "first" or "second".');
    end
    
    for i = 2:N-1
        A(i,i-1) = 1; % Subdiagonal
        A(i,i)   = 4; % Main diagonal
        A(i,i+1) = 1; % Superdiagonal
    end
    alpha = A \ b;

end
