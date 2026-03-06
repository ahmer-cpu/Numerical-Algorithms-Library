function [divdiffs, local_meshes] = local_divdiffs(f, global_mesh, s, node_type, hermite_flag, fprime)
%LOCAL_DIVDIFFS Computes local Newton or Hermite divided differences for piecewise interpolation
%
%   [divdiffs, local_meshes] = local_divdiffs(f, global_mesh, s, node_type, hermite_flag, fprime)
%
%   Inputs:
%       f            - Function handle to interpolate
%       global_mesh  - Vector defining subinterval endpoints (length = num_intervals + 1)
%       s            - Degree of local Newton interpolants (ignored for Hermite, always cubic)
%       node_type    - 'chebyshev2' or 'uniform'
%       hermite_flag - Boolean: true = Hermite interpolation, false = standard Newton
%       fprime       - Function handle for f'(x), only used if hermite_flag = true
%
%   Outputs:
%       divdiffs     - Matrix, each row is local divided differences
%       local_meshes - Matrix, each row is the local interpolation nodes

    num_intervals = length(global_mesh) - 1;
    if hermite_flag
        % Hermite cubic: always 4 divided differences from 2 nodes
        divdiffs = zeros(num_intervals, 4);
        local_meshes = zeros(num_intervals, 2);  % Only 2 actual nodes (a and b)
    else
        % Newton: s+1 nodes and divided differences
        divdiffs = zeros(num_intervals, s+1);
        local_meshes = zeros(num_intervals, s+1);
    end
    
    for i = 1:num_intervals
        a = global_mesh(i);
        b = global_mesh(i+1);

        if hermite_flag
            % Use just the endpoints for Hermite interpolation
            local_nodes = [a, b];
            fvals_local = f(local_nodes);
            fprime_local = fprime(local_nodes);
            n = 4;  % Hermite cubic (with node duplication)

            div_row = nddcoeffs(local_nodes, fvals_local, fprime_local, n, true);
        else
            % Choose local interpolation nodes based on node_type
            switch lower(node_type)
                case 'chebyshev2'
                    local_nodes = chebyshev2(a, b, s);
                case 'uniform'
                    local_nodes = linspace(a, b, s+1);
                otherwise
                    error('Unsupported node type. Use "chebyshev2" or "uniform".');
            end

            local_nodes = sort(local_nodes); % Ensure ascending order
            fvals_local = f(local_nodes);
            n = s + 1;
            div_row = nddcoeffs(local_nodes, fvals_local, [], n, false);
        end

        divdiffs(i, :) = div_row;
        local_meshes(i, :) = local_nodes;
    end
end
