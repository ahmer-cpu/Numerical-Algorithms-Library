function [divdiffs] = nddcoeffs(nodes, fvals, f_prime_vals, n, hermite_flag)
%NDDCOEFFS Compute divided difference table
%   The function takes as input function y values and x values and
%   outputs the top row/left column of the divided differences table for
%   interpolation assuming the x values are indexed from top to bottom in the array.
%   If hermite_flag is true, the function uses Hermite interpolation by incorporating
%   first derivative values into the divided difference computation.

    if hermite_flag
        % **Hermite Interpolation: Duplicate nodes and function values**
        nodes = reshape([nodes; nodes], 1, []);
        fvals = reshape([fvals; fvals], 1, []);

        % Initialize divided differences array
        divdiffs = fvals;

        % **First iteration: Use f' values for repeated nodes**
        i = 1;
        divdiffs(i+1:end) = divdiffs(i+1:end) - divdiffs(i:end-1);
        D = nodes(i+1:end) - nodes(1:end-i);

        % Replace zero denominators with derivative values
        zero_idx = find(D == 0);
        if ~isempty(zero_idx)
            divdiffs(zero_idx + i) = f_prime_vals((zero_idx + 1) / 2);
            D(zero_idx) = 1;  % Prevent division by zero (not needed but ensures safety)
        end

        divdiffs(i+1:end) = divdiffs(i+1:end) ./ D;

        % **Subsequent iterations: Standard divided difference**
        for i = 2:n-1
            divdiffs(i+1:end) = divdiffs(i+1:end) - divdiffs(i:end-1);
            D = nodes(i+1:end) - nodes(1:end-i);
            divdiffs(i+1:end) = divdiffs(i+1:end) ./ D;
        end

    else
        % **Standard Newton Interpolation**
        divdiffs = fvals;
        for i = 1:n-1
            divdiffs(i+1:end) = divdiffs(i+1:end) - divdiffs(i:end-1);
            D = nodes(i+1:end) - nodes(1:end-i);
            divdiffs(i+1:end) = divdiffs(i+1:end) ./ D;
        end
    end
end
