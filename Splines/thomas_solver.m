function x = thomas_solver(a, b, c, d)
%THOMAS_SOLVER Solves a tridiagonal system A * x = d using the Thomas algorithm
%
%   x = thomas_solver(a, b, c, d)
%
%   Inputs:
%       a - (n-1)x1 subdiagonal vector (a(2) to a(n))
%       b - (n)x1 main diagonal vector
%       c - (n-1)x1 superdiagonal vector (c(1) to c(n-1))
%       d - (n)x1 right-hand side vector
%
%   Output:
%       x - (n)x1 solution vector

    n = length(b);

    % Ensure column vectors
    a = a(:);
    b = b(:);
    c = c(:);
    d = d(:);

    % Forward elimination
    for i = 2:n
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end

    % Back substitution
    x = zeros(n, 1);
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end
