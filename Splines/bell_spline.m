function B_val = bell_spline(xcheck, x_ref, h)
% BELL_SPLINE Evaluates the B-spline centered at x_ref at a given xcheck.
%
% Inputs:
%   xcheck - The point where we evaluate the B-spline.
%   x_ref - The reference center point x_i.
%   h - Uniform interval spacing.
%
% Output:
%   B_val - The computed B-spline value at xcheck.

    % Compute the five reference points
    x_m2 = x_ref - 2*h; % x_{i-2}
    x_m1 = x_ref - h;   % x_{i-1}
    x_0  = x_ref;       % x_i (center)
    x_1  = x_ref + h;   % x_{i+1}
    x_2  = x_ref + 2*h; % x_{i+2}

    % Initialize B-spline value
    B_val = 0;

    % Evaluate the B-spline based on the interval conditions
    if x_m2 <= xcheck && xcheck < x_m1
        B_val = (1/h^3) * (xcheck - x_m2)^3;

    elseif x_m1 <= xcheck && xcheck < x_0
        B_val = (1/h^3) * (h^3 + 3*h^2*(xcheck - x_m1) + 3*h*(xcheck - x_m1)^2 - 3*(xcheck - x_m1)^3);

    elseif x_0 <= xcheck && xcheck < x_1
        B_val = (1/h^3) * (h^3 + 3*h^2*(x_1 - xcheck) + 3*h*(x_1 - xcheck)^2 - 3*(x_1 - xcheck)^3);

    elseif x_1 <= xcheck && xcheck <= x_2
        B_val = (1/h^3) * (x_2 - xcheck)^3;
    end

end
