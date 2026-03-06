function fvals = function4(xvals)
%FUNCTION4 Outputs the values of function 4
%   This routine takes in the xvals array and outputs fvals = 1/(1+25*x**2)
fvals = (xvals).^2;
fvals = fvals *25;
fvals = fvals +1;
fvals = 1 ./ fvals;
end

