function y = lineartransform(x,a,b,c,d)
%LINEARTRANSFORM Linearly transform into the interval c,d from the interval
% a,b
%   This routine takes in the array of x values in the interval [a,b] and
%   applies a linear transform to all the values into the interval [c,d]
if a == b
        error('Original interval [a, b] must have a != b');
end
y = (x - a)*(d - c)/ (b - a) + c;
end

