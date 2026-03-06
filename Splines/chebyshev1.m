function [x] = chebyshev1(a,b,n)

x = zeros(1,n+1);
for j = 0:n
    x(j+1) = cos( (2.0*j + 1.0)*pi / (2.0*n + 2) );
end

x = lineartransform(x,-1,1,a,b);
end
