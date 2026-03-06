function [x] = chebyshev2(a,b,n)
x = zeros(1,n+1);
for j = 0:n
    x(j+1) = cos((j * pi) / n);
end

x = lineartransform(x,-1,1,a,b);
end