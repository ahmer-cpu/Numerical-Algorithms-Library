function x = chebyshev2(a,b,n)
j = 0:n-1;
x = cos((j * pi) / n);
x = lineartransform(x,-1,1,a,b);
end