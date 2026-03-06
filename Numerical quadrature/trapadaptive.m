function [integrals, errors, estimates] = trapadaptive(func,mesh,H,n,true_int, ref_num)
%TRAPADAPTIVE Summary of this function goes here
%   Detailed explanation goes here

% Initializing results
integrals = zeros(1,ref_num+1);
errors = zeros(1,ref_num+1);
estimates = zeros(1,ref_num);

% Starting recurrence
[int_start, error_start] = trapezoidal(func,mesh,H,n,true_int);
integrals(1) = int_start;
errors(1) = error_start;

int = int_start;
h_shift = H;
mesh_current = mesh;

for i = 2:ref_num+1
    h_old = h_shift;
    
    % New mesh values (same as evaluation values since H = h)
    h_shift = h_shift/2.0;
    mesh_new = mesh_current(1:n) + h_shift;
    
    % Complete re-use - only evaluating at new values
    fvals_new = func(mesh_new);
    
    % Recurrence
    int_old = int;
    int  = (1/2.0)*(int + h_old*(sum(fvals_new)));
    
    % Update Global Mesh
    mesh_temp = zeros(1,(2*n+1));
    mesh_temp(1:2:end) = mesh_current;
    mesh_temp(2:2:end) = mesh_new;
    mesh_current = mesh_temp;
    n = 2*n;
    
    % Store results
    integrals(i) = int;
    errors(i) = true_int - int;
    % Error estimate recurrence
    estimates(i-1) = (4/3)*(int - int_old);
end
end

