function [integrals, errors, estimates] = midadaptive(func,mesh,H,n,true_int, ref_num)
%MIDADAPTIVE Summary of this function goes here
%   Detailed explanation goes here

% Initializing results
integrals = zeros(1,ref_num+1);
errors = zeros(1,ref_num+1);
estimates = zeros(1,ref_num);

% Starting recurrence
[int_start, error_start] = midpoint(func,mesh,H,n,true_int);
integrals(1) = int_start;
errors(1) = error_start;

int = int_start;
h_shift = H;
mesh_current = mesh;

for i = 2:ref_num+1
    h_old = h_shift;
    
    % New mesh values
    h_shift = h_shift/3.0;
    mesh_left = mesh_current(1:n) + h_shift;
    mesh_right = mesh_current(2:n+1) - h_shift;
    
    % New values for evaluation (since h = H/2)
    h_local = h_shift/2.0;
    eval_left = mesh_left - h_local;
    eval_right = mesh_right + h_local;
    
    % Complete re-use - only evaluating at new values
    fvals_new_left = func(eval_left);
    fvals_new_right = func(eval_right); 
    
    % Recurrence
    int_old = int;
    int = (1/3.0)*(int + h_old*(sum(fvals_new_left) + sum(fvals_new_right)));
    
    % Update Global Mesh
    mesh_temp = zeros(1, 3*n + 1);
    mesh_temp(1:3:end) = mesh_current;
    mesh_temp(2:3:end) = mesh_left;
    mesh_temp(3:3:end) = mesh_right;
    mesh_current = mesh_temp;
    % I_{3n)
    n = 3*n;
    
    % Store results
    integrals(i) = int;
    errors(i) = true_int - int;
    % Error estimate recurrence
    estimates(i-1) = (9/8)*(int - int_old);
end
end