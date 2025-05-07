% Clear the workspace and initialize the environment
clc; 
clear;   
close all;    

% Define the symbolic function and its derivatives
syms f(x, y);                    % Define a symbolic function f with variables x and y
f(x, y) = (x^5) * exp(-x^2 - y^2);  % Function definition: a product of a polynomial and exponential decay
h = hessian(f);                  % Compute the Hessian matrix 
gradf = gradient(f);             % Compute the gradient (vector of first partial derivatives)

%% Visualization of the Function
% Create a 3D surface plot to visualize the behavior of the function
figure("Name", "Function plot"); % Open a new figure and set a title
fsurf(f);                        % Plot the surface of the function in 3D
xlabel("x");                     % Label the x-axis
ylabel("y");                     % Label the y-axis
zlabel("f(x, y)");               % Label the z-axis (function value)
hold on;                         % Allow additional plots on the same figure

%% Initialization of Parameters
epsilon = 0.001;                 % Convergence threshold 
step = 0.05;                     % Constant step size 

% Optimization using a constant step size
% Starting from different initial points (0, 0), (-1, 1), and (1, -1)
%% Start from (0, 0)
f_point = [0 0];                   % Initial point (0, 0)
points_zero_const = zeros([1 2]);  % Initialize array to store the optimization path
points_zero_const(1, :) = f_point; % Store the initial point
k = 1;                             % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the constant step size
    points_zero_const(end + 1, :) = f_point;        % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient at the new point
    k = k + 1; % Increment iteration counter
end

%% Start from (-1, 1)
f_point = [-1 1];                  % Initial point (-1, 1)
points_one_const = zeros([1 2]);   % Initialize array to store the optimization path
points_one_const(1, :) = f_point;  % Store the initial point
k = 1;                             % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the constant step size
    points_one_const(end + 1, :) = f_point;         % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient at the new point
    k = k + 1; % Increment iteration counter
end

%% Start from (1, -1)
f_point = [1 -1];                  % Initial point (1, -1)
points_two_const = zeros([1 2]);   % Initialize array to store the optimization path
points_two_const(1, :) = f_point;  % Store the initial point
k = 1;                             % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the constant step size
    points_two_const(end + 1, :) = f_point;         % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient at the new point
    k = k + 1; % Increment iteration counter
end


%% Plot results for starting point (0, 0)
% Contour plot with optimization path
figure();
fcontour(f);
hold on;
title('Contour plot and optimization path, starting point = (0, 0)');
xlabel('x');
ylabel('y');
plot3(points_zero_const(:, 1), points_zero_const(:, 2), f(points_zero_const(:, 1), points_zero_const(:, 2)), "r*");
hold on;

% Plot function values at each iteration
figure("Name", "f value at each iteration, constant step = 0.005, starting point = (0, 0)");
size_zr = size(points_zero_const);
plot(1:size_zr(1), f(points_zero_const(:, 1), points_zero_const(:, 2)));
title('Function value at each iteration, starting point = (0, 0)');
xlabel("iteration no");
ylabel("function value");
hold on;
lm_const_step_zero_iter = size_zr;


%% Plot results for starting point (-1, 1)
figure();
fcontour(f);
hold on;
title('Contour plot and optimization path, starting point = (-1, 1)');
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*");
hold on;

% Plot function values at each iteration
figure("Name", "f value at each iteration, constant step = 0.005, starting point = (-1, 1)");
size_zr = size(points_one_const);
plot(1:size_zr(1), f(points_one_const(:, 1), points_one_const(:, 2)));
title('Function value at each iteration, starting point = (-1, 1)');
xlabel("iteration no");
ylabel("function value");
hold on;
lm_const_step_one_iter = size_zr;

figure("Name", "Contour plot with points and connecting lines, starting point = (-1, 1)");
fsurf(f, [-2, 2, -2, 2]); % Create a transparent surface plot for the function
xlabel("x");
ylabel("y");
zlabel("f(x, y)");
hold on;
title('Contour plot and optimization path, starting point = (-1, 1)');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*", 'MarkerSize', 8);
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'r-', 'LineWidth', 1.5);
grid on;
alpha 0.3; % Set surface transparency
hold on;


%% Plot results for starting point (1, -1)
figure();
fcontour(f);
hold on;
title('Contour plot and optimization path, starting point = (1, -1)');
xlabel('x');
ylabel('y');
plot3(points_two_const(:, 1), points_two_const(:, 2), f(points_two_const(:, 1), points_two_const(:, 2)), "r*");
hold on;

% Plot function values at each iteration
figure("Name", "f value at each iteration, constant step = 0.005, starting point = (1, -1)");
size_zr = size(points_two_const);
plot(1:size_zr(1), f(points_two_const(:, 1), points_two_const(:, 2)));
title('Function value at each iteration, starting point = (1, -1)');
xlabel("iteration no");
ylabel("function value");
hold on;
lm_const_step_two_iter = size_zr;

