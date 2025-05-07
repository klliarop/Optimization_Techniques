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

%% Case 1: Starting from (0, 0)
% Initialize the starting point and storage for visited points
f_point = [0 0];                 % Initial point (x0, y0)
points_zero_const = zeros([1 2]); % Array to store the points visited
points_zero_const(1, :) = f_point; % Store the initial point
k = 1;                           % Iteration counter

% Compute the gradient at the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent with constant step size
while norm(gradf_value) >= epsilon % Stop when the gradient norm is less than epsilon
    % Update the point using the gradient and Hessian
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value')';
    points_zero_const(end + 1, :) = f_point; % Store the updated point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Update the gradient
    k = k + 1;                       % Increment the iteration counter
end

%% Case 2: Starting from (-1, 1)
% Initialize the starting point and storage for visited points
f_point = [-1 1];                % Initial point (x0, y0)
points_one_const = zeros([1 2]); % Array to store the points visited
points_one_const(1, :) = f_point; % Store the initial point
k = 1;                           % Iteration counter

% Compute the gradient at the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent with constant step size
while norm(gradf_value) >= epsilon % Stop when the gradient norm is less than epsilon
    % Update the point using the gradient and Hessian
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value')';
    points_one_const(end + 1, :) = f_point; % Store the updated point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Update the gradient
    k = k + 1;                       % Increment the iteration counter
end

%% Case 3: Starting from (1, -1)
% Initialize the starting point and storage for visited points
f_point = [1 -1];                % Initial point (x0, y0)
points_two_const = zeros([1 2]); % Array to store the points visited
points_two_const(1, :) = f_point; % Store the initial point
k = 1;                           % Iteration counter

% Compute the gradient at the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent with constant step size
while norm(gradf_value) >= epsilon % Stop when the gradient norm is less than epsilon
    % Update the point using the gradient and Hessian
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value')';
    points_two_const(end + 1, :) = f_point; % Store the updated point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Update the gradient
    k = k + 1;                       % Increment the iteration counter
end

%% Results Visualization
% Case 1: Plot the trajectory for starting point (0, 0)
figure();                        % Create a new figure
fcontour(f);                     % Plot contour lines of the function
hold on;                         % Allow overlaying additional plots
title('contour values of f, constant step starting point = (0, 0)')
xlabel('x');                     % Label the x-axis
ylabel('y');                     % Label the y-axis
plot3(points_zero_const(:, 1), points_zero_const(:, 2),f(points_zero_const(:, 1), points_zero_const(:, 2)), "r*");

% Plot the function value at each iteration
figure("Name", "f value at each iteration, constant step = 0.05, starting point = (0, 0)");
size_zr = size(points_zero_const); % Get the number of iterations
plot(1:size_zr(1), f(points_zero_const(:, 1), points_zero_const(:, 2))); % Plot function values
title('f value at each iteration, constant step, starting point = (0, 0)')
xlabel("iteration number");        % Label the x-axis
ylabel("function value");          % Label the y-axis

% Case 2: Plot the trajectory for starting point (-1, 1)
figure();                        % Create a new figure
fcontour(f);                     % Plot contour lines of the function
hold on;                         % Allow overlaying additional plots
title('contour values of f, constant step starting point = (-1, 1)')
xlabel('x');                     % Label the x-axis
ylabel('y');                     % Label the y-axis
plot3(points_one_const(:, 1), points_one_const(:, 2), ...
      f(points_one_const(:, 1), points_one_const(:, 2)), "r*");

% Plot the function value at each iteration
figure("Name", "f value at each iteration, constant step = 0.05, starting point = (-1, 1)");
size_zr = size(points_one_const); % Get the number of iterations
plot(1:size_zr(1), f(points_one_const(:, 1), points_one_const(:, 2))); % Plot function values
title('f value at each iteration, constant step, starting point = (-1, 1)')
xlabel("iteration number");        % Label the x-axis
ylabel("function value");          % Label the y-axis

% Case 3: Plot the trajectory for starting point (1, -1)
figure();                        % Create a new figure
fcontour(f);                     % Plot contour lines of the function
hold on;                         % Allow overlaying additional plots
title('contour values of f, constant step starting point = (1, -1)')
xlabel('x');                     % Label the x-axis
ylabel('y');                     % Label the y-axis
plot3(points_two_const(:, 1), points_two_const(:, 2), f(points_two_const(:, 1), points_two_const(:, 2)), "r*");

% Plot the function value at each iteration
figure("Name", "f value at each iteration, constant step = 0.05, starting point = (1, -1)");
size_zr = size(points_two_const); % Get the number of iterations
plot(1:size_zr(1), f(points_two_const(:, 1), points_two_const(:, 2))); % Plot function values
title('f value at each iteration, constant step, starting point = (1, -1)')
xlabel("iteration number");        % Label the x-axis
ylabel("function value");          % Label the y-axis
