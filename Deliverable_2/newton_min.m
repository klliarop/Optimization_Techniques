% Clear the workspace and initialize the environment
clc;                % Clear the command window
clear;              % Remove all variables from the workspace
close all;          % Close all open figures

% Define the symbolic function and compute its derivatives
syms f(x, y);                      % Define a symbolic function f with variables x and y
f(x, y) = (x^5) * exp(-x^2 - y^2); % Function definition: product of a polynomial and exponential decay
h = hessian(f);                    % Compute the Hessian matrix (second partial derivatives)
gradf = gradient(f);               % Compute the gradient (vector of first partial derivatives)

%% Visualization of the Function
% Create a 3D surface plot of the function for visualization
figure("Name", "Function plot");   % Open a new figure and name it
fsurf(f);                          % Plot the surface of the function in 3D
xlabel("x");                       % Label the x-axis
ylabel("y");                       % Label the y-axis
zlabel("f(x, y)");                 % Label the z-axis (function value)
hold on;                           % Allow additional plots on the same figure

%% Initialize Parameters
epsilon = 0.001;                   % Convergence threshold for the gradient norm
step = 0.05;                       % Initial constant step size (will be overwritten in the loop)

%% Gradient Descent: Starting from (0, 0)
% Initialize starting point and storage for visited points
f_point = [0 0];                   % Starting point (x0, y0)
points_zero_min = zeros([1 2]);    % Preallocate an array to store visited points
points_zero_min(1, :) = f_point;   % Store the initial point
k = 1;                             % Iteration counter

% Compute the gradient at the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent with minimization of step size
while norm(gradf_value) >= epsilon % Stop when the gradient norm is below the threshold
    % Define a symbolic function g(x) representing f along the search direction
    syms g(x); 
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2)); 
    dg = diff(g);                  % Compute the derivative of g(x)

    % Bisection method to find the step size that minimizes g(x)
    upper = 1;                     % Upper bound for the step size
    lower = 0;                     % Lower bound for the step size
    l = 0.1;                       % Desired precision of the bisection method

    while upper - lower >= l       % Continue until the interval width is less than precision
        midpoint = (upper + lower) / 2; % Compute the midpoint of the interval
        der_g = dg(midpoint);          % Evaluate the derivative of g(x) at the midpoint
        if der_g > 0                  % If derivative is positive, minimize in the left interval
            upper = midpoint;
        elseif der_g < 0              % If derivative is negative, minimize in the right interval
            lower = midpoint;
        else                          % If derivative is zero, we've found the minimum
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper;                    % Update the step size to the midpoint of the final interval

    % Update the point using the gradient descent formula
    f_point = f_point - step * gradf_value; 
    points_zero_min(end + 1, :) = f_point; % Store the updated point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Update the gradient
    k = k + 1;                        % Increment the iteration counter
end

%% Gradient Descent: Starting from (-1, 1)
f_point = [-1 1];                % Initial point (x0, y0)
points_one_min = zeros([1 2]);   % Preallocate array to store visited points
points_one_min(1, :) = f_point;  % Store the initial point
k = 1;                           % Iteration counter

gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    syms g(x); 
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2));
    dg = diff(g);

    upper = 1;
    lower = 0;
    l = 0.1;

    while upper - lower >= l
        midpoint = (upper + lower) / 2;
        der_g = dg(midpoint);
        if der_g > 0
            upper = midpoint;
        elseif der_g < 0
            lower = midpoint;
        else
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper;

    f_point = f_point - step * gradf_value;
    points_one_min(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
    k = k + 1;
end

%% Gradient Descent: Starting from (1, -1)
f_point = [1 -1];
points_two_min = zeros([1 2]);
points_two_min(1, :) = f_point;
k = 1;

gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    syms g(x);
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2));
    dg = diff(g);

    upper = 1;
    lower = 0;
    l = 0.1;

    while upper - lower >= l
        midpoint = (upper + lower) / 2;
        der_g = dg(midpoint);
        if der_g > 0
            upper = midpoint;
        elseif der_g < 0
            lower = midpoint;
        else
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper;

    f_point = f_point - step * gradf_value;
    points_two_min(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
    k = k + 1;
end

%% Results Visualization
% Plot results for each starting point

% Case 1: Starting from (0, 0)
figure();
fcontour(f); % Plot the contour lines of the function
hold on;
title('contour values of f, minimizing step, starting point = (0, 0)')
xlabel('x');
ylabel('y');
plot3(points_zero_min(:, 1), points_zero_min(:, 2), ...
      f(points_zero_min(:, 1), points_zero_min(:, 2)), "r*");

figure("Name", "f value at each iteration, minimizing step, starting point = (0, 0)");
size_zr = size(points_zero_min);
plot(1:size_zr(1), f(points_zero_min(:, 1), points_zero_min(:, 2)));title('f value at each iteration, minimizing step, starting point = (0, 0)')
xlabel("iteration number");
ylabel("function value");

% Case 2: Starting from (-1, 1)
figure();
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_min(:, 1), points_one_min(:, 2), ...
      f(points_one_min(:, 1), points_one_min(:, 2)), "r*");

figure("Name", "f value at each iteration, minimizing step, starting point = (-1, 1)");
size_zr = size(points_one_min);
plot(1:size_zr(1), f(points_one_min(:, 1), points_one_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (-1, 1)')
xlabel("iteration number");
ylabel("function value");

% Case 3: Starting from (1, -1)
figure();
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (1, -1)')
xlabel('x');
ylabel('y');
plot3(points_two_min(:, 1), points_two_min(:, 2), f(points_two_min(:, 1), points_two_min(:, 2)), "r*");

figure("Name", "f value at each iteration, minimizing step, starting point = (1, -1)");
size_zr = size(points_two_min);
plot(1:size_zr(1), f(points_two_min(:, 1), points_two_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (1, -1)')
xlabel("iteration number");
ylabel("function value");
