% Clear the workspace and initialize the environment
clc;               % Clear the command window
clear;             % Remove all variables from the workspace
close all;         % Close all open figures

% Define the symbolic function and compute derivatives
syms f(x, y);                          % Define symbolic function f with variables x and y
f(x, y) = (x^5) * exp(-x^2 - y^2);     % Define the function: polynomial with exponential decay
h = hessian(f);                        % Compute the Hessian matrix (second partial derivatives)
gradf = gradient(f);                   % Compute the gradient (vector of first partial derivatives)

%% Function Visualization
% Create a 3D surface plot of the function to understand its shape
figure("Name", "Function plot");       % Open a new figure and name it
fsurf(f);                              % Plot the surface of the function
xlabel("x");                           % Label the x-axis
ylabel("y");                           % Label the y-axis
zlabel("f(x, y)");                     % Label the z-axis (function value)
hold on;                               % Allow additional plots on the same figure

%% Initialize Parameters for Armijo Rule
epsilon = 0.001;                       % Convergence threshold for the gradient norm
a = 0.001;                             % Parameter controlling the decrease in function value
b = 0.2;                               % Scaling factor for step size reduction
s = 1;                                 % Initial step size for Armijo's rule
step = 1;                              % Initial step size for gradient descent


%% Gradient Descent with Armijo's Rule: Starting Point (0, 0)
% Initialize starting point and array to store visited points
f_point = [0 0];                       % Starting point (x0, y0)
points_zero_arm = zeros([1, 2]);       % Preallocate array to store visited points
points_zero_arm(1, :) = f_point;       % Store the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))'; % Compute the initial gradient

% Perform gradient descent using Armijo's rule
while norm(gradf_value) >= epsilon     % Stop when the gradient norm is below the threshold
    f_point_old = f_point;             % Store the current point for comparison
    % Update the point using Newton's method with Armijo step size
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2))) ) * gradf_value')';
    mk = 1;                            % Counter for backtracking step size
    
    % Armijo's condition to adjust the step size
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value'))
        mk = mk + 1;                   % Increase backtracking counter
    end
    step = s * b^mk;                   % Update step size using backtracking

    % Store the new point and update the gradient
    points_zero_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

%% Gradient Descent with Armijo's Rule: Starting Point (-1, 1)
% Initialize starting point and array to store visited points
f_point = [-1 1];                      % Starting point (x0, y0)
points_one_arm = zeros([1, 2]);        % Preallocate array to store visited points
points_one_arm(1, :) = f_point;        % Store the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent using Armijo's rule
while norm(gradf_value) >= epsilon
    f_point_old = f_point;             % Store the current point for comparison
    % Update the point using Newton's method with Armijo step size
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2))) ) * gradf_value')';
    mk = 1;                            % Counter for backtracking step size
    
    % Armijo's condition to adjust the step size
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value'))
        mk = mk + 1;                   % Increase backtracking counter
    end
    step = s * b^mk;                   % Update step size using backtracking

    % Store the new point and update the gradient
    points_one_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

%% Gradient Descent with Armijo's Rule: Starting Point (1, -1)
% Initialize starting point and array to store visited points
f_point = [1 -1];                      % Starting point (x0, y0)
points_two_arm = zeros([1, 2]);        % Preallocate array to store visited points
points_two_arm(1, :) = f_point;        % Store the initial point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Perform gradient descent using Armijo's rule
while norm(gradf_value) >= epsilon
    f_point_old = f_point;             % Store the current point for comparison
    % Update the point using Newton's method with Armijo step size
    f_point = f_point - step * (double(inv(h(f_point(1), f_point(2))) ) * gradf_value')';
    mk = 1;                            % Counter for backtracking step size
    
    % Armijo's condition to adjust the step size
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * (double(inv(h(f_point(1), f_point(2)))) * gradf_value'))
        mk = mk + 1;                   % Increase backtracking counter
    end
    step = s * b^mk;                   % Update step size using backtracking

    % Store the new point and update the gradient
    points_two_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end


%% Visualization of Results
% Case 1: Starting point (0, 0)
figure();
fcontour(f); % Plot the contour lines of the function
hold on;
xlabel('x');
ylabel('y');
plot3(points_zero_arm(:, 1), points_zero_arm(:, 2), ...
      f(points_zero_arm(:, 1), points_zero_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, Armijo's rule, starting point = (0, 0)");
size_zr = size(points_zero_arm);
plot(1:size_zr(1), f(points_zero_arm(:, 1), points_zero_arm(:, 2)));
xlabel("iteration number");
ylabel("function value");

% Case 2: Starting point (-1, 1)
figure();
fcontour(f);
hold on;
xlabel('x');
ylabel('y');
plot3(points_one_arm(:, 1), points_one_arm(:, 2), ...
      f(points_one_arm(:, 1), points_one_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, Armijo's rule, starting point = (-1, 1)");
size_zr = size(points_one_arm);
plot(1:size_zr(1), f(points_one_arm(:, 1), points_one_arm(:, 2)));
xlabel("iteration number");
ylabel("function value");

% Case 3: Starting point (1, -1)
figure();
fcontour(f);
hold on;
xlabel('x');
ylabel('y');
plot3(points_two_arm(:, 1), points_two_arm(:, 2), ...
      f(points_two_arm(:, 1), points_two_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, Armijo's rule, starting point = (1, -1)");
size_zr = size(points_two_arm);
plot(1:size_zr(1), f(points_two_arm(:, 1), points_two_arm(:, 2)));
xlabel("iteration number");
ylabel("function value");
