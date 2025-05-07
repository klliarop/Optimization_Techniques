% Clear the workspace and initialize the environment
clc; 
clear;   
close all;    

% Define the symbolic function and its derivatives
syms f(x, y);                    % Define a symbolic function f with variables x and y
f(x, y) = (x^5) * exp(-x^2 - y^2);  % Define the function: a product of polynomial and exponential decay
h = hessian(f);                  % Compute the Hessian matrix of f (second partial derivatives)
gradf = gradient(f);             % Compute the gradient of f (vector of first partial derivatives)

%% Visualization of the Function
% Create a 3D surface plot to visualize the function behavior
figure("Name", "Function plot"); % Open a new figure and set its title
fsurf(f);                        % Plot the 3D surface of the function
xlabel("x");                     % Label the x-axis
ylabel("y");                     % Label the y-axis
zlabel("f(x, y)");               % Label the z-axis (function value)
hold on;                         % Allow additional plots on the same figure

%% Initialization of Parameters
epsilon = 0.001;                 % Convergence threshold for the optimization
step = 0.05;                     % Initial constant step size 

% Initialize arrays to store optimization paths
f_point = [0 0];                 % Starting point of optimization
points_zero_min = zeros([1 2]);  % Array to store path of points for the first optimization
points_zero_min(1, :) = f_point; % Add the initial point
k = 1;                           % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

%% Main optimization loop for starting point (0, 0)
while norm(gradf_value) >= epsilon
    % Use 1D optimization to find the optimal step size
    upper = 1;                   % Initial upper bound for the step size
    lower = 0;                   % Initial lower bound
    l = 0.1;                     % Precision for the line search

    syms g(x);                   % Define a 1D symbolic function for the line search
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2)); % f along the gradient direction
    dg = diff(g);                % Compute the derivative of g with respect to x

    % Bisection method to find the optimal step size
    while upper - lower >= l
        midpoint = (upper + lower) / 2; % Calculate the midpoint
        der_g = dg(midpoint);           % Evaluate the derivative at the midpoint
        if der_g > 0
            upper = midpoint;           % Update the upper bound
        elseif der_g < 0
            lower = midpoint;           % Update the lower bound
        else
            % Derivative is zero; step size is optimal
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper; % Choose the final step size (upper == lower)

    % Ensure Hessian matrix is positive definite
    mk = 0; % Initialize regularization parameter
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2)))) <= 0
        mk = mk + 1; % Increase mk until the modified Hessian is positive definite
    end

    % Compute the descent direction using modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2))); 
    f_point = f_point - step * (dk * gradf_value')'; % Update the current point
    points_zero_min(end + 1, :) = f_point;          % Store the new point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient

    k = k + 1; % Increment iteration count
end

%% Repeat optimization for starting point (-1, 1)
% Initialize starting point and storage for the path
f_point = [-1 1];
points_one_min = zeros([1 2]);        % Array to store optimization path
points_one_min(1, :) = f_point;      % Add the initial point
k = 1;                               % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon
    % Use 1D optimization to find the optimal step size
    upper = 1;                       % Initial upper bound for the step size
    lower = 0;                       % Initial lower bound
    l = 0.1;                         % Precision for the line search

    syms g(x);                       % Define a 1D symbolic function for the line search
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2)); % f along the gradient direction
    dg = diff(g);                    % Compute the derivative of g with respect to x

    % Bisection method to find the optimal step size
    while upper - lower >= l
        midpoint = (upper + lower) / 2; % Calculate the midpoint
        der_g = dg(midpoint);           % Evaluate the derivative at the midpoint
        if der_g > 0
            upper = midpoint;           % Update the upper bound
        elseif der_g < 0
            lower = midpoint;           % Update the lower bound
        else
            % Derivative is zero; step size is optimal
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper; % Choose the final step size (upper == lower)

    % Ensure Hessian matrix is positive definite
    mk = 0; % Initialize regularization parameter
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2)))) <= 0
        mk = mk + 1; % Increase mk until the modified Hessian is positive definite
    end

    % Compute the descent direction using modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2))); 
    f_point = f_point - step * (dk * gradf_value')'; % Update the current point
    points_one_min(end + 1, :) = f_point;           % Store the new point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient

    k = k + 1; % Increment iteration count
end

%% Repeat optimization for starting point (1, -1)
% Initialize starting point and storage for the path
f_point = [1 -1];
points_two_min = zeros([1 2]);       % Array to store optimization path
points_two_min(1, :) = f_point;     % Add the initial point
k = 1;                              % Iteration counter

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon
    % Use 1D optimization to find the optimal step size
    upper = 1;                       % Initial upper bound for the step size
    lower = 0;                       % Initial lower bound
    l = 0.1;                         % Precision for the line search

    syms g(x);                       % Define a 1D symbolic function for the line search
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2)); % f along the gradient direction
    dg = diff(g);                    % Compute the derivative of g with respect to x

    % Bisection method to find the optimal step size
    while upper - lower >= l
        midpoint = (upper + lower) / 2; % Calculate the midpoint
        der_g = dg(midpoint);           % Evaluate the derivative at the midpoint
        if der_g > 0
            upper = midpoint;           % Update the upper bound
        elseif der_g < 0
            lower = midpoint;           % Update the lower bound
        else
            % Derivative is zero; step size is optimal
            upper = midpoint;
            lower = midpoint;
        end
    end

    step = upper; % Choose the final step size (upper == lower)

    % Ensure Hessian matrix is positive definite
    mk = 0; % Initialize regularization parameter
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2)))) <= 0
        mk = mk + 1; % Increase mk until the modified Hessian is positive definite
    end

    % Compute the descent direction using modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2))); 
    f_point = f_point - step * (dk * gradf_value')'; % Update the current point
    points_two_min(end + 1, :) = f_point;           % Store the new point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient

    k = k + 1; % Increment iteration count
end

%% Plot results for starting point (0, 0)
% Contour plot with optimization path
figure();
fcontour(f); % Plot contour lines of the function
hold on;
title('Contour plot and optimization path, starting point = (0, 0)');
xlabel('x');
ylabel('y');
plot3(points_zero_min(:, 1), points_zero_min(:, 2), f(points_zero_min(:, 1), points_zero_min(:, 2)), "r*"); % Plot optimization path
hold on;

% Plot function values at each iteration
figure("Name", "Function value vs iteration, starting point = (0, 0)");
size_zr = size(points_zero_min); % Determine the number of iterations
plot(1:size_zr(1), f(points_zero_min(:, 1), points_zero_min(:, 2))); % Plot function values
title('Function value at each iteration, starting point = (0, 0)');
xlabel("Iteration number");
ylabel("Function value");
hold on;
lm_min_step_zero_iter = size_zr; % Store the number of iterations

%% Plot results for starting point (-1, 1)

figure();
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_min(:, 1), points_one_min(:, 2), f(points_one_min(:, 1), points_one_min(:, 2)), "r*");
hold on;

figure("Name", "f value at each iteration, minimizing step, starting point = (-1, 1)");
size_zr = size(points_one_min);
plot(1:size_zr(1), f(points_one_min(:, 1), points_one_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (-1, 1)')
xlabel("iteration no");
ylabel("function value");
hold on;
lm_min_step_one_iter = size_zr;


figure("Name", "Contour plot with points and connecting lines, starting point = (-1, 1)");
fsurf(f, [-2, 2, -2, 2]); % Create a transparent surface plot for the function
xlabel("x");
ylabel("y");
zlabel("f(x, y)");
hold on;
title('Contour plot and optimization path, starting point = (-1, 1)');
plot3(points_one_min(:, 1), points_one_min(:, 2), f(points_one_min(:, 1), points_one_min(:, 2)), "r*", 'MarkerSize', 8);
plot3(points_one_min(:, 1), points_one_min(:, 2), f(points_one_min(:, 1), points_one_min(:, 2)), 'r-', 'LineWidth', 1.5);
grid on;
alpha 0.3; % Set surface transparency
hold on;

%% Plot results for starting point (1, -1)
figure();
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (1, -1)')
xlabel('x');
ylabel('y');
plot3(points_two_min(:, 1), points_two_min(:, 2), f(points_two_min(:, 1), points_two_min(:, 2)), "r*");
hold on;

figure("Name", "f value at each iteration, minimizing step, starting point = (1, -1)");
size_zr = size(points_two_min);
plot(1:size_zr(1), f(points_two_min(:, 1), points_two_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (1, -1)')
xlabel("iteration no");
ylabel("function value");
hold on;
lm_min_step_two_iter = size_zr;

