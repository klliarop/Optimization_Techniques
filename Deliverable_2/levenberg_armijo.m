% Clear workspace
clc;              
clear;           
close all;        

% Define the symbolic function and its gradient
syms f(x, y);
f(x, y) = (x^5) * exp(-x^2 - y^2);  % Function definition
h = hessian(f);                     % Compute Hessian
gradf = gradient(f);                % Compute the gradient of the function

%% Function Visualization
% Create a 3D surface plot of the function
figure("Name", "Function plot");
fsurf(f);                     % Plot the surface of the function
xlabel("x");                  % Label for x-axis
ylabel("y");                  % Label for y-axis
zlabel("f(x, y)");            % Label for z-axis
hold on;                      % Allow adding additional plots to the figure

%% Initialize Parameters
epsilon = 0.001;              % Convergence threshold 
step = 0.05;                  % Constant step size 

% Optimization using Armijo rule for step size adjustment
% Starting from different initial points (0, 0), (-1, -1), and (1, 1)

%% Arbitrary constants for Armijo rule
a = 0.001; % Parameter to control sufficient decrease condition
b = 0.2;   % Reduction factor for step size
s = 1;     % Initial step size scaling factor

%% Start from (0, 0)
f_point = [0 0];                  % Initial point (0, 0)
points_zero_arm = zeros([1 2]);   % Initialize array to store the optimization path
points_zero_arm(1, :) = f_point;  % Store the initial point
k = 1;                            % Iteration counter
step = 1;                         % Initial step size

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    f_point_old = f_point; % Store the current point for step size calculation
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the current step size
    
    % Recompute step size using the Armijo rule
    mk = 1; % Initialize step size reduction counter
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * dk * gradf_value')
        mk = mk + 1; % Increment mk to reduce step size
    end
    step = s * b^mk; % Update step size

    points_zero_arm(end + 1, :) = f_point; % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient
    k = k + 1; % Increment iteration counter
end

%% Start from (-1, -1)
f_point = [-1 -1];                % Initial point (-1, -1)
points_one_arm = zeros([1 2]);    % Initialize array to store the optimization path
points_one_arm(1, :) = f_point;   % Store the initial point
k = 1;                            % Iteration counter
step = 1;                         % Initial step size

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    f_point_old = f_point; % Store the current point for step size calculation
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the current step size
    
    % Recompute step size using the Armijo rule
    mk = 1; % Initialize step size reduction counter
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * dk * gradf_value')
        mk = mk + 1; % Increment mk to reduce step size
    end
    step = s * b^mk; % Update step size

    points_one_arm(end + 1, :) = f_point; % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient
    k = k + 1; % Increment iteration counter
end

%% Start from (1, 1)
f_point = [1 1];                  % Initial point (1, 1)
points_two_arm = zeros([1 2]);    % Initialize array to store the optimization path
points_two_arm(1, :) = f_point;   % Store the initial point
k = 1;                            % Iteration counter
step = 1;                         % Initial step size

% Compute the initial gradient at the starting point
gradf_value = double(gradf(f_point(1), f_point(2)))';

% Main optimization loop
while norm(gradf_value) >= epsilon % Continue until the gradient norm is below the threshold
    f_point_old = f_point; % Store the current point for step size calculation
    mk = 0; % Regularization parameter to ensure the Hessian is positive definite
    
    % Adjust mk until the modified Hessian is positive definite
    while min(double(eig(h(f_point(1), f_point(2)) + mk * eye(2))) > 0) == 0
        mk = mk + 1; % Increment mk to ensure positive definiteness
    end

    % Compute the descent direction using the modified Hessian
    dk = double(inv(h(f_point(1), f_point(2)) + mk * eye(2)));
    f_point = f_point - step * (dk * gradf_value')'; % Update the point using the current step size
    
    % Recompute step size using the Armijo rule
    mk = 1; % Initialize step size reduction counter
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            -a * b^mk * s * (step * dk * gradf_value')
        mk = mk + 1; % Increment mk to reduce step size
    end
    step = s * b^mk; % Update step size

    points_two_arm(end + 1, :) = f_point; % Store the new point in the path
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute the gradient
    k = k + 1; % Increment iteration counter
end


%% Plot results for starting point (0, 0)
figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (0, 0)')
xlabel('x');
ylabel('y');
plot3(points_zero_arm(:, 1), points_zero_arm(:, 2), f(points_zero_arm(:, 1), points_zero_arm(:, 2)), "r*");

hold on;
figure("Name", "f value at each iteration, armijo's rule step, starting point = (0, 0)");
size_zr = size(points_zero_arm);
plot(1:size_zr(1), f(points_zero_arm(:, 1), points_zero_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (0, 0)')
xlabel("iteration number");
ylabel("function value");
hold on;
lm_arm_step_zero_iter = size_zr;

%% Plot results for starting point (-1, 1)

figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_arm(:, 1), points_one_arm(:, 2), f(points_one_arm(:, 1), points_one_arm(:, 2)), "r*");

hold on;
figure("Name", "f value at each iteration, armijo's rule step, starting point = (-1, 1)");
size_zr = size(points_one_arm);
plot(1:size_zr(1), f(points_one_arm(:, 1), points_one_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (-1, 1)')
xlabel("iteration number");
ylabel("function value");
hold on;
lm_arm_step_one_iter = size_zr;

figure("Name", "Contour plot with points and connecting lines, starting point = (-1, 1)");
fsurf(f, [-2, 2, -2, 2]); % Create a transparent surface plot for the function
xlabel("x");
ylabel("y");
zlabel("f(x, y)");
hold on;
title('Contour plot and optimization path, starting point = (-1, 1)');
plot3(points_one_arm(:, 1), points_one_arm(:, 2), f(points_one_arm(:, 1), points_one_arm(:, 2)), "r*", 'MarkerSize', 8);
plot3(points_one_arm(:, 1), points_one_arm(:, 2), f(points_one_arm(:, 1), points_one_arm(:, 2)), 'r-', 'LineWidth', 1.5);
grid on;
alpha 0.3; % Set surface transparency
hold on;

%% Plot results for starting point (1, -1)

figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (1, -1)')
xlabel('x');
ylabel('y');
plot3(points_two_arm(:, 1), points_two_arm(:, 2), f(points_two_arm(:, 1), points_two_arm(:, 2)), "r*");
hold on;

figure("Name", "f value at each iteration, armijo's rule step, starting point = (1, -1)");
size_zr = size(points_two_arm);
plot(1:size_zr(1), f(points_two_arm(:, 1), points_two_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (1, -1)')
xlabel("iteration number");
ylabel("function value");
hold on;
lm_arm_step_two_iter = size_zr;
