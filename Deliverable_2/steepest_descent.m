% Clear workspace
clc;              
clear;           
close all;        

% Define the symbolic function and its gradient
syms f(x, y);
f(x, y) = (x^5) * exp(-x^2 - y^2);  % Function definition
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
epsilon = 0.001;              % Convergence threshold for gradient norm
step = 0.05;                  % Constant step size for gradient descent

%% Gradient Descent with Constant Step Size

% Start from the initial point (0, 0)
f_point = [0 0];                            % Initial point
points_zero_const = zeros([1 2]);           % Store all points visited
points_zero_const(1, :) = f_point;          % Add initial point to the list
gradf_value = double(gradf(f_point(1), f_point(2)))';  % Evaluate gradient

while norm(gradf_value) >= epsilon
    f_point = f_point - step * gradf_value;             % Update using gradient descent
    points_zero_const(end + 1, :) = f_point;            % Append the new point
    gradf_value = double(gradf(f_point(1), f_point(2)))'; % Recompute gradient
end

% Start from the initial point (-1, 1)
f_point = [-1 1];
points_one_const = zeros([1 2]);
points_one_const(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    f_point = f_point - step * gradf_value;
    points_one_const(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

% Start from the initial point (1, -1)
f_point = [1 -1];
points_two_const = zeros([1 2]);
points_two_const(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    f_point = f_point - step * gradf_value;
    points_two_const(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

%% Gradient Descent with Minimizing Step Size
% Start from (0, 0)
f_point = [0 0];
points_zero_min = zeros([1 2]);
points_zero_min(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    % Define a 1D function g(x) for line search along the gradient direction
    syms g(x);
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2));
    dg = diff(g);  % Compute derivative of g(x)

    % Perform line search using binary search to find optimal step size
    upper = 1; lower = 0; l = 0.1;  % Bounds for binary search and precision
    while upper - lower >= l
        midpoint = (upper + lower) / 2;
        if dg(midpoint) > 0
            upper = midpoint;
        else
            lower = midpoint;
        end
    end
    step = upper; % Use the optimal step size from the binary search
    
    % Update the point using gradient descent
    f_point = f_point - step * gradf_value;
    points_zero_min(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

% Repeat the process for (-1, 1)
f_point = [-1 1];
points_one_min = zeros([1 2]);
points_one_min(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    syms g(x);
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2));
    dg = diff(g);

    upper = 1; lower = 0; l = 0.1;
    while upper - lower >= l
        midpoint = (upper + lower) / 2;
        if dg(midpoint) > 0
            upper = midpoint;
        else
            lower = midpoint;
        end
    end
    step = upper;

    f_point = f_point - step * gradf_value;
    points_one_min(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

% Repeat the process for (1, -1)
f_point = [1 -1];
points_two_min = zeros([1 2]);
points_two_min(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    syms g(x);
    g(x) = f(f_point(1) + x * gradf_value(1), f_point(2) + x * gradf_value(2));
    dg = diff(g);

    upper = 1; lower = 0; l = 0.1;
    while upper - lower >= l
        midpoint = (upper + lower) / 2;
        if dg(midpoint) > 0
            upper = midpoint;
        else
            lower = midpoint;
        end
    end
    step = upper;

    f_point = f_point - step * gradf_value;
    points_two_min(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end

%% Gradient Descent with Armijo Rule
% Initialize Armijo parameters
a = 0.001; b = 0.2; s = 1;

% Start from (0, 0)
f_point = [0 0];
points_zero_arm = zeros([1 2]);
points_zero_arm(1, :) = f_point;
gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    f_point_old = f_point;
    f_point = f_point - step * gradf_value;

    % Armijo rule adjustment for step size
    mk = 1;
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
          -a * b^mk * s * (step * norm(gradf_value)^2)
        mk = mk + 1;
    end
    step = s * b^mk;

    points_zero_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
end


% Start from (-1, 1)
f_point = [-1 1];
points_one_arm = zeros([1 2]);
points_one_arm(1, :) = f_point;
k = 1;
step = 1; 

gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    f_point_old = f_point;
    f_point = f_point - step * gradf_value;
    mk = 1;
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            - a * b^mk * s * (step * gradf_value * (gradf_value'))
        mk = mk + 1;
    end
    step = s * b^mk;
    points_one_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
    k = k + 1;
end

% Start from (1, -1)
f_point = [1 -1];
points_two_arm = zeros([1 2]);
points_two_arm(1, :) = f_point;
k = 1;
step = 1;

gradf_value = double(gradf(f_point(1), f_point(2)))';

while norm(gradf_value) >= epsilon
    f_point_old = f_point;
    f_point = f_point - step * gradf_value;
    mk = 1;
    while f(f_point_old(1), f_point_old(2)) - f(f_point(1), f_point(2)) < ...
            - a * b^mk * s * (step * gradf_value * (gradf_value'))
        mk = mk + 1;
    end
    step = s * b^mk;
    points_two_arm(end + 1, :) = f_point;
    gradf_value = double(gradf(f_point(1), f_point(2)))';
    k = k + 1;
end


%% Plot Results
% Generate plots for each descent method (constant step, minimizing step, Armijo)
% Each plot will show the path on the contour of the function and function values at each iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constant step

figure(); 
fcontour(f);
hold on;
title('contour values of f, constant step starting point = (0, 0)')
xlabel('x');
ylabel('y');
plot3(points_zero_const(:, 1), points_zero_const(:, 2), f(points_zero_const(:, 1), points_zero_const(:, 2)), "r*");

figure("Name", "f value at each iteration, constant step = 0.005, starting point = (0, 0)");
size_zr = size(points_zero_const);
plot(1:size_zr(1), f(points_zero_const(:, 1), points_zero_const(:, 2)));
title('f value at each iteration, constant step, starting point = (0, 0)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_const_step_zero_iter = size_zr;


figure(); 
fcontour(f);
hold on;
title('contour values of f, constant step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*");

figure("Name", "f value at each iteration, constant step = 0.005, starting point = (-1, 1)");
size_zr = size(points_one_const);
plot(1:size_zr(1), f(points_one_const(:, 1), points_one_const(:, 2)));
title('f value at each iteration, constant step, starting point = (-1, 1)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_const_step_one_iter = size_zr;
    

figure();
fcontour(f);
hold on;
title('contour values of f, constant step, starting point = (1, -1)')
xlabel('x');
ylabel('y');
plot3(points_two_const(:, 1), points_two_const(:, 2), f(points_two_const(:, 1), points_two_const(:, 2)), "r*");

figure("Name", "f value at each iteration, constant step = 0.005, starting point = (1, -1)");
size_zr = size(points_two_const);
plot(1:size_zr(1), f(points_two_const(:, 1), points_two_const(:, 2)));
title('f value at each iteration, constant step, starting point = (1, -1)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_const_step_two_iter = size_zr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimizing step

figure();
% Contour lines of the function
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (0, 0)')
xlabel('x');
ylabel('y');
plot3(points_zero_min(:, 1), points_zero_min(:, 2), f(points_zero_min(:, 1), points_zero_min(:, 2)), "r*");

figure("Name", "f value at each iteration, minimizing step, starting point = (0, 0)");
size_zr = size(points_zero_min);
plot(1:size_zr(1), f(points_zero_min(:, 1), points_zero_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (0, 0)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_min_step_zero_iter = size_zr;
  

figure();
% Contour lines of the function
fcontour(f);
hold on;
title('contour values of f, minimizing step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_min(:, 1), points_one_min(:, 2), f(points_one_min(:, 1), points_one_min(:, 2)), "r*");

figure("Name", "f value at each iteration, minimizing step, starting point = (-1, 1)");
size_zr = size(points_one_min);
plot(1:size_zr(1), f(points_one_min(:, 1), points_one_min(:, 2)));
title('f value at each iteration, minimizing step, starting point = (-1, 1)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_min_step_one_iter = size_zr;
  

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
hold on;

sd_min_step_two_iter = size_zr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Armijo's rule

figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (0, 0)')
xlabel('x');
ylabel('y');
plot3(points_zero_arm(:, 1), points_zero_arm(:, 2), f(points_zero_arm(:, 1), points_zero_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, armijo's rule step, starting point = (0, 0)");
size_zr = size(points_zero_arm);
plot(1:size_zr(1), f(points_zero_arm(:, 1), points_zero_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (0, 0)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_arm_step_zero_iter = size_zr;


figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (-1, 1)')
xlabel('x');
ylabel('y');
plot3(points_one_arm(:, 1), points_one_arm(:, 2), f(points_one_arm(:, 1), points_one_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, armijo's rule step, starting point = (-1, 1)");
size_zr = size(points_one_arm);
plot(1:size_zr(1), f(points_one_arm(:, 1), points_one_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (-1, 1)')
xlabel("iteration number");
ylabel("function value");
hold on;

sd_arm_step_one_iter = size_zr;


figure();
fcontour(f);
hold on;
title('contour values of f, armijo rule step, starting point = (1, -1)')
xlabel('x');
ylabel('y');
plot3(points_two_arm(:, 1), points_two_arm(:, 2), f(points_two_arm(:, 1), points_two_arm(:, 2)), "r*");

figure("Name", "f value at each iteration, armijo's rule step, starting point = (1, -1)");
size_zr = size(points_two_arm);
plot(1:size_zr(1), f(points_two_arm(:, 1), points_two_arm(:, 2)));
title('f value at each iteration, armijo rule step, starting point = (1, -1)')
xlabel("iteration number");
ylabel("Function value");
hold on;

sd_arm_step_two_iter = size_zr;
