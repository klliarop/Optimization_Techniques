clc;
clear;
close all;

syms f(x, y);
f(x, y) = (1/3) * x^2 + 3 * y^2;
gradf = gradient(f);

% plot
figure("Name", "Function plot");
fsurf(f);
xlabel("x");
ylabel("y");
zlabel("f(x, y)");
hold on;

e = 0.001;

%% γk = 0.1
step = 0.1;

f_point = [-5 -6];   %starting point
points_one_const = zeros([1 2]);
points_one_const(1, :) = f_point;
k = 1;
gradf_val = double(gradf(f_point(1), f_point(2)))';
while norm(gradf_val) >= e
    f_point = f_point - step * gradf_val;
    points_one_const(end + 1, :) = f_point;
    gradf_val = double(gradf(f_point(1), f_point(2)))';
end

% Plotting the contour and path of the gradient descent
figure("Name", "γ = 0.1"); 
fcontour(f);
hold on;
title("γ = 0.1");
title("γκ = 0.1")
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'g-', 'LineWidth', 2); % Connected path
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'r*'); % Mark points

%% γk = 0.3
step = 0.3;

f_point =  [-5 -6];    %starting point
points_one_const = zeros([1 2]);
points_one_const(1, :) = f_point;
k = 1;
gradf_val = double(gradf(f_point(1), f_point(2)))';
while norm(gradf_val) >= e
    f_point = f_point - step * gradf_val;
    points_one_const(end + 1, :) = f_point;
    gradf_val = double(gradf(f_point(1), f_point(2)))';
end

% Contour lines of the function
figure("Name", "γ = 0.3"); 
fcontour(f);
hold on;
title("γ = 0.3");
title("γκ = 0.3")
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'g-', 'LineWidth', 2); % Connected path
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*");

%% γk = 3
step = 3;

f_point =  [-5 -6];    %starting point
points_one_const = zeros([1 2]);
points_one_const(1, :) = f_point;
k = 1;
gradf_val = double(gradf(f_point(1), f_point(2)))';
while norm(gradf_val) >= e
    f_point = f_point - step * gradf_val;
    points_one_const(end + 1, :) = f_point;
    gradf_val = double(gradf(f_point(1), f_point(2)))';
end

% Contour lines of the function
figure("Name", "γ = 3"); 
fcontour(f);
hold on;
title("γκ = 3")
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'g-', 'LineWidth', 2); % Connected path
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*");

%% γk = 5
step = 5;

f_point =  [-5 -6];    %starting point
points_one_const = zeros([1 2]);
points_one_const(1, :) = f_point;
k = 1;
gradf_val = double(gradf(f_point(1), f_point(2)))';
while norm(gradf_val) >= e
    f_point = f_point - step * gradf_val;
    points_one_const(end + 1, :) = f_point;
    gradf_val = double(gradf(f_point(1), f_point(2)))';
end

% Contour lines of the function
figure(); 
fcontour(f);
hold on;
title("γκ = 5")
xlabel('x');
ylabel('y');
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), 'g-', 'LineWidth', 2); % Connected path
plot3(points_one_const(:, 1), points_one_const(:, 2), f(points_one_const(:, 1), points_one_const(:, 2)), "r*");
