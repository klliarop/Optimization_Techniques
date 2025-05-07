clc;
clear;
close all;

syms f(x, y);
f(x, y) = (1/3) * x^2 + 3 * y^2;
gradf = gradient(f);

% plot function
figure("Name", "Plot");
fsurf(f);
xlabel("x");
ylabel("y");
zlabel("f(x, y)");
hold on;

% set limits
xlim = [-10 5];
ylim = [-8 12];

e = 0.01;
step = 0.2;
sk = 0.1;
start = [8 -10];

arr = zeros([1 2]);
arr(1, :) = start;

f_point = start;
k = 1;
gradf_val = double(gradf(f_point(1), f_point(2)))';
while (norm(gradf_val) >= e) 
    xbarint = f_point - sk * gradf_val;
    
    xold = xbarint(1);
    yold = xbarint(2);
    xnew = xold;
    if xold >= xlim(2)
        xnew = xlim(2);
    end
    if xold <= xlim(1)
        xnew = xlim(1);
    end
    ynew = yold;
    if yold >= ylim(2)
        ynew = ylim(2);
    end
    if yold <= ylim(1)
        ynew = ylim(1);
    end
    
    [xbar, ybar] = deal(xnew, ynew);

    f_point = f_point + step * ([xbar, ybar] - f_point);
    arr(end + 1, :) = f_point;
    gradf_val = double(gradf(f_point(1), f_point(2)))';
    k = k + 1;
end

%% plots
figure();
fcontour(f);
hold on;
title("Starting point (8, -10)");
xlabel("x");
ylabel("y");
plot3(arr(:, 1), arr(:, 2), f(arr(:, 1), arr(:, 2)), 'g-', 'LineWidth', 2); % Connected path
plot3(arr(:, 1), arr(:, 2), f(arr(:, 1), arr(:, 2)), "r*");
hold on;

figure();
size_arr = size(arr(:, :));
plot(1:size_arr(1), f(arr(:, 1), arr(:, 2)));
xlabel("iteration");
ylabel("f(x,y)");
hold on;
