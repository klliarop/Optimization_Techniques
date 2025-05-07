% Clear workspace, close all figures, and clear command window
close all;
clear all;
clc;

% Define symbolic functions f1(x), f2(x), and f3(x)
syms f1(x);
syms f2(x);
syms f3(x);
f1(x) = (x - 2)^2 + x*log(x + 3);
f2(x) = exp(-2*x) + (x - 2)^2;
f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x);

% Define bounds for x
lower = -1;
upper = 3;

% Plot of function f1(x)
figure("Name", "f1(x)")
fplot(f1, [lower, upper]); % Plot f1 in the range [lower, upper]
title('f1(x) = (x - 2)^2 + x*log(x + 3)');
xlabel("x");
ylabel("y");
hold on;

% Plot of function f2(x)
figure("Name", "f2(x)")
fplot(f2, [lower, upper]); % Plot f2 in the range [lower, upper]
title('f2(x) = exp(-2*x) + (x - 2)^2');
xlabel("x");
ylabel("y");
hold on;

% Plot of function f3(x)
figure("Name", "f3(x)")
fplot(f3, [lower, upper]); % Plot f3 in the range [lower, upper]
title('f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x)');
xlabel("x");
ylabel("y");

% The following line adjusts the y-axis for better visibility of details
%ylim([-5, 15]);
yline(0, 'k'); % Add horizontal line at y=0
hold on;

% Set up constant l and varying e parameters for bisection method
l = 0.01; % constant tolerance l for stopping criterion
n_steps = 100; % number of steps for varying e
e_first = 0.00001; % initial value of e
e_last = (l / 2) - e_first; % final value of e for linspace

E = linspace(e_first, e_last, n_steps); % Create array of e values

A = zeros([n_steps 3]); % Initialize matrix to store iteration counts for each function

% Perform bisection method for each value of e
for i = 1:n_steps
    e = E(i); % Get current e value from the array

    % Initialize counters for each function
    reps_1 = 1;
    reps_2 = 1;
    reps_3 = 1;

    % Apply bisection method to f1(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        % Update interval based on function evaluations
        if f1(x_1k) < f1(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_1 = reps_1 + 1; % Increment counter for f1
    end

    % Apply bisection method to f2(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        % Update interval based on function evaluations
        if f2(x_1k) < f2(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_2 = reps_2 + 1; % Increment counter for f2
    end

    % Apply bisection method to f3(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        % Update interval based on function evaluations
        if f3(x_1k) < f3(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_3 = reps_3 + 1; % Increment counter for f3
    end

    % Store iteration counts for current e in matrix A
    A(i, :) = [reps_1 reps_2 reps_3];
end

% Plot the variation of function evaluations with respect to e
figure("Name", " για l = 0.01 και μεταβαλλόμενο e")
plot(E, (A(:, 1) - 1) * 2, E, (A(:, 2) - 1) * 2, E, (A(:, 3) - 1) * 2);
title('Μεταβολή υπολογισμών αντικειμενικής συνάρτησης ');
subtitle('l = 0.01 και μεταβαλλόμενο e');
xlabel("e");
ylabel("Μεταβολή των υπολογισμών της f(x)");
legend("f1(x)", "f2(x)", "f3(x)");
hold on;

% Set up constant e and varying l parameters
l_steps = 100; % number of steps for varying l
l_first = 0.005; % initial value of l
l_last = 0.1; % final value of l
e = 0.001; % constant tolerance e for this part
ls = linspace(l_first, l_last, l_steps); % Create array of l values

B = zeros([l_steps 3]); % Initialize matrix to store iteration counts for each function
low_lim = zeros([l_steps 1 3]);
high_lim = zeros([l_steps 1 3]);

% Perform bisection method for each value of l
for i = 1:l_steps
    l = ls(i); % Get current l value from the array

    % Initialize counters for each function
    reps_1 = 1;
    reps_2 = 1;
    reps_3 = 1;

    % Set initial boundaries for visualization
    low_lim(i, 1, :) = [0 0 0];
    high_lim(i, 1, :) = [3 3 3];

    % Apply bisection method to f1(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        % Update interval based on function evaluations
        if f1(x_1k) < f1(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_1 = reps_1 + 1;

        % Record boundaries for each iteration
        low_lim(i, reps_1, 1) = a;
        high_lim(i, reps_1, 1) = b;
    end

    % Repeat the process for f2(x) and f3(x), storing boundaries
    % Apply bisection method to f2(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        if f2(x_1k) < f2(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_2 = reps_2 + 1;
        low_lim(i, reps_2, 2) = a;
        high_lim(i, reps_2, 2) = b;
    end

    % Apply bisection method to f3(x)
    a = lower;
    b = upper;
    while b - a >= l
        x_1k = (a + b) / 2 - e;
        x_2k = (a + b) / 2 + e;

        if f3(x_1k) < f3(x_2k)
            b = x_2k;
        else
            a = x_1k;
        end

        reps_3 = reps_3 + 1;
        low_lim(i, reps_3, 3) = a;
        high_lim(i, reps_3, 3) = b;
    end

    % Store iteration counts for current l in matrix B
    B(i, :) = [reps_1 reps_2 reps_3];
end

% Plot the variation of function evaluations with respect to l
figure("Name", "e = 0.001 και μεταβαλλόμενο l")
plot(ls, (B(:, 1) - 1) * 2, ls, (B(:, 2) - 1) * 2 , ls, (B(:, 3) - 1) * 2);
title('Μεταβολή υπολογισμών αντικειμενικής συνάρτησης ');
subtitle('e = 0.001 και μεταβαλλόμενο l');
xlabel("l");
ylabel("Μεταβολή των υπολογισμών της f(x)");
legend("f1(x)", "f2(x)", "f3(x)");
hold on;

% Visualization of interval reduction process for each function at varying l values

% Plot upper and lower limits for f1(x) over iterations with different l values
var = 1;
figure("Name", "f1(x)")
plot(1:B(var, 1), low_lim(var, 1:B(var, 1), 1), ...
    1:B(var, 1), high_lim(var, 1:B(var, 1), 1));
hold on;
var =  l_steps / 2;
plot(1:B(var, 1), low_lim(var, 1:B(var, 1), 1), ...
    1:B(var, 1), high_lim(var, 1:B(var, 1), 1));
hold on;
var =  l_steps;
plot(1:B(var, 1), low_lim(var, 1:B(var, 1), 1), ...
    1:B(var, 1), high_lim(var, 1:B(var, 1), 1));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f1) l = 0.005", "upper(f1) l = 0.005", "lower(f1) l = 0.01 " , "upper(f1) l = 0.01", "lower(f1) l = 0.1 ", "upper(f1) l = 0.1" );

% Plot upper and lower limits for f2(x) over iterations with different l values
var = 1;
figure("Name", "f2(x)")
plot(1:B(var, 2), low_lim(var, 1:B(var, 2), 2), ...
    1:B(var, 2), high_lim(var, 1:B(var, 2), 2));
hold on;
var =  l_steps / 2;
plot(1:B(var, 2), low_lim(var, 1:B(var, 2), 2), ...
    1:B(var, 2), high_lim(var, 1:B(var, 2), 2));
hold on;
var =  l_steps;
plot(1:B(var, 2), low_lim(var, 1:B(var, 2), 2), ...
    1:B(var, 2), high_lim(var, 1:B(var, 2), 2));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f2) l = 0.005", "upper(f2) l = 0.005", "lower(f2) l = 0.01 " , "upper(f2) l = 0.01", "lower(f2) l = 0.1 ", "upper(f2) l = 0.1" );

% Plot upper and lower limits for f3(x) over iterations with different l values
var = 1;
figure("Name", "f3(x)")
plot(1:B(var, 3), low_lim(var, 1:B(var, 3), 3), ...
    1:B(var, 3), high_lim(var, 1:B(var, 3), 3));
hold on;
var =  l_steps / 2;
plot(1:B(var, 3), low_lim(var, 1:B(var, 3), 3), ...
    1:B(var, 3), high_lim(var, 1:B(var, 3), 3));
hold on;
var =  l_steps;
plot(1:B(var, 3), low_lim(var, 1:B(var, 3), 3), ...
    1:B(var, 3), high_lim(var, 1:B(var, 3), 3));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f3) l = 0.005", "upper(f3) l = 0.005", "lower(f3) l = 0.01 " , "upper(f3) l = 0.01", "lower(f3) l = 0.1 ", "upper(f3) l = 0.1" );
