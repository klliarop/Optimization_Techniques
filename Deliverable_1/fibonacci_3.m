% Clear workspace, close all figures, and clear command window
close all;
clear all;
clc;

% Define symbolic functions f1(x), f2(x), and f3(x)
syms f1(x);
syms f2(x);
syms f3(x);
f1(x) = (x - 2)^2 + x*log(x + 3);
f2(x) = exp(-2*x) + (x - 2)^2 ; 
f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x);

% Define bounds for x
lower = -1;
upper = 3;

% Plot of function f1(x)

figure("Name", "f1(x)")
fplot(f1, [lower, upper]);
title('f1(x) = (x - 2)^2 + x*log(x + 3)');
xlabel("x");
ylabel("y");
hold on;

% Plot of function f2(x)
figure("Name", "f2(x)")
fplot(f2, [lower, upper]);
title('f2(x) = exp(-2*x) + (x - 2)^2');
xlabel("x");
ylabel("y");
hold on;

% Plot of function f3(x)
figure("Name", "f3(x)")
fplot(f3, [lower, upper]);
title('f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x)');
xlabel("x");
ylabel("y");
% Η παρακάτω γραμμή απεικονίζει καλύτερα την αρχική συμπεριφορά της
% συνάρτησης καθώς διαφορετικά δεν φαίνονται οι λεπτομέρειες
%ylim([-5, 15]);
yline(0, 'k'); % Horizontal line at y=0
hold on;

% Set up constant l and varying e parameters for method
n_l = 100;
first_l = 0.001;
last_l = 0.1;

ls = linspace(first_l, last_l, n_l);

reps = zeros([n_l 3]);
results = zeros([n_l 3]);
lower_limits = zeros([n_l 1 3]);
upper_limits = zeros([n_l 1 3]);


% Perform fibonacci method for each value of e
for i = 1:n_l
   
    l = ls(i);
    epsilon = 0.0001;

    compared = (upper - lower) / l;
    n_reps = 1;
    while fibonacci(n_reps) <= compared
        n_reps = n_reps + 1;
    end

    
    a = lower;
    upper_1 = upper;

    lower_2 = lower;
    upper_2 = upper;

    lower_3 = lower;
    upper_3 = upper;

   
    lower_limits(i, 1, :) = [0 0 0];
    upper_limits(i, 1, :) = [3 3 3];

    
    % f1
    x_k_1 = a + (fibonacci(n_reps - 2) / fibonacci(n_reps)) * (upper_1 - a);
    x_k_2 = a + (fibonacci(n_reps - 1) / fibonacci(n_reps)) * (upper_1 - a);
    lower_value_1 = f1(x_k_1);
    upper_value_1 = f1(x_k_2);
    iterations_1 = 1;
    while iterations_1 < n_reps - 2
        if lower_value_1 > upper_value_1
            a = x_k_1;
            x_k_1 = x_k_2;
            lower_value_1 = upper_value_1;
            x_k_2 = a + (fibonacci(n_reps - iterations_1 - 1) / fibonacci(n_reps - iterations_1)) * (upper_1 - a);
            upper_value_1 = f1(x_k_2);
        else
            upper_1 = x_k_2;
            x_k_2 = x_k_1;
            upper_value_1 = lower_value_1;
            x_k_1 = a + (fibonacci(n_reps - iterations_1 - 2) / fibonacci(n_reps - iterations_1)) * (upper_1 - a);
            lower_value_1 = f1(x_k_1);
        end
        iterations_1 = iterations_1 + 1;
        lower_limits(i, iterations_1, 1) = a;
        upper_limits(i, iterations_1, 1) = upper_1;
    end
    if lower_value_1 > upper_value_1
        a = x_k_1;
        x_k_1 = x_k_2;
        lower_value_1 = upper_value_1;
%         x_k_2 = lower_1 + (fibonacci(n_iterations - iterations_1 - 1) / fibonacci(n_iterations - iterations_1)) * (upper_1 - lower_1);
    else
        upper_1 = x_k_2;
%         x_k_2 = x_k_1;
        upper_value_1 = lower_value_1;
        x_k_1 = a + (fibonacci(n_reps - iterations_1 - 2) / fibonacci(n_reps - iterations_1)) * (upper_1 - a);
    end
    lower_limits(i, iterations_1 + 1, 1) = a;
    upper_limits(i, iterations_1 + 1, 1) = upper_1;
    if f1(x_k_1) > f1(x_k_1 + epsilon)
        a = x_k_1;
    else
        upper_1 = x_k_1 + epsilon;
    end
    lower_limits(i, iterations_1 + 2, 1) = a;
    upper_limits(i, iterations_1 + 2, 1) = upper_1;
    
    % f2
    x_k_1 = lower_2 + (fibonacci(n_reps - 2) / fibonacci(n_reps)) * (upper_2 - lower_2);
    x_k_2 = lower_2 + (fibonacci(n_reps - 1) / fibonacci(n_reps)) * (upper_2 - lower_2);
    lower_value_2 = f2(x_k_1);
    upper_value_2 = f2(x_k_2);
    iterations_2 = 1;
    while iterations_2 < n_reps - 2
        if lower_value_2 > upper_value_2
            lower_2 = x_k_1;
            x_k_1 = x_k_2;
            lower_value_2 = upper_value_2;
            x_k_2 = lower_2 + (fibonacci(n_reps - iterations_2 - 1) / fibonacci(n_reps - iterations_2)) * (upper_2 - lower_2);
            upper_value_2 = f2(x_k_2);
        else
            upper_2 = x_k_2;
            x_k_2 = x_k_1;
            upper_value_2 = lower_value_2;
            x_k_1 = lower_2 + (fibonacci(n_reps - iterations_2 - 2) / fibonacci(n_reps - iterations_2)) * (upper_2 - lower_2);
            lower_value_2 = f2(x_k_1);
        end
        iterations_2 = iterations_2 + 1;
        lower_limits(i, iterations_2, 2) = lower_2;
        upper_limits(i, iterations_2, 2) = upper_2;
    end
    if lower_value_2 > upper_value_2
        lower_2 = x_k_1;
        x_k_1 = x_k_2;
        lower_value_2 = upper_value_2;
    else
        upper_2 = x_k_2;
        upper_value_2 = lower_value_2;
        x_k_1 = lower_2 + (fibonacci(n_reps - iterations_2 - 2) / fibonacci(n_reps - iterations_2)) * (upper_2 - lower_2);
    end
    lower_limits(i, iterations_2 + 1, 2) = lower_2;
    upper_limits(i, iterations_2 + 1, 2) = upper_2;
    if f2(x_k_1) > f2(x_k_1 + epsilon)
        lower_2 = x_k_1;
    else
        upper_2 = x_k_1 + epsilon;
    end
    lower_limits(i, iterations_2 + 2, 2) = lower_2;
    upper_limits(i, iterations_2 + 2, 2) = upper_2;

    % f3
    x_k_1 = lower_3 + (fibonacci(n_reps - 2) / fibonacci(n_reps)) * (upper_3 - lower_3);
    x_k_2 = lower_3 + (fibonacci(n_reps - 1) / fibonacci(n_reps)) * (upper_3 - lower_3);
    lower_value_3 = f3(x_k_1);
    upper_value_3 = f3(x_k_2);
    iterations_3 = 1;
    while iterations_3 < n_reps - 2
        if lower_value_3 > upper_value_3
            lower_3 = x_k_1;
            x_k_1 = x_k_2;
            lower_value_3 = upper_value_3;
            x_k_2 = lower_3 + (fibonacci(n_reps - iterations_3 - 1) / fibonacci(n_reps - iterations_3)) * (upper_3 - lower_3);
            upper_value_3 = f1(x_k_2);
        else
            upper_3 = x_k_2;
            x_k_2 = x_k_1;
            upper_value_3 = lower_value_3;
            x_k_1 = lower_3 + (fibonacci(n_reps - iterations_3 - 2) / fibonacci(n_reps - iterations_3)) * (upper_3 - lower_3);
            lower_value_3 = f3(x_k_1);
        end
        iterations_3 = iterations_3 + 1;
        lower_limits(i, iterations_3, 3) = lower_3;
        upper_limits(i, iterations_3, 3) = upper_3;
    end
    if lower_value_3 > upper_value_3
        lower_3 = x_k_1;
        x_k_1 = x_k_2;
        lower_value_3 = upper_value_3;
        x_k_2 = lower_3 + (fibonacci(n_reps - iterations_3 - 1) / fibonacci(n_reps - iterations_3)) * (upper_3 - lower_3);
    else
        upper_3 = x_k_2;
        x_k_2 = x_k_1;
        upper_value_3 = lower_value_3;
        x_k_1 = lower_3 + (fibonacci(n_reps - iterations_3 - 2) / fibonacci(n_reps - iterations_3)) * (upper_3 - lower_3);
    end
    lower_limits(i, iterations_3 + 1, 3) = lower_3;
    upper_limits(i, iterations_3 + 1, 3) = upper_3;
    if f3(x_k_1) > f3(x_k_1 + epsilon)
        lower_3 = x_k_1;
    else
        upper_3 = x_k_1 + epsilon;
    end
    lower_limits(i, iterations_3 + 2, 3) = lower_3;
    upper_limits(i, iterations_3 + 2, 3) = upper_3;

    % save the result
    reps(i, :) = [n_reps n_reps n_reps];
    results(i, :) = [n_reps n_reps n_reps];
end


figure("Name", "")
plot(ls, results(:, 1), ls, results(:, 2), ls, results(:, 3) );
title('Μεταβολή υπολογισμών αντικειμενικής συνάρτησης ');
subtitle('e = 0.001 και μεταβαλλόμενο l');
xlabel("l");
ylabel("Μεταβολή των υπολογισμών της f(x)");
legend("f1(x)", "f2(x)", "f3(x)");
hold on;

% Plotting the upper and lower limits for f1, f2, and f3

var = 1;
figure("Name", "f1(x)")
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
    1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
hold on;
var =  n_l / 2;
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
    1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
hold on;
var =  n_l;
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
    1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f1) l = 0.005", "upper(f1) l = 0.005", "lower(f1) l = 0.01 " , "upper(f1) l = 0.01", "lower(f1) l = 0.1 ", "upper(f1) l = 0.1" );




var = 1;
figure("Name", "f2(x)")
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
    1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
hold on;
var =  n_l / 2;
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
    1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
hold on;
var =  n_l;
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
    1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f2) l = 0.005", "upper(f2) l = 0.005", "lower(f2) l = 0.01 " , "upper(f2) l = 0.01", "lower(f2) l = 0.1 ", "upper(f2) l = 0.1" );


var = 1;
figure("Name", "f3(x)")
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
    1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
hold on;
var =  n_l / 2;
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
    1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
hold on;
var =  n_l;
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
    1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
hold on;

xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f3) l = 0.005", "upper(f3) l = 0.005", "lower(f3) l = 0.01 " , "upper(f3) l = 0.01", "lower(f3) l = 0.1 ", "upper(f3) l = 0.1" );
