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


% Compute derivatives
df1 = diff(f1, x);
df2 = diff(f2, x);
df3 = diff(f3, x);

% Convert symbolic derivatives to function handles
df1_func = matlabFunction(df1);
df2_func = matlabFunction(df2);
df3_func = matlabFunction(df3);


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


%varying l parameters
l_steps = 100;
l_first = 0.001;
l_last = 0.1;

ls = linspace(l_first, l_last, l_steps);

B = zeros([l_steps 3]);
results = zeros([l_steps 3]);
low_lim = zeros([l_steps 1 3]);
high_lim = zeros([l_steps 1 3]);

for i = 1:l_steps
    
    l = ls(i);

    
    iterations_1 = 1;
    iterations_2 = 1;
    iterations_3 = 1;

   
    calculations_1 = 0;
    calculations_2 = 0;
    calculations_3 = 0;

    
    lower_1 = lower;
    upper_1 = upper;

    lower_2 = lower;
    upper_2 = upper;

    lower_3 = lower;
    upper_3 = upper;

    % initialize limit storing
    low_lim(i, 1, :) = [-1 -1 -1];
    high_lim(i, 1, :) = [3 3 3];


    % f1
    while upper_1 - lower_1 >= l
        midpoint = (upper_1 + lower_1) / 2;
        der_f1 = df1_func(midpoint);
        if der_f1 > 0
            upper_1 = midpoint;
        elseif der_f1 < 0
            lower_1 = midpoint;
        else
           % break the loop
           upper_1 = midpoint;
           lower_1 = midpoint;
        end

        iterations_1 = iterations_1 + 1;
        calculations_1 = calculations_1 + 1;

        low_lim(i, iterations_1, 1) = lower_1;
        high_lim(i, iterations_1, 1) = upper_1;
    end

    % f2
    while upper_2 - lower_2 >= l
        midpoint = (upper_2 + lower_2) / 2;
        der_f2 = df2_func(midpoint);
        if der_f2 > 0
           upper_2 = midpoint;
        elseif der_f2 < 0
           lower_2 = midpoint;
        else
            % break the loop
            upper_2 = midpoint;
            lower_2 = midpoint;
        end

        iterations_2 = iterations_2 + 1;
        calculations_2 = calculations_2 + 1;

        low_lim(i, iterations_2, 2) = lower_2;
        high_lim(i, iterations_2, 2) = upper_2;
    end

    % f3
    while upper_3 - lower_3 >= l
        midpoint = (upper_3 + lower_3) / 2;
        der_f3 = df3_func(midpoint);
        if der_f3 > 0
           upper_3 = midpoint;
        elseif der_f3 < 0
           lower_3 = midpoint;
        else
            % break the loop
            upper_3 = midpoint;
            lower_3 = midpoint;
        end

        iterations_3 = iterations_3 + 1;
        calculations_3 = calculations_3 + 1;

        low_lim(i, iterations_3, 3) = lower_3;
        high_lim(i, iterations_3, 3) = upper_3;
    end

    % save the result
    B(i, :) = [iterations_1 iterations_2 iterations_3];
    results(i, :) = [calculations_1 calculations_2 calculations_3];
end


figure("Name", "Μεταβολή των υπολογισμών")
plot(ls, results(:, 1), ls, results(:, 2),ls, results(:, 3) );
title('Μεταβολή υπολογισμών αντικειμενικής συνάρτησης ');
subtitle('e = 0.001 και μεταβαλλόμενο l');
xlabel("l");
ylabel("Μεταβολή των υπολογισμών της f(x)");
legend("f1(x)", "f2(x)", "f3(x)");
hold on;



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


