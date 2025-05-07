close all;
clc;
clear;

V = 85:115;
%upper and lower bounds for x
lower = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
upper = [54.13, 21.56, 34.08, 49.19, 33.05, 21.84, 29.96, 24.87, ...
      47.24, 33.97, 26.89, 32.76, 39.98, 37.12, 53.83, 61.65, 59.73];
upper = upper - 0.002;

% Constraints
A = [-1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
       0, 0, 0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0; ...
       0, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
       1, 0, 0, 0, -1, -1, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0; ...
       0, 0, 1, 0, 0, 0, 0, 1, 1, 0, -1, -1, -1, 0, 0, 0, 0; ...
       0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0; ...
       0, 0, 0,0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, -1; ...
       0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0; ...
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1;];

x = zeros([length(V) 17]);
fval = zeros([length(V) 1]);
generations = zeros([length(V) 1]);

for k = 1:length(V)
    Vt = V(k);
    fprintf("Solving problem for V = %d\n", Vt)
    b = [-Vt; 0; 0; 0; 0; 0; 0; 0; Vt];
    objective_fun = @(x) objective_function(x, Vt);
    options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
    [xt, fvalt, exitflag, output, population, scores] = ...
        ga(objective_fun, 17, [], [], A, b, lower, upper, [], options);
    x(k, :) = xt(:);
    fval(k) = fvalt;
    generations(k) = output.generations;
end

%% Results

results = [V', generations, fval]; % Concatenate Volume, Generations, and Fitness
%split in +,- 15% for better visualization
V_low = V(V <= 100); 
results_low = results(V <= 100, :);
V_high = V(V > 100);
results_high = results(V > 100, :); 

%Volumes 85-100
figure;
set(gcf, 'Position', [100, 100, 600, 600]); 
axis off;
title('Genetic Algorithm Results for V = 85 to V = 100', 'FontSize', 14, 'FontWeight', 'bold');
hold on;

start_y = 0.9;
spacing = 0.05;
for i = 1:length(V_low)
    text(0.1, start_y - (i * spacing), sprintf('V = %.2f           Generations = %d           Best Fitness = %.3f', results_low(i, 1), results_low(i, 2), results_low(i, 3)), 'FontSize', 10, 'Interpreter', 'tex');
end
hold off;

% Volumes 101-115
figure;
set(gcf, 'Position', [100, 100, 600, 600]); % Set figure size
axis off;
title('Genetic Algorithm Results for V = 101 to V = 115', 'FontSize', 14, 'FontWeight', 'bold');
hold on;

start_y = 0.9;
spacing = 0.05;

for i = 1:length(V_high)
    text(0.1, start_y - (i * spacing), sprintf('V = %.2f           Generations = %d           Best Fitness = %.3f', results_high(i, 1), results_high(i, 2), results_high(i, 3)), 'FontSize', 10, 'Interpreter', 'tex');
end
hold off;

% Relation of V and fval
figure;
plot(V, fval, 'b-', 'LineWidth', 1.5);
xlabel('Volume V');
ylabel('Best Fitness Value');
title('Best Fitness Value for each V');
grid on;

%% Objective function

function time = objective_function(x, V)
    sum = 0;
    c = [54.13, 21.56, 34.08, 49.19, 33.05, 21.84, 29.96, 24.87, 47.24, 33.97, 26.89, 32.76, 39.98, 37.12, 53.83, 61.65, 59.73];
    ai = [1.25, 1.25, 1.25, 1.25, 1.25, 1.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1];
    ti = ones([1 17]) * 5;

    for k = 1:17
        sum = sum + x(k) * (ti(k) + ai(k) * x(k) / (1 - x(k)/c(k)));
    end

    time = sum / V;
end
