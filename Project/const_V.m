close all;
clear; 
clc;

V = 100;
% Upper and lower bounds for x
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
b = [-V; 0; 0; 0; 0; 0; 0; 0; V];

% Solve the problem
objective_fun = @(x) objective_function(x, V);
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf);
[x, fval, exitflag, output, population, scores] = ...
    ga(objective_fun, 17, [], [], A, b, lower, upper, [], options);

%% Results
best_fitness = min(scores);  
mean_fitness = mean(scores);
figure;
set(gcf, 'Position', [100, 100, 600, 600]);
axis off;
title('Genetic Algorithm Results', 'FontSize', 14, 'FontWeight', 'bold');
hold on;
text(0.1, 0.9, sprintf('Best Fitness: %.3f', best_fitness), 'FontSize', 12, 'FontWeight', 'bold');
text(0.55, 0.9, sprintf('Mean Fitness: %.3f', mean_fitness), 'FontSize', 12, 'FontWeight', 'bold'); % Adjusted to the right
start_y = 0.8; 
spacing = 0.04; 
for i = 1:length(x)
    text(0.1, start_y - (i * spacing), sprintf('x_{%d} = %.3f', i, x(i)), 'FontSize', 10, 'Interpreter', 'tex');
end
hold off;

disp(['Number of generations: ', num2str(output.generations)]);

%% Plot of Fitness Value per generation
% Because of some disproportionally big starting values 
% a logarithmic scale is set to the vertical axis 

scores = scores(:)';
figure;
plot(scores, 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
xlabel('Generation');
ylabel('Fitness Value');
title('Fitness Value per Generation');
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
