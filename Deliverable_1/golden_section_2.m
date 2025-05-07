close all;     % Κλείσιμο όλων των παραθύρων γραφημάτων
clear all;     % Καθαρισμός όλων των μεταβλητών από το workspace
clc;           % Καθαρισμός της κονσόλας

syms f1(x);    % Ορισμός της f1(x) ως συμβολική συνάρτηση
syms f2(x);    % Ορισμός της f2(x) ως συμβολική συνάρτηση
syms f3(x);    % Ορισμός της f3(x) ως συμβολική συνάρτηση

% Ορισμός των συναρτήσεων f1, f2 και f3
f1(x) = (x - 2)^2 + x*log(x + 3);
f2(x) = exp(-2*x) + (x - 2)^2;
f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x);

lower = -1;    % Κάτω όριο του x για το διάστημα αναζήτησης
upper = 3;     % Άνω όριο του x για το διάστημα αναζήτησης

% Διαγράμματα των συναρτήσεων f1, f2 και f3
figure("Name", "f1(x)");
fplot(f1, [lower, upper]);
title('f1(x) = (x - 2)^2 + x*log(x + 3)');
xlabel("x");
ylabel("y");
hold on;

figure("Name", "f2(x)");
fplot(f2, [lower, upper]);
title('f2(x) = exp(-2*x) + (x - 2)^2');
xlabel("x");
ylabel("y");
hold on;

figure("Name", "f3(x)");
fplot(f3, [lower, upper]);
title('f3(x) = exp(x)*(x^3 - 1) + (x - 1)*sin(x)');
xlabel("x");
ylabel("y");
% Η επόμενη γραμμή για καλύτερη απεικόνιση αρχικής συμπεριφοράς
%ylim([-5, 15]);
yline(0, 'k'); % Οριζόντια γραμμή στο y=0
hold on;

% Μεταβλητές για τα βήματα του l
l_steps = 100;
l_first = 0.005;
l_last = 0.1;

ls = linspace(l_first, l_last, l_steps);

% Πίνακες για αποθήκευση επαναλήψεων και υπολογισμών
reps = zeros([l_steps 3]);
calculations = zeros([l_steps 3]);
lower_limits = zeros([l_steps 1 3]);
upper_limits = zeros([l_steps 1 3]);

for i = 1:l_steps
    l = ls(i); % Επιλογή του τρέχοντος l από τη σειρά ls

    % Αρχικοποίηση μεταβλητών για κάθε συνάρτηση
    reps_1 = 1; reps_2 = 1; reps_3 = 1;
    calculations_1 = 0; calculations_2 = 0; calculations_3 = 0;
    gamma = 0.618;

    lower_limits(i, 1, :) = [0 0 0];
    upper_limits(i, 1, :) = [3 3 3];

    % Εύρεση βέλτιστου διαστήματος για f1(x)
    a = lower; b = upper;
    lower_value_1 = f1(a + (1 - gamma) * (b - a));
    upper_value_1 = f1(a + gamma * (b - a));
    while b - a >= l
        if lower_value_1 > upper_value_1
           a = a + (1 - gamma) * (b - a);
           lower_value_1 = upper_value_1;
           upper_value_1 = f1(a + gamma * (b - a));
        else
            b = a + gamma * (b - a);
            upper_value_1 = lower_value_1;
            lower_value_1 = f1(a + (1 - gamma) * (b - a));
        end    
        reps_1 = reps_1 + 1;
        calculations_1 = calculations_1 + 1;

        lower_limits(i, reps_1, 1) = a;
        upper_limits(i, reps_1, 1) = b;
    end

    % Εύρεση βέλτιστου διαστήματος για f2(x)
    a = lower; b = upper;
    lower_value_2 = f2(a + (1 - gamma) * (b - a));
    upper_value_2 = f2(a + gamma * (b - a));
    while b - a >= l
        if lower_value_2 > upper_value_2
           a = a + (1 - gamma) * (b - a);
           lower_value_2 = upper_value_2;
           upper_value_2 = f2(a + gamma * (b - a));
        else
            b = a + gamma * (b - a);
            upper_value_2 = lower_value_2;
            lower_value_2 = f2(a + (1 - gamma) * (b - a));
        end    
        reps_2 = reps_2 + 1;
        calculations_2 = calculations_2 + 1;

        lower_limits(i, reps_2, 2) = a;
        upper_limits(i, reps_2, 2) = b;
    end

    % Εύρεση βέλτιστου διαστήματος για f3(x)
    a = lower; b = upper;
    lower_value_3 = f3(a + (1 - gamma) * (b - a));
    upper_value_3 = f3(a + gamma * (b - a));
    while b - a >= l
        if lower_value_3 > upper_value_3
           a = a + (1 - gamma) * (b - a);
           lower_value_3 = upper_value_3;
           upper_value_3 = f3(a + gamma * (b - a));
        else
            b = a + gamma * (b - a);
            upper_value_3 = lower_value_3;
            lower_value_3 = f3(a + (1 - gamma) * (b - a));
        end    
        reps_3 = reps_3 + 1;
        calculations_3 = calculations_3 + 1;

        lower_limits(i, reps_3, 3) = a;
        upper_limits(i, reps_3, 3) = b;
    end

    % Αποθήκευση αποτελεσμάτων για κάθε τιμή του l
    reps(i, :) = [reps_1 reps_2 reps_3];
    calculations(i, :) = [calculations_1 calculations_2 calculations_3];
end

% Γράφημα υπολογισμών αντικειμενικής συνάρτησης
figure("Name", "Μεταβολή των υπολογισμών αντικειμενικής συνάρτησης");
plot(ls, calculations(:, 1), ls, calculations(:, 2), ls, calculations(:, 3));
title('Μεταβολή υπολογισμών αντικειμενικής συνάρτησης');
subtitle('e = 0.001 και μεταβαλλόμενο l');
xlabel("l");
ylabel("calculations");
legend("f1", "f2", "f3");
hold on;

% Διαδικασία απεικόνισης διαστήματος για κάθε συνάρτηση με διαφορετικά l

% Διάγραμμα για f1(x)
var = 1;
figure("Name", "f1(x)");
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
     1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
hold on;
var = l_steps / 2;
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
     1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
hold on;
var = l_steps;
plot(1:reps(var, 1), lower_limits(var, 1:reps(var, 1), 1), ...
     1:reps(var, 1), upper_limits(var, 1:reps(var, 1), 1));
xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f1) l = 0.005", "upper(f1) l = 0.005", ...
       "lower(f1) l = 0.01", "upper(f1) l = 0.01", ...
       "lower(f1) l = 0.1", "upper(f1) l = 0.1");

% Διάγραμμα για f2(x)
var = 1;
figure("Name", "f2(x)");
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
     1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
hold on;
var = l_steps / 2;
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
     1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
hold on;
var = l_steps;
plot(1:reps(var, 2), lower_limits(var, 1:reps(var, 2), 2), ...
     1:reps(var, 2), upper_limits(var, 1:reps(var, 2), 2));
xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f2) l = 0.005", "upper(f2) l = 0.005", ...
       "lower(f2) l = 0.01", "upper(f2) l = 0.01", ...
       "lower(f2) l = 0.1", "upper(f2) l = 0.1");

% Διάγραμμα για f3(x)
var = 1;
figure("Name", "f3(x)");
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
     1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
hold on;
var = l_steps / 2;
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
     1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
hold on;
var = l_steps;
plot(1:reps(var, 3), lower_limits(var, 1:reps(var, 3), 3), ...
     1:reps(var, 3), upper_limits(var, 1:reps(var, 3), 3));
xlabel("k");
ylabel("Upper and Lower limit");
legend("lower(f3) l = 0.005", "upper(f3) l = 0.005", ...
       "lower(f3) l = 0.01", "upper(f3) l = 0.01", ...
       "lower(f3) l = 0.1", "upper(f3) l = 0.1");
