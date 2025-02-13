close all;
clc


% Define a range of epsilon values, including values close to zero
epsilons = logspace(-4, -1, 50);  % Range from 1e-6 to 1

s = tf('s');
H_infinity_norms1 = zeros(1, numel(epsilons));
H_infinity_norms2 = zeros(1, numel(epsilons));
H_infinity_norms3 = zeros(1, numel(epsilons));
H_infinity_norms4 = zeros(1, numel(epsilons));

opts = 10^(-8);
% Compute the H_infinity norm for different epsilons
for i = 1:numel(epsilons)
    epsilon = epsilons(i);
    G1 = 1 / (s^2 + (1 - epsilon^2) * s + epsilon^5) * [s - epsilon^2, -epsilon^2 - epsilon^5; ...
                                                     -epsilon^2 - epsilon^5, - (epsilon^2 + epsilon^5) * (s + 1)];
    
    G2 = 1 / (s^2 + (1 - epsilon^2) * s + epsilon^4) * [s - epsilon^2, -epsilon^2 - epsilon^4; ...
                                                     -epsilon^2 - epsilon^4, - (epsilon^2 + epsilon^4) * (s + 1)];
    G3 = 1 / (s^2 + (1 - epsilon^2) * s + epsilon^3) * [s - epsilon^2, -epsilon^2 - epsilon^3; ...
                                                     -epsilon^2 - epsilon^3, - (epsilon^2 + epsilon^3) * (s + 1)];
    
    G4 = 1 / (s^2 + (1 - epsilon^2) * s + epsilon^2) * [s - epsilon^2, -epsilon^2 - epsilon^2; ...
                                                     -epsilon^2 - epsilon^2, - (epsilon^2 + epsilon^2) * (s + 1)];
                                            
    % Compute the H_infinity norm with specified options
    H_infinity_norms1(i) = norm(G1, inf, opts);
    H_infinity_norms2(i) = norm(G2, inf, opts);
    H_infinity_norms3(i) = norm(G3, inf, opts);
    H_infinity_norms4(i) = norm(G4, inf, opts);
end

%% estimation 
% Load or define the data from the previous computation
% epsilon values and corresponding H_infinity norm values

% Perform linear regression to estimate the relationship
X = log(epsilons);  % Log-transform of epsilon
Y1 = log(H_infinity_norms1);  % Log-transform of H_infinity norm
coefficients1 = polyfit(X, Y1, 1);  % Fit a linear model y = mx + c

Y2 = log(H_infinity_norms2);  % Log-transform of H_infinity norm
coefficients2 = polyfit(X, Y2, 1);  % Fit a linear model y = mx + c

Y3 = log(H_infinity_norms3);  % Log-transform of H_infinity norm
coefficients3 = polyfit(X, Y3, 1);  % Fit a linear model y = mx + c

Y4 = log(H_infinity_norms4);  % Log-transform of H_infinity norm
coefficients4 = polyfit(X, Y4, 1);  % Fit a linear model y = mx + c


% Extract the slope (m) and intercept (c) from the coefficients
m1 = coefficients1(1); c1 = coefficients1(2);
m2 = coefficients2(1); c2 = coefficients2(2);
m3 = coefficients3(1); c3 = coefficients3(2);
m4 = coefficients4(1); c4 = coefficients4(2);

% Compute the estimated function y = m*x + c in log-log space
Y_estimated1 = m1 * X + c1;
Y_estimated2 = m2 * X + c2;
Y_estimated3 = m3 * X + c3;
Y_estimated4 = m4 * X + c4;

% Plot the original data and the linear regression estimate
FontSize = 12;
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [10, 6];


hold on;
loglog(epsilons(1:5:end), H_infinity_norms1(1:5:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
l1 = loglog(epsilons, exp(Y_estimated1), 'r', 'LineWidth', 2);
%hold on;

hold on
loglog(epsilons(1:5:end), H_infinity_norms2(1:5:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
l2 = loglog(epsilons, exp(Y_estimated2), 'k', 'LineWidth', 2);

hold on
loglog(epsilons(1:5:end), H_infinity_norms3(1:5:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
l3 = loglog(epsilons, exp(Y_estimated3), 'g', 'LineWidth', 2);

%hold on
%loglog(epsilons(1:5:end), H_infinity_norms4(1:5:end), 'o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
%loglog(epsilons, exp(Y_estimated4), 'b', 'LineWidth', 2);

set(gca,'YScale','log','XScale','log');
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

xlabel('$\epsilon$','Interpreter','latex','FontSize',FontSize);
ylabel('$\mathcal{H}_\infty$ norm','Interpreter','latex','FontSize',FontSize);

h = legend([l1,l2,l3],'$d = 5$','$d = 4$','$d = 3$','Interpreter','latex','FontSize',FontSize);
set(h,'box','off')
grid on;

print(gcf,'Hinf_boundary_2.eps','-depsc2','-r300');


