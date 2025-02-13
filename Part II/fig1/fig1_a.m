 
% Fig 1(a) in "Benign Nonconvex Landscapes in Optimal and Robust Control,
%           Part II: Extended Convex Lifting"
%       By Yang Zheng, Chih-Fan Pai, Yujie Tang

% Define the range for x and y values
x = linspace(-5, 5, 500);
y = linspace(-5, 5, 500);

% Create a grid of x and y values
[X, Y] = meshgrid(x, y);

% Evaluate the inequality xy < 1
Z = X .* Y < 1;

% Plot the set where the inequality holds true
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [7, 6]; %[-k,-k] achieve min
FontSize = 11;
contourf(X, Y, Z, [0 1], 'k--','LineWidth', 2);
xlabel('$k_1$','Interpreter','latex','FontSize',FontSize);
ylabel('$k_2$','Interpreter','latex','FontSize',FontSize);
% colormap([0 0 0; 1 1 1]); % Custom colormap for better visualization
custom_colormap = [
    1 1 1;   % White
    0 0.85 0.9;   % Blue
];
colormap(custom_colormap);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');
