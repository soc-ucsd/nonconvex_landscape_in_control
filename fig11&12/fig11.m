%% Figure 11 in the paper

close all;

% Create a grid of epsilon and omega values
e_values = linspace(0.0001, 0.002, 200); % Adjust the range as needed
o_values = linspace(-0.000001, 0.000001, 1000);   % Adjust the range as needed
[e_grid, o_grid] = meshgrid(e_values, o_values);

% Compute the function values on the grid using vectorized operations
sigma_max = (3*e_grid.^4 + o_grid.^2 + e_grid.^4.*o_grid.^2 + ...
    sqrt(e_grid.^8.*o_grid.^4 + 6*e_grid.^8.*o_grid.^2 + 5*e_grid.^8 +...
    8*e_grid.^6.*o_grid.^2 - 2*e_grid.^4.*o_grid.^4 +...
    2*e_grid.^4.*o_grid.^2 + o_grid.^4)) ./ (2*e_grid.^4 - 4*e_grid.^2.*o_grid.^2 +...
    2*o_grid.^4 + 2*o_grid.^2);

o_grid1 = 0;
sigma_max0 = (3*e_values.^4 + o_grid1.^2 + e_values.^4.*o_grid1.^2 + ...
    sqrt(e_values.^8.*o_grid1.^4 + 6*e_values.^8.*o_grid1.^2 + 5*e_values.^8 +...
    8*e_values.^6.*o_grid1.^2 - 2*e_values.^4.*o_grid1.^4 +...
    2*e_values.^4.*o_grid1.^2 + o_grid1.^4)) ./ (2*e_values.^4 - 4*e_values.^2.*o_grid1.^2 +...
    2*o_grid1.^4 + 2*o_grid1.^2);

% Create a 3D surface plot
FontSize = 11;
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [12, 6];

%mesh(e_grid, o_grid, sigma_max);
surf(e_grid, o_grid, sigma_max,'EdgeColor', 'none'); hold on
plot3(e_values,zeros(length(e_values),1),sigma_max0,'-','Color','m','linewidth',1.2);
xlabel('$\epsilon$','Interpreter','latex','FontSize',FontSize);
ylabel('$\omega$','Interpreter','latex','FontSize',FontSize);
zlabel('$\sigma_{\max}^2$','Interpreter','latex','FontSize',FontSize);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

%title('Shape of sigma_max Function');
grid off
colorbar;
view(-60,30)
print(gcf,'Hinf_sigma.eps','-depsc2','-r300');

