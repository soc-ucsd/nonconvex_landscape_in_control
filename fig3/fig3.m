% Figure 3: saddle point

close all

load LQGdata1.mat

fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [8, 6];
FontSize = 11;

surf(BK_grid, CK_grid, J_val,'EdgeColor', 'none'); hold on;
plot3(BK_1,CK_1,J_line1,'--','Color','m','linewidth',1.2)
plot3(BK_2,CK_2,J_line2,'--','Color','c','linewidth',1.2)
plot3(0,0,J_center,'o','Color','r','MarkerSize',6,'MarkerFaceColor','r')


xlabel('$B_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
ylabel('$C_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
zlabel('$J_{\texttt{LQG},1}(\mathsf{K})$','Interpreter','latex','FontSize',FontSize);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

view(20,20)
print(gcf,'LQG_saddle_xy.eps','-depsc2','-r300');


