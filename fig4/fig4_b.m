
close all

load hinfb.mat

fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [7, 6];
FontSize = 11;


surf(Bk_grid, Ck_grid, J_val,'EdgeColor', 'none','FaceAlpha',0.9); 
hold on;
plot3(Bk_vals_cur, Ck_vals_cur, J_val_cur,'-','Color','r','linewidth',1.2)
plot3(Bk_vals_cur2, Ck_vals_cur2, J_val_cur2,'-','Color','r','linewidth',1.2)
%plot3(0,0,J_center,'o','Color','r','MarkerSize',6,'MarkerFaceColor','r')


xlabel('$B_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
ylabel('$C_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
zlabel('$J_{\infty,1}(\mathsf{K})$','Interpreter','latex','FontSize',FontSize);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

view(-80,25)
print(gcf,'Hinf_b.eps','-depsc2','-r300');


