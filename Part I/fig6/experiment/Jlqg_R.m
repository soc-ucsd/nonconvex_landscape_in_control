clear all; close all; clc;

R = [0.01,0.1,1,10];
J_opt = zeros(1,length(R));
J_hifoo = zeros(1,length(R));
J_GD = zeros(1,length(R));
J_mincx = zeros(1,length(R));
J_solver = zeros(1,length(R));

for i = 1:length(R)
    [J_opt(i), J_hifoo(i), J_GD(i), J_mincx(i), J_solver(i)] ...
        = Jopt_LQG_with_diff_R(R(i));
end

% plot(R,J_opt,R,J_solver,R,J_mincx)
% 
% xlabel('R','FontSize',20,'Interpreter','latex')
% %ylabel('$J_{\texttt{lqg},2}(\mathsf{K})-J^\star)$','FontSize',20,'Interpreter','latex')
% grid on;
% set(gca,'FontSize',20)
% %title('HIFOO');
% legend('Opt','Mosek','mincx')