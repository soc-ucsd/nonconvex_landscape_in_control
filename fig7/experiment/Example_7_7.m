clear all; close all; clc;

rho = [0.01,0.1,1,10];
J_hifoo = zeros(1,length(rho));
J_riccati = zeros(1,length(rho));
J_lmi = zeros(1,length(rho));
J_mincx = zeros(1,length(rho));
J_solver = zeros(1,length(rho));
J_cvx = zeros(1,length(rho));

for i = 1:length(rho)
    [J_hifoo(i), J_riccati(i), J_lmi(i), J_mincx(i), J_solver(i), J_cvx(i)] = Jopt_Hinf_with_diff_rho(rho(i));
end

%%
plot(rho,J_hifoo,rho,J_riccati,rho,J_lmi,rho,J_mincx);

xlabel('\rho','FontSize',20,'Interpreter','latex')
%ylabel('$J_{\texttt{lqg},2}(\mathsf{K})-J^\star)$','FontSize',20,'Interpreter','latex')
grid on;
set(gca,'FontSize',20)
%title('HIFOO');
legend('hifoo','Ric','lmi','mincx')