

close all; clear all; clc;
%%
a = -1; % -1,0,1
gamma = sqrt(a^2+2)+a;
DK_star = -(a+sqrt(a^2+1-gamma^(-2)))/(1-gamma^(-2));

A = a;
B2 = 1;
C2 = 1;

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

B1 = [eye(nx) zeros(nx, ny)];
C1 = [eye(nx); zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); eye(nu)];
D21 = [zeros(ny, nx) eye(ny)];

sys_pl = ss(A, [B1 B2], [C1; C2], [D11 D12; D21 zeros(ny, nu)]);
%[K_opt, sys_cl_opt, J_opt] = hinfsyn(sys_pl, ny, nu);


DK_vals = DK_star-0.6 : 5e-2 : DK_star+1;
n_ax = length(DK_vals);

eps = 1e-10;
J_val = zeros(1, n_ax);
for i = 1:n_ax
    AK = -1;
    BK = 0;
    CK = 0;
    DK = DK_vals(i); % -0.7321

    Acl = [A + B2*DK*C2, B2*CK;
        BK*C2, AK];
    Bcl = [B1 + B2*DK*D21; BK*D21];
    Ccl = [C1 + D12*DK*C2, D12*CK];
    Dcl = D11 + D12*DK*D21;

    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
    J_val(i) = hinfnorm(sys_cl, eps);
end
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [6.5, 6];
FontSize = 11;


plot(DK_vals, J_val,'m','linewidth',1.2); hold on;grid on;
plot(DK_star,gamma,'o','Color','r','MarkerSize',4,'MarkerFaceColor','r')
xlabel('$D_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize); 
ylabel('$J_{\infty,0}(D_{\mathsf{K}})$','Interpreter','latex','FontSize',FontSize);
ylim([0.6,1.5]); xlim([-1.5,0.5])
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');
print(gcf,'Hinf_a.eps','-depsc2','-r300');
