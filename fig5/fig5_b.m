close all; clear all; clc;
%%
A = -1;
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

eps_vals = 1e-3:1e-2:1;
n_eps = length(eps_vals);
n_ay = n_eps;

Bk_vals = -eps_vals;
Ck_vals = eps_vals;
[Bk_grid,Ck_grid] = meshgrid(Bk_vals, Ck_vals);

eps = 1e-9;
J_val = zeros(n_ay,1);
for jj = 1:n_ay
    AK = 0;
    BK = Bk_vals(jj);
    CK = Ck_vals(jj);
    DK = 0; % -0.2

    Acl = [A + B2*DK*C2, B2*CK;
        BK*C2, AK];
    Bcl = [B1 + B2*DK*D21; BK*D21];
    Ccl = [C1 + D12*DK*C2, D12*CK];
    Dcl = D11 + D12*DK*D21;

    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
    J_val(jj,1) = hinfnorm(sys_cl, eps);
end
%%
% figure;
% plot(eps_vals,J_val);
% xlabel('$\epsilon$','Interpreter','latex');
% title('$J(K_\epsilon)$','Interpreter','latex');
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [6.5, 6];
FontSize = 11;

plot(eps_vals,J_val,'m','linewidth',1.2); hold on;grid on;
% plot(DK_star,gamma,'o','Color','r','MarkerSize',4,'MarkerFaceColor','r')
xlabel('$\epsilon$','Interpreter','latex','FontSize',FontSize); 
ylabel('$J_{\infty,1}(\mathsf{K}_1(\epsilon))$','Interpreter','latex','FontSize',FontSize);
ylim([1.55,2.4]); xlim([0,1])
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');
print(gcf,'fig5_b.eps','-depsc2','-r300');
