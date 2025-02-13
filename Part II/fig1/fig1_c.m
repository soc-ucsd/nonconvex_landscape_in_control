% Fig 1(c) in "Benign Nonconvex Landscapes in Optimal and Robust Control,
%           Part II: Extended Convex Lifting"
%       By Yang Zheng, Chih-Fan Pai, Yujie Tang


close all; clear all; clc;
%%
a = 1; % -1,0,1
gamma_opt = sqrt(a^2+2)+a;
DK_star = -(a+sqrt(a^2+1-gamma_opt^(-2)))/(1-gamma_opt^(-2));

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

n_ax = 101;
n_ay = 101;

Bk_vals = linspace(-0.65,0.65,n_ax); % around BK=0
Ck_vals = linspace(-0.65,0.65,n_ay); % around CK=0
[Bk_grid,Ck_grid] = meshgrid(Bk_vals, Ck_vals);

eps = 1e-9;
J_val = zeros(n_ax, n_ay);
for ii = 1:n_ax
    for jj = 1:n_ay
        AK = -1;
        BK = Bk_grid(ii, jj);
        CK = Ck_grid(ii, jj);
        DK = DK_star;

        Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
        J_val(ii, jj) = hinfnorm(sys_cl, eps);
    end
end
%%
Bk_vals_cur = linspace(-0.65,0.65,n_ax);
Ck_vals_cur = zeros(length(Bk_vals_cur),1);

n_eps_cur = length(Bk_vals_cur);

eps = 1e-9;
J_val_cur = zeros(n_eps_cur,1);
for ii = 1:n_eps_cur
    AK = -1;
    BK = Bk_vals_cur(ii);
    CK = Ck_vals_cur(ii);
    DK = DK_star;

    Acl = [A + B2*DK*C2, B2*CK;
        BK*C2, AK];
    Bcl = [B1 + B2*DK*D21; BK*D21];
    Ccl = [C1 + D12*DK*C2, D12*CK];
    Dcl = D11 + D12*DK*D21;

    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
    J_val_cur(ii) = hinfnorm(sys_cl, eps);
end

%%
Ck_vals_cur2 = linspace(-0.65,0.65,n_ax);
Bk_vals_cur2 = zeros(length(Ck_vals_cur2),1);

n_eps_cur = length(Bk_vals_cur2);

eps = 1e-9;
J_val_cur2 = zeros(n_eps_cur,1);
for ii = 1:n_eps_cur
    AK = -1;
    BK = Bk_vals_cur2(ii);
    CK = Ck_vals_cur2(ii);
    DK = DK_star;

    Acl = [A + B2*DK*C2, B2*CK;
        BK*C2, AK];
    Bcl = [B1 + B2*DK*D21; BK*D21];
    Ccl = [C1 + D12*DK*C2, D12*CK];
    Dcl = D11 + D12*DK*D21;

    sys_cl = ss(Acl, Bcl, Ccl, Dcl);
    J_val_cur2(ii) = hinfnorm(sys_cl, eps);
end
%%
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
%print(gcf,'Hinf_b.eps','-depsc2','-r300');
