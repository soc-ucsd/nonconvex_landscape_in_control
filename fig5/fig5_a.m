close all; clear all; clc;
%%
A = -1; B2 = 1; C2 = 1;

[nx, nu] = size(B2);
[ny, ~] = size(C2);
nw = nx + ny;
nz = nx + nu;

B1 = [eye(nx) zeros(nx, ny)];
C1 = [eye(nx); zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); eye(nu)];
D21 = [zeros(ny, nx) eye(ny)];

eps_vals = linspace(-1,1,100);
Bk_vals = eps_vals;
Ck_vals = eps_vals;

n_eps = length(eps_vals);
n_ax = n_eps;
n_ay = n_eps;

[Bk_grid,Ck_grid] = meshgrid(Bk_vals, Ck_vals);

eps = 1e-9;

J_val = zeros(n_ax, n_ay);
fpeak = zeros(n_ax, n_ay);
for ii = 1:n_ax
    for jj = 1:n_ay
        AK = 0;
        BK = Bk_grid(ii, jj);
        CK = Ck_grid(ii, jj);
        DK = 0; % -0.2
        % AK = -924.2308;
        % BK = Bk_grid(ii, jj);
        % CK = Ck_grid(ii, jj);
        % DK = -0.0005;

        Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
        [J_val(ii, jj),fpeak(ii,jj)] = hinfnorm(sys_cl, eps);
    end
end
%%
figure;
surf(Bk_grid, Ck_grid, J_val);
colorbar;
xlabel('$B_k$','Interpreter','latex'); ylabel('$C_k$','Interpreter','latex');
title('$J(K)$','Interpreter','latex');
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [7, 6];
FontSize = 11;

surf(Bk_grid, Ck_grid, J_val,'EdgeColor', 'none','FaceAlpha',0.9); 
hold on;

xlabel('$B_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
ylabel('$C_{\mathsf{K}}$','Interpreter','latex','FontSize',FontSize);
zlabel('$J_{\infty,1}(\mathsf{K})$','Interpreter','latex','FontSize',FontSize);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

view(-80,25)
print(gcf,'fig5_a.eps','-depsc2','-r300');