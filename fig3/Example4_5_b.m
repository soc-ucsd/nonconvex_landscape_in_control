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

t1 = linspace(-0.1,0.1,100);
t2 = linspace(-0.1,0.1,100);
[t1_grid,t2_grid] = meshgrid(t1, t2);

J_val = zeros(length(t1),length(t2));
for ii = 1 : length(t1)
    for jj = 1 : length(t2)
        AK = -1;
        BK = t1_grid(ii,jj)-t2_grid(ii,jj);
        CK = t1_grid(ii,jj)+t2_grid(ii, jj);
        DK = 0;

        Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
        J_val(ii,jj) = norm(sys_cl, 2)^2;
    end
end

%% one point 
AK = -1; BK = 0; CK = 0;  DK = 0;
Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
J_center = norm(sys_cl, 2);


%% two line
t1 = linspace(-0.1,0.1,100);
t2 = t1;

J_line1 = zeros(length(t1),1);
for ii = 1 : length(BK_1)
        AK = -1;
        BK = t1_grid(ii,jj)-t2_grid(ii,jj);
        CK = t1_grid(ii,jj)+t2_grid(ii, jj);
        DK = 0;

        Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
        J_line1(ii) = norm(sys_cl, 2);
end

t11 = linspace(-0.1,0.1,100);
t2 = t1;

J_line2 = zeros(length(t11),1);
for ii = 1 : length(BK_1)
        AK = -1;
        BK = t1_grid(ii,jj)-t2_grid(ii,jj);
        CK = t1_grid(ii,jj)+t2_grid(ii, jj);
        DK = 0;

        Acl = [A + B2*DK*C2, B2*CK;
            BK*C2, AK];
        Bcl = [B1 + B2*DK*D21; BK*D21];
        Ccl = [C1 + D12*DK*C2, D12*CK];
        Dcl = D11 + D12*DK*D21;

        sys_cl = ss(Acl, Bcl, Ccl, Dcl);
        J_line2(ii) = norm(sys_cl, 2);
end
save LQGdata2