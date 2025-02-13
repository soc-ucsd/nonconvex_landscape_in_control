function [K,J] = recover_K_from_Fcvx_LQG(element_Fcvx,A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

X = element_Fcvx.X;
Y = element_Fcvx.Y;
M = element_Fcvx.M;
H = element_Fcvx.H;
F = element_Fcvx.F;

E = eye(nx,nx); % any invertible
PI = -inv(E)*(Y-inv(X))*X;
K_mat = inv([eye(nu,nu) zeros(nu,nx); Y*B2 E])*...
    [zeros(nu,ny) F; H M-Y*A*X]*...
    inv([eye(ny,ny) C2*X; zeros(nx,ny), PI]);

K = struct;
K.AK = K_mat(nu+1:nu+nx,ny+1:ny+nx);
K.BK = K_mat(nu+1:nu+nx,1:ny);
K.CK = K_mat(1:nu,ny+1:ny+nx);
K.DK = K_mat(1:nu,1:ny);

% compute J_lqg
AK = K.AK; BK = K.BK;
CK = K.CK; DK = K.DK;

Acl = [A + B2*DK*C2, B2*CK;
    BK*C2, AK];
Bcl = [B1 + B2*DK*D21; BK*D21];
Ccl = [C1 + D12*DK*C2, D12*CK];
Dcl = D11 + D12*DK*D21;
sys_cl = ss(Acl, Bcl, Ccl, Dcl);

J = norm(sys_cl, 2);