function [K,J] = recover_K_from_Fcvx(element_Fcvx,A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

X = element_Fcvx.X;
Y = element_Fcvx.Y;
M = element_Fcvx.M;
H = element_Fcvx.H;
F = element_Fcvx.F;
G = element_Fcvx.G;
gamma = element_Fcvx.gamma; 

E = eye(nx,nx); % any invertible
PI = -inv(E)*(Y-inv(X))*X;
K_mat = inv([eye(nu,nu) zeros(nu,nx); Y*B2 E])*...
    [G F; H M-Y*A*X]*...
    inv([eye(ny,ny) C2*X; zeros(nx,ny), PI]);

K = struct;
K.AK = K_mat(nu+1:nu+nx,ny+1:ny+nx);
K.BK = K_mat(nu+1:nu+nx,1:ny);
K.CK = K_mat(1:nu,ny+1:ny+nx);
K.DK = K_mat(1:nu,1:ny);

J = hinfnorm_sys_cl(K.AK,K.BK,K.CK,K.DK,A,B2,C2,B1,C1,D11,D12,D21);