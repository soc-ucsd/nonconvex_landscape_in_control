clear all; close all; clc;
%%
A = [1 1;0 1]; B2 = [0;1]; C2 = [1 1;1 0];

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

B1 = [eye(nx) zeros(nx, ny)];
C1 = [eye(nx); zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); eye(nu)];
D21 = [zeros(ny, nx) eye(ny)];
%% min Fcvx (LMI)
[Fcvx_min] = min_Fcvx_hinf_solver(A,B2,C2,B1,C1,D11,D12,D21);
%%
% returned point by min Fcvx
X = Fcvx_min.X;
Y = Fcvx_min.Y;
M = Fcvx_min.M;
H = Fcvx_min.H;
F = Fcvx_min.F;
G = Fcvx_min.G;
gamma = Fcvx_min.gamma;

% check if feasible (yes)
% tf_Fcvx_min = isFeasible_Fcvx_Hinf(Fcvx_min,A,B2,C2,B1,C1,D11,D12,D21);
AA = [X eye(nx,nx) ; eye(nx,nx) Y];
BB = [A*X+B2*F+(A*X+B2*F).' M.'+A+B2*G*C2 B1+B2*G*D21 (C1*X+D12*F).';...
    M+(A+B2*G*C2).' Y*A+H*C2+(Y*A+H*C2).' Y*B1+H*D21 (C1+D12*G*C2).';...
    (B1+B2*G*D21).' (Y*B1+H*D21).' -gamma*eye(nx+ny,nx+ny) (D11+D12*G*D21).';...
    C1*X+D12*F C1+D12*G*C2 D11+D12*G*D21 -gamma*eye(nz,nz)];
%%
% reduce a bit J and see if still feasible (no)
gamma_ = gamma - 0.001;
AAA = [X eye(nx,nx) ; eye(nx,nx) Y];
BBB = [A*X+B2*F+(A*X+B2*F).' M.'+A+B2*G*C2 B1+B2*G*D21 (C1*X+D12*F).';...
    M+(A+B2*G*C2).' Y*A+H*C2+(Y*A+H*C2).' Y*B1+H*D21 (C1+D12*G*C2).';...
    (B1+B2*G*D21).' (Y*B1+H*D21).' -gamma_*eye(nx+ny,nx+ny) (D11+D12*G*D21).';...
    C1*X+D12*F C1+D12*G*C2 D11+D12*G*D21 -gamma_*eye(nz,nz)];
%%
% recover controller K 
% [K_LMI_Fcvx,J_LMI_Fcvx] = recover_K_from_Fcvx_Hinf(Fcvx_min,A,B2,C2,B1,C1,D11,D12,D21);

E = eye(nx,nx); % P12_min;  % any invertible
PI = -inv(E)*(Y-inv(X))*X;

K_mat = ([eye(nu,nu) zeros(nu,nx); Y*B2 E])\...
    [G F; H M-Y*A*X]/...
    ([eye(ny,ny) C2*X; zeros(nx,ny), PI]);

K = struct;
K.AK = K_mat(nu+1:nu+nx,ny+1:ny+nx);
K.BK = K_mat(nu+1:nu+nx,1:ny);
K.CK = K_mat(1:nu,ny+1:ny+nx);
K.DK = K_mat(1:nu,1:ny);

J = hinfnorm_sys_cl(K.AK,K.BK,K.CK,K.DK,A,B2,C2,B1,C1,D11,D12,D21);

% check if K nondegenerate (yes)
[P_min,P12_min] = findP_LMI_Hinf_SDP(K,J,A,B2,C2,B1,C1,D11,D12,D21);

% get Phi(K,J,P_min)
Fcvx_min_2 = Phi_Hinf(K,P_min,J,A,B1,B2,C1,C2,D11,D12,D21);

% another feasible point in Fcvx with smaller cost J_LMI_Fcvx
tf_Fcvx_min_2 = isFeasible_Fcvx_Hinf(Fcvx_min_2,A,B2,C2,B1,C1,D11,D12,D21);

% check if (X,Y,M,H,F,G,gamma) feasible (yes, but with smaller cost J)
X2 = Fcvx_min_2.X;
Y2 = Fcvx_min_2.Y;
M2 = Fcvx_min_2.M;
H2 = Fcvx_min_2.H;
F2 = Fcvx_min_2.F;
G2 = Fcvx_min_2.G;
gamma2 = Fcvx_min_2.gamma;

AAAA = [X2 eye(nx,nx) ; eye(nx,nx) Y2];
BBBB = [A*X2+B2*F2+(A*X2+B2*F2).' M2.'+A+B2*G2*C2 B1+B2*G2*D21 (C1*X+D12*F2).';...
    M2+(A+B2*G2*C2).' Y2*A+H2*C2+(Y2*A+H2*C2).' Y2*B1+H2*D21 (C1+D12*G2*C2).';...
    (B1+B2*G2*D21).' (Y2*B1+H2*D21).' -J*eye(nx+ny,nx+ny) (D11+D12*G2*D21).';...
    C1*X2+D12*F2 C1+D12*G2*C2 D11+D12*G2*D21 -J*eye(nz,nz)];
%%
J_min = Fcvx_min.gamma;
J_recover_min = Fcvx_min_2.gamma;
%%
E = P12_min; % any invertible
PI = -inv(E)*(Y2-inv(X2))*X2;

K2_mat = ([eye(nu,nu) zeros(nu,nx); Y2*B2 E])\...
    [G2 F2; H2 M2-Y2*A*X2]/...
    ([eye(ny,ny) C2*X2; zeros(nx,ny), PI]);

K2 = struct;
K2.AK = K2_mat(nu+1:nu+nx,ny+1:ny+nx);
K2.BK = K2_mat(nu+1:nu+nx,1:ny);
K2.CK = K2_mat(1:nu,ny+1:ny+nx);
K2.DK = K2_mat(1:nu,1:ny);

JJ = hinfnorm_sys_cl(K2.AK,K2.BK,K2.CK,K2.DK,A,B2,C2,B1,C1,D11,D12,D21);

%%
Fcvx_min_3 = Phi_Hinf(K2,P_min,JJ,A,B1,B2,C1,C2,D11,D12,D21);

X3 = Fcvx_min_3.X;
Y3 = Fcvx_min_3.Y;
M3 = Fcvx_min_3.M;
H3 = Fcvx_min_3.H;
F3 = Fcvx_min_3.F;
G3 = Fcvx_min_3.G;
gamma3 = Fcvx_min_3.gamma;