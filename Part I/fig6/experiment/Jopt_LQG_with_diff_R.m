function [J_opt, J_hifoo, J1, J_our_lmi_mincx, J_our_lmi_solver] = Jopt_LQG_with_diff_R(R)
%%

A = [1 1;0 1]; B2 = [0;1]; C2 = [1 -1];

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

W = 1e-3*eye(2);
V = 1e-3;
Q = eye(2);
% R = 1:10;

B1 = [W^0.5 zeros(nx, ny)];
C1 = [Q^0.5; zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); R^0.5];
D21 = [zeros(ny, nx) V^0.5];

%%
B = [0;1]; C = [1 -1]; Qc = eye(2);
Num = 1; % number of initial points
% ---------------------------------------------
% Gradient descent 
% ---------------------------------------------
info_full = cell(Num,1);
info_cano = cell(Num,1);
K_full    = cell(Num,1);  % controller
K_cano    = cell(Num,1);  % controller
J1 = 0;
% for ind = 1:Num
%     p = rand(nx,1)-2;   % a random number between (-2,-1)
%     K = place(A,B,p);
%     L = place(A',C',p); L = L';
%     Ak = A - B*K- L*C;
%     Bk = L;
%     Ck = -K;
%     K0.Ak = Ak; K0.Bk = Bk;K0.Ck = Ck;
% 
%     opts.tol      = 1e-3;
%     opts.maxIter  = 5e3;
%     opts.Disp     = 100;
% 
%     % full gradient
%     opts.stepsize = 5e-5;
%     [K1,J1,info1] = LQG_gd_full(A,B,C,Qc,R,W,V,K0,opts);
%     info_full{ind} = info1;
%     K_full{ind}    = K1;
%     % Hess1 = LQGhessfull(A,B,C,K1,Qc,R,W,V);  % hessian
% 
%     % gradient over canonical form
%     % opts.stepsize = 5e-2;
%     % [K2,J2,info2] = LQG_gd_cano(A,B,C,Qc,R,W,V,K0,opts);
%     % info_cano{ind} = info2;
%     % K_cano{ind}    = K2;
%     % Hess2 = LQGhessfull(A,B,C,K2,Qc,R,W,V);  % hessian
% end
%%
sys_struct = struct;
sys_struct.A = A; sys_struct.B1 = B1; sys_struct.B2 = B2; 
sys_struct.C1 = C1; sys_struct.C2 = C2; 
sys_struct.D11 = D11; sys_struct.D12 = D12; sys_struct.D21 = D21;

[K_hifoo, J_hifoo] = hifoo(sys_struct, nx, 't');

%% optimal controller from matlab, compute optimal cost using formula
sys_pl = ss(A, B2, C2, 0);
QXU = blkdiag(Q,R);
QWV = blkdiag(W,V);
[reg,info] = lqg(sys_pl, QXU, QWV);

AK = reg.A; BK = reg.B; CK = reg.C;

Acl = [A, B2*CK; BK*C2, AK];
Bcl = [B1; BK*D21];
Ccl = [C1, D12*CK];
Dcl = D11;

Lc = -mat(inv(kron(eye(4),Acl) + kron(Acl,eye(4)))*vec(Bcl*Bcl.'));
J_opt = sqrt(trace(Ccl*Lc*Ccl.'));

K_opt = struct;
Kopt.AK = AK; Kopt.BK = BK; Kopt.CK = CK;

% feasible --> K_opt nondegenerate
% [P,P12,Gammaopt] = findP_LMI_LQG_SDP(AK,BK,CK,J_opt,A,B2,C2,B1,C1,D11,D12,D21);

%% solve minimization of gamma over Fcvx (LMI) (using mincx)

% return "argmin_Fcvx_mincx" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_mincx = min_Fcvx_mincx_LQG(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_mincx = argmin_Fcvx_mincx.gamma;

% [K_mincx, J_recovered_mincx] = recover_K_from_Fcvx_LQG...
%     (argmin_Fcvx_mincx, A, B2, C2, B1, C1, D11, D12, D21);

%% solve minimization of gamma over Fcvx (LMI) (using solver, e.g., Mosek)

% return "argmin_Fcvx_solver" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_solver = min_Fcvx_solver_LQG(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_solver = argmin_Fcvx_solver.gamma;

% This recovered K has large element
% also, its cost "J_recovered_solver" is lower than "J_our_lmi_solver"
% [K_solver, J_recovered_solver] = recover_K_from_Fcvx_LQG...
%     (argmin_Fcvx_solver, A, B2, C2, B1, C1, D11, D12, D21);

end