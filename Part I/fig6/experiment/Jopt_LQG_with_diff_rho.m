function [J_opt, J_hifoo, J_mincx, J_our_lmi_solver] = Jopt_LQG_with_diff_rho(rho)
%%

A = [1 1 1;0 1 0;1 0 0];
B2 = [1 0;0 1;0 0];
C2 = [0 0 1;1 0 0;0 1 2];

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

W = 1e-1*eye(nx); V = 1e-1*eye(ny);
Q = eye(nx); R = rho*eye(nu);

B1 = [W^0.5 zeros(nx, ny)];
C1 = [Q^0.5; zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); R^0.5];
D21 = [zeros(ny, nx) V^0.5];


%% HIFOO
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

Lc = -mat( inv(kron(eye(2*nx),Acl) + ...
    kron(Acl, eye(2*nx))) * vec(Bcl*Bcl.') );
J_opt = sqrt(trace(Ccl*Lc*Ccl.'));

K_opt = struct;
K_opt.AK = AK; K_opt.BK = BK; K_opt.CK = CK;

% feasible --> K_opt nondegenerate
[P,P12,Gammaopt] = findP_LMI_LQG_SDP(AK,BK,CK,J_opt,A,B2,C2,B1,C1,D11,D12,D21);

%% solve minimization of gamma over Fcvx (LMI) (using mincx)

% return "argmin_Fcvx_mincx" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_mincx = min_Fcvx_mincx_LQG(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_mincx = argmin_Fcvx_mincx.gamma;

% [K_mincx, J_recovered_mincx] = recover_K_from_Fcvx_LQG...
%    (argmin_Fcvx_mincx, A, B2, C2, B1, C1, D11, D12, D21);

%% solve minimization of gamma over Fcvx (LMI) (using solver, e.g., Mosek)

% return "argmin_Fcvx_solver" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_solver = min_Fcvx_solver_LQG(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_solver = argmin_Fcvx_solver.gamma;

% This recovered K has large element
% also, its cost "J_recovered_solver" is lower than "J_our_lmi_solver"
[K_solver, J_recovered_solver] = recover_K_from_Fcvx_LQG...
    (argmin_Fcvx_solver, A, B2, C2, B1, C1, D11, D12, D21);

end