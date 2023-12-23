clear all; close all; clc;
%%
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib('REA2');

%%
argmin_Fcvx_cvx = min_Fcvx_cvx(A, B, C, B1, C1, D11, D12, D21);

% take out the gamma
J_cvx = argmin_Fcvx_cvx.gamma;

[K_cvx, J_recovered_cvx] = recover_K_from_Fcvx...
    (argmin_Fcvx_cvx, A, B, C, B1, C1, D11, D12, D21);
%%
% specify system using structure
sys_struct = struct;
sys_struct.A = A; sys_struct.B1 = B1; sys_struct.B2 = B; 
sys_struct.C1 = C1; sys_struct.C2 = C; 
sys_struct.D11 = D11; sys_struct.D12 = D12; sys_struct.D21 = D21;

% print info of hifoo iterates
% 0 (no printing), 1 (minimal, default), 2 (includes output from HANSO), 3 (verbose)
options.prtlevel = 3;

[K_Hifoo, J_hifoo] = hifoo(sys_struct, nx, options); % order nx

% save the controller as a structure variable
K_hifoo = struct;
K_hifoo.AK = K_Hifoo.a;
K_hifoo.BK = K_Hifoo.b;
K_hifoo.CK = K_Hifoo.c;
K_hifoo.DK = K_Hifoo.d;

%%
argmin_Fcvx_mincx = min_Fcvx_mincx(A, B, C, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_mincx = argmin_Fcvx_mincx.gamma;

% This recovered K also has large element
% but now "J_recovered_mincx" is equal to "J_our_lmi_mincx"
[K_mincx, J_recovered_mincx] = recover_K_from_Fcvx...
    (argmin_Fcvx_mincx, A, B, C, B1, C1, D11, D12, D21);

%% 2nd method: (default Riccati-based) hinfsyn (from Matlab robust control toolbox)

% specify system using SS object
sys_ss = ss(A, [B1 B], [C1; C], [D11 D12; D21 zeros(ny, nu)]);

% return a order 2 controller and its cost
[K_hinfsyn_riccati, sys_cl_riccati, J_hinfsyn_riccati] = hinfsyn(sys_ss, ny, nu);

% save the controller as a structure variable
K_riccati = struct;
K_riccati.AK = K_hinfsyn_riccati.A;
K_riccati.BK = K_hinfsyn_riccati.B;
K_riccati.CK = K_hinfsyn_riccati.C;
K_riccati.DK = K_hinfsyn_riccati.D;

%% 3rd method: (LMI-based) hinfsyn 

% return a static (?) controller and its cost (can we check Phi(K) in Fcax?)
[K_hinfsyn_lmi, sys_cl_hinfsyn_lmi, J_hinfsyn_lmi] = hinfsyn(sys_ss, ny, nu, 'method', 'lmi');

%% 4th method: solve minimization of gamma over Fcvx (LMI) (using solver, e.g., Mosek)

% return "argmin_Fcvx_solver" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_solver = min_Fcvx_solver(A, B, C, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_solver = argmin_Fcvx_solver.gamma;

% This recovered K has large element
% also, its cost "J_recovered_solver" is lower than "J_our_lmi_solver"
[K_solver, J_recovered_solver] = recover_K_from_Fcvx...
    (argmin_Fcvx_solver, A, B, C, B1, C1, D11, D12, D21);