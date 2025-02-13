%% Compare 5 methods of hinf

% J_hifoo < J_hinfsyn_lmi ~ J_our_lmi_mincx < J_hinfsyn_riccati <
% J_our_lmi_solver

clear all; close all; clc;

%% problem data
A = [1 1 1;0 1 0;1 0 0];
B2 = [1 0;0 1;0 0];
C2 = [0 0 1;1 0 0;0 1 2];

[nx, nu] = size(B2);
[ny, ~] = size(C2);

% nw = nx + ny;
% nz = nx + nu;

B1 = [eye(nx) zeros(nx, ny)];
C1 = [eye(nx); zeros(nu, nx)];

D12 = [zeros(nx, nu); eye(nu)];
D21 = [zeros(ny, nx) eye(ny)];


[nx, nw]=size(B1); [nz, nx]=size(C1);

D11 = zeros(nz, nw);

%%
argmin_Fcvx_cvx = min_Fcvx_cvx(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_cvx = argmin_Fcvx_cvx.gamma;

[K_cvx, J_recovered_cvx] = recover_K_from_Fcvx...
    (argmin_Fcvx_cvx, A, B2, C2, B1, C1, D11, D12, D21);
%% 1st method: hifoo

% load the saved K_hifoo and J_hifoo from hifoo
% load("hifoo_3.mat"); 

% specify system using structure
sys_struct = struct;
sys_struct.A = A; sys_struct.B1 = B1; sys_struct.B2 = B2; 
sys_struct.C1 = C1; sys_struct.C2 = C2; 
sys_struct.D11 = D11; sys_struct.D12 = D12; sys_struct.D21 = D21;

% print info of hifoo iterates
% 0 (no printing), 1 (minimal, default), 2 (includes output from HANSO), 3 (verbose)
options.prtlevel = 3;

[K_Hifoo, J_hifoo] = hifoo(sys_struct, 3, options); % order 3

% save the controller as a structure variable
K_hifoo = struct;
K_hifoo.AK = K_Hifoo.a;
K_hifoo.BK = K_Hifoo.b;
K_hifoo.CK = K_Hifoo.c;
K_hifoo.DK = K_Hifoo.d;

%% 2nd method: (default Riccati-based) hinfsyn (from Matlab robust control toolbox)

% specify system using SS object
sys_ss = ss(A, [B1 B2], [C1; C2], [D11 D12; D21 zeros(ny, nu)]);

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
argmin_Fcvx_solver = min_Fcvx_solver(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_solver = argmin_Fcvx_solver.gamma;

% This recovered K has large element
% also, its cost "J_recovered_solver" is lower than "J_our_lmi_solver"
[K_solver, J_recovered_solver] = recover_K_from_Fcvx...
    (argmin_Fcvx_solver, A, B2, C2, B1, C1, D11, D12, D21);

%% 5th method: solve minimization of gamma over Fcvx (LMI) (using mincx)

% mincx returns better cost than method 4 (other solvers)

% return "argmin_Fcvx_mincx" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_mincx = min_Fcvx_mincx(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_mincx = argmin_Fcvx_mincx.gamma;

% This recovered K also has large element
% but now "J_recovered_mincx" is equal to "J_our_lmi_mincx"
[K_mincx, J_recovered_mincx] = recover_K_from_Fcvx...
    (argmin_Fcvx_mincx, A, B2, C2, B1, C1, D11, D12, D21);