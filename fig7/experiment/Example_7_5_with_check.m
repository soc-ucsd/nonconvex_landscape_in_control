%% Compare 5 methods of hinf

% Compared to "Compare_hinf_sol_dim_2.m", this file checks the non-degeneracy or feasibility of the
% of the returned solutions

clear all; close all; clc;

%% problem data
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

%% 1st method: hifoo

% load the saved K_hifoo and J_hifoo from hifoo
load("hifoo_2.mat"); 

% % specify system using structure
% sys_struct = struct;
% sys_struct.A = A; sys_struct.B1 = B1; sys_struct.B2 = B2; 
% sys_struct.C1 = C1; sys_struct.C2 = C2; 
% sys_struct.D11 = D11; sys_struct.D12 = D12; sys_struct.D21 = D21;
% 
% % print info of hifoo iterates
% % 0 (no printing), 1 (minimal, default), 2 (includes output from HANSO), 3 (verbose)
% options.prtlevel = 0;
% 
% [K_Hifoo, J_hifoo] = hifoo(sys_struct, 2, options); % order 2
% 
% % save the controller as a structure variable
% K_hifoo = struct;
% K_hifoo.AK = K_Hifoo.a;
% K_hifoo.BK = K_Hifoo.b;
% K_hifoo.CK = K_Hifoo.c;
% K_hifoo.DK = K_Hifoo.d;

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

[P, tmin] = find_P_LMI(K_riccati,J_hinfsyn_riccati,A,B2,C2,B1,C1,D11,D12,D21);
%% 3rd method: (LMI-based) hinfsyn 

% return a static (?) controller and its cost (can we check Phi(K) in Fcax?)
[K_hinfsyn_lmi, sys_cl_hinfsyn_lmi, J_hinfsyn_lmi] = hinfsyn(sys_ss, ny, nu, 'method', 'lmi');

%% 4th method: solve minimization of gamma over Fcvx (LMI) (using solver, e.g., Mosek)

% return "argmin_Fcvx_solver" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_solver = min_Fcvx_solver(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_solver = argmin_Fcvx_solver.gamma;

% This recovered K has very large element
% also, its cost "J_recovered_solver" is lower than "J_our_lmi_solver"
[K_solver, J_recovered_solver] = recover_K_from_Fcvx...
    (argmin_Fcvx_solver, A, B2, C2, B1, C1, D11, D12, D21);

%[P, tmin] = find_P_LMI(K_solver,J_recovered_solver,A,B2,C2,B1,C1,D11,D12,D21);
%% 5th method: solve minimization of gamma over Fcvx (LMI) (using mincx)
% mincx returns better cost than other solvers

% return "argmin_Fcvx_mincx" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_mincx = min_Fcvx_mincx(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_our_lmi_mincx = argmin_Fcvx_mincx.gamma;

% This recovered K has very large element
% but now "J_recovered_mincx" is equal to "J_our_lmi_mincx"
[K_mincx, J_recovered_mincx] = recover_K_from_Fcvx...
    (argmin_Fcvx_mincx, A, B2, C2, B1, C1, D11, D12, D21);

%[P, tmin] = find_P_LMI(K_mincx,J_our_lmi_mincx,A,B2,C2,B1,C1,D11,D12,D21);
%% verify if "Phi(K_hifoo, P_hifoo, J_hifoo)" is feasible in Fcvx

% find the certificate P for hinfnorm evaluation of K
% infeasible
[P_hifoo, P12_hifoo] = ...
    find_P_KYP(K_hifoo, J_hifoo, A, B2, C2, B1, C1, D11, D12, D21);

% we can make it feasible by perturbing(increasing) J_hifoo a little bit
[P_hifoo_pert, P12_hifoo_pert] = ...
    find_P_KYP(K_hifoo, J_hifoo + 1e-2, A, B2, C2, B1, C1, D11, D12, D21);

% seems these are due to numerical issues when solving LMI?
% Can we conclude K_hifoo is degenerate

% Find P by LMI-based approach
[P, tmin] = find_P_LMI(K_hifoo,J_hifoo,A,B2,C2,B1,C1,D11,D12,D21);

Acl = [A + B2*K_hifoo.DK*C2, B2*K_hifoo.CK; K_hifoo.BK*C2, K_hifoo.AK];
Bcl = [B1 + B2*K_hifoo.DK*D21; K_hifoo.BK*D21];
Ccl = [C1 + D12*K_hifoo.DK*C2, D12*K_hifoo.CK];
Dcl = D11 + D12*K_hifoo.DK*D21;

test_M = [Acl'*P + P*Acl, P*Bcl, Ccl';
    Bcl'*P, -J_hifoo*eye(nw), Dcl';
    Ccl, Dcl, -J_hifoo*eye(nz)];

% minP < 1e-4 --> P may not be positive definite
% maxM >= 1e-6 --> may not exist P to certify
minP = min(eig(P));
maxM = max(eig(test_M));


%% verify if "K_riccati after Phi" is feasible in Fcvx

% feasible, K_riccati non-degenerate
[P_riccati, P12_riccati] = find_P_KYP...
    (K_riccati, J_hinfsyn_riccati, A, B2, C2, B1, C1, D11, D12, D21);

Phi_riccati = Phi_Hinf...
    (K_riccati, P_riccati, J_hinfsyn_riccati, A, B1, B2, C1, C2, D11, D12, D21);

is_Phi_riccati_in_Fcvx = is_in_Fcvx_hinf(Phi_riccati, A, B2, C2, B1, C1, D11, D12, D21);

%% verify if "argmin_Fcvx_solver" (returned by method 4) is feasible in Fcvx

% feasible
is_argmin_Fcvx_solver_in_Fcvx = is_in_Fcvx_hinf...
    (argmin_Fcvx_solver, A, B2, C2, B1, C1, D11, D12, D21);

% reduce the cost a bit
argmin_Fcvx_solver_copy = argmin_Fcvx_solver;
argmin_Fcvx_solver_copy.gamma = J_our_lmi_solver - 1e-2;

% argmin_Fcvx_solver_copy (with reduced gamma) becomes infeasible
is_argmin_Fcvx_solver_copy_in_Fcvx = is_in_Fcvx_hinf...
    (argmin_Fcvx_solver_copy, A, B2, C2, B1, C1, D11, D12, D21);

% check if K_solver nondegenerate
[P_solver, P12_solver] = find_P_KYP(K_solver, J_recovered_solver, A, B2, C2, B1, C1, D11, D12, D21);

% get Phi(K_solver, J_our_lmi_solver, P_solver)
Phi_solver = Phi_Hinf(K_solver, P_solver, J_recovered_solver, A, B1, B2, C1, C2, D11, D12, D21);

% another feasible point in Fcvx with a smaller cost J_recovered_solver
is_Phi_solver_in_Fcvx = is_in_Fcvx_hinf(Phi_solver, A, B2, C2, B1, C1, D11, D12, D21);