function [J_hifoo, J_riccati, J_lmi, J_mincx, J_solver, J_cvx] = Jopt_Hinf_with_diff_rho(rho)
%%

A = [1 1;0 1]; B2 = [0;1]; C2 = [1 -1];

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

W = 1e-1*eye(nx);
V = 1e-1*eye(ny);
Q = eye(nx);
R = rho*eye(nu);

B1 = [W^0.5 zeros(nx, ny)];
C1 = [Q^0.5; zeros(nu, nx)];
D11 = zeros(nz, nw);
D12 = [zeros(nx, nu); R^0.5];
D21 = [zeros(ny, nx) V^0.5];

%%
argmin_Fcvx_cvx = min_Fcvx_cvx(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_cvx = argmin_Fcvx_cvx.gamma;

% [K_cvx, J_recovered_cvx] = recover_K_from_Fcvx...
%     (argmin_Fcvx_cvx, A, B2, C2, B1, C1, D11, D12, D21);
%% 1st method: hifoo

% specify system using structure
sys_struct = struct;
sys_struct.A = A; sys_struct.B1 = B1; sys_struct.B2 = B2; 
sys_struct.C1 = C1; sys_struct.C2 = C2; 
sys_struct.D11 = D11; sys_struct.D12 = D12; sys_struct.D21 = D21;

% print info of hifoo iterates
% 0 (no printing), 1 (minimal, default), 2 (includes output from HANSO), 3 (verbose)
options.prtlevel = 0;

[~, J_hifoo] = hifoo(sys_struct, nx, options); % order 2

%%
% specify system using SS object
sys_ss = ss(A, [B1 B2], [C1; C2], [D11 D12; D21 zeros(ny, nu)]);

% return a order 2 controller and its cost
[~, ~, J_riccati] = hinfsyn(sys_ss, ny, nu);

%% 3rd method: (LMI-based) hinfsyn 

% return a static (?) controller and its cost (can we check Phi(K) in Fcax?)
[~, ~, J_lmi] = hinfsyn(sys_ss, ny, nu, 'method', 'lmi');

%% 4th method: solve minimization of gamma over Fcvx (LMI) (using solver, e.g., Mosek)
% 
% % return "argmin_Fcvx_solver" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_solver = min_Fcvx_solver(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_solver = argmin_Fcvx_solver.gamma;


%% 5th method: solve minimization of gamma over Fcvx (LMI) (using mincx)

% mincx returns better cost than method 4 (other solvers)

% return "argmin_Fcvx_mincx" as a structure containing (X,Y,M,H,F,G,gamma)
argmin_Fcvx_mincx = min_Fcvx_mincx(A, B2, C2, B1, C1, D11, D12, D21);

% take out the gamma
J_mincx = argmin_Fcvx_mincx.gamma;

end