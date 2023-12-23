function [argmin_Fcvx_solver] = min_Fcvx_solver(A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

% nw = nx + ny;
% nz = nx + nu;
[nx,nw]=size(B1); [nz,nx]=size(C1);

X = sdpvar(nx,nx);
Y = sdpvar(nx,nx);
M = sdpvar(nx,nx,'full');
H = sdpvar(nx,ny,'full');
F = sdpvar(nu,nx,'full');
G = sdpvar(nu,ny,'full');
gamma = sdpvar(1,1);

eps = 1e-4;

Constraints = [[X, eye(nx,nx) ; eye(nx,nx), Y] >= eps*eye(2*nx,2*nx),...
    [A*X+B2*F+(A*X+B2*F).', M.'+A+B2*G*C2, B1+B2*G*D21, (C1*X+D12*F).';...
    M+(A+B2*G*C2).', Y*A+H*C2+(Y*A+H*C2).', Y*B1+H*D21, (C1+D12*G*C2).';...
    (B1+B2*G*D21).', (Y*B1+H*D21).' -gamma*eye(nw,nw) (D11+D12*G*D21).';...
    C1*X+D12*F, C1+D12*G*C2, D11+D12*G*D21, -gamma*eye(nz,nz)] <= 0];

ops = sdpsettings('solver','mosek','verbose',1,'debug',1);
% ops = sdpsettings('solver','sdpt3','verbose',1,'debug',1);

% [model,~]=export(Constraints,1*gamma,sdpsettings('solver','mosek'));
% prob = model.prob;
% param = model.param;
% param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1.0e-20;
% param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1.0e-20;
% param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-20;
% param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = 1.0e-20;
% [r,res] = mosekopt('minimize',prob,param);

optimize(Constraints,1*gamma,ops);

% check(Constraints);

argmin_Fcvx_solver = struct;

argmin_Fcvx_solver.X = value(X);
argmin_Fcvx_solver.Y = value(Y);
argmin_Fcvx_solver.M = value(M);
argmin_Fcvx_solver.H = value(H);
argmin_Fcvx_solver.F = value(F);
argmin_Fcvx_solver.G = value(G);
argmin_Fcvx_solver.gamma = value(gamma);

end