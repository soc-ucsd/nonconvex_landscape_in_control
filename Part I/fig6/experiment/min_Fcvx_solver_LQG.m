function [argmin_Fcvx_solver] = min_Fcvx_solver_LQG(A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

X = sdpvar(nx,nx);
Y = sdpvar(nx,nx);
M = sdpvar(nx,nx,'full');
H = sdpvar(nx,ny,'full');
F = sdpvar(nu,nx,'full');
gamma = sdpvar(1,1);
Gamma = sdpvar(nx+nu,nx+nu);

eps = 1e-4;
Constraints = [[X eye(nx,nx);eye(nx,nx) Y] >= eps*eye(2*nx,2*nx),...
    [X eye(nx,nx) (C1*X+D12*F).';...
    eye(nx,nx) Y C1.';(C1*X+D12*F) C1 Gamma] >= 0,...
    [A*X+B2*F+(A*X+B2*F).' M.'+A B1 ;...
    M+A.' Y*A+H*C2+(Y*A+H*C2).' Y*B1+H*D21;...
    (B1).' (Y*B1+H*D21).' -gamma*eye(nx+ny,nx+ny)] <= 0,...
    trace(Gamma) <= gamma];

% ops = sdpsettings('solver','cdcs');
% ops.cdcs.maxIter = 10000;
% ops.cdcs.relTol  = 10^-9;
%ops.maxiter = 100000;
%ops = sdpsettings('solver','cdcs','verbose',1,'debug',1,'maxIteration',50000);

optimize(Constraints,gamma);

Xfeasible = value(X);
Yfeasible = value(Y);
Mfeasible = value(M);
Hfeasible = value(H);
Ffeasible = value(F);
gammafeasible = value(gamma);
Gammafeasible = value(Gamma);

argmin_Fcvx_solver = struct;
argmin_Fcvx_solver.X = Xfeasible;
argmin_Fcvx_solver.Y = Yfeasible;
argmin_Fcvx_solver.M = Mfeasible;
argmin_Fcvx_solver.H = Hfeasible;
argmin_Fcvx_solver.F = Ffeasible;
argmin_Fcvx_solver.gamma = gammafeasible;
argmin_Fcvx_solver.Gamma = Gammafeasible;

end