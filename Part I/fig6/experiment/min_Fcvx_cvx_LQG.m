function [argmin_Fcvx_solver] = min_Fcvx_cvx_LQG(A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

% nw = nx + ny;
% nz = nx + nu;
[nx,nw]=size(B1); [nz,nx]=size(C1);

% cvx_solver mosek
% cvx_solver SeDuMi
% cvx_solver SDPT3

cvx_begin sdp
    variable X(nx,nx) symmetric
    variable Y(nx,nx) symmetric
    variable M(nx,nx)
    variable H(nx,ny)
    variable F(nu,nx)
    variable Gamma(nz,nz)
    variable gam
    minimize (gam)
    [X eye(nx); eye(nx) Y] >=  1e-4*eye(2*nx)
    [X eye(nx,nx) (C1*X+D12*F).';...
    eye(nx,nx) Y C1.';(C1*X+D12*F) C1 Gamma] >= 0
    [A*X+B2*F+(A*X+B2*F).' M.'+A B1 ;...
    M+A.' Y*A+H*C2+(Y*A+H*C2).' Y*B1+H*D21;...
    (B1).' (Y*B1+H*D21).' -gam*eye(nw,nw)] <= 0
    trace(Gamma) <= gam
cvx_end

argmin_Fcvx_solver = struct;

argmin_Fcvx_solver.X = X;
argmin_Fcvx_solver.Y = Y;
argmin_Fcvx_solver.M = M;
argmin_Fcvx_solver.H = H;
argmin_Fcvx_solver.F = F;
argmin_Fcvx_solver.Gamma = Gamma;
argmin_Fcvx_solver.gamma = gam;

end