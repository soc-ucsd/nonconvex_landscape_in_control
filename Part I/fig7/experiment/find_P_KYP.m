function [Pfeasible,P12feasible] = find_P_KYP(K,gamma,A,B2,C2,B1,C1,D11,D12,D21)

AK = K.AK; BK = K.BK;
CK = K.CK; DK = K.DK;

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

Acl = [A + B2*DK*C2, B2*CK; BK*C2, AK];
Bcl = [B1 + B2*DK*D21; BK*D21];
Ccl = [C1 + D12*DK*C2, D12*CK];
Dcl = D11 + D12*DK*D21;

P = sdpvar(2*nx,2*nx);

Constraints =[[Acl.'*P + P*Acl, P*Bcl, Ccl.';...
              Bcl.'*P, -gamma*eye(nw), Dcl.';...
              Ccl, Dcl, -gamma*eye(nz)] <= 0,...
              P >= (1e-2)*eye(2*nx)]; % (1e-5)*eye(2*nx)
optimize(Constraints);

Pfeasible = value(P);

P12feasible = Pfeasible(1:nx,nx+1:2*nx);