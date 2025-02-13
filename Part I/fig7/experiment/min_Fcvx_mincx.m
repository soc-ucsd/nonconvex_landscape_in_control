function [argmin_Fcvx_mincx] = min_Fcvx_mincx(A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

setlmis([]);

Xvar = lmivar(1,[nx,1]);
Yvar = lmivar(1,[nx,1]);
Mvar = lmivar(2,[nx,nx]);
Hvar = lmivar(2,[nx,ny]);
Fvar = lmivar(2,[nu,nx]);
Gvar = lmivar(2,[nu,ny]);
gammavar = lmivar(1,[1,1]);

%eps = 1e-3;
lmiterm([-1 1 1 Xvar],1,1);
lmiterm([-1 1 2 0],1);
lmiterm([-1 2 2 Yvar],1,1);
%lmiterm([1 1 1 0],eps);

lmiterm([2 1 1 Xvar],A,1,'s');
lmiterm([2 1 1 Fvar],B2,1,'s');
lmiterm([2 1 2 -Mvar],1,1);
lmiterm([2 1 2 0],A);
lmiterm([2 1 2 Gvar],B2,C2);
lmiterm([2 1 3 0],B1);
lmiterm([2 1 3 Gvar],B2,D21);
lmiterm([2 1 4 -Xvar],1,C1.');
lmiterm([2 1 4 -Fvar],1,D12.');

lmiterm([2 2 2 Yvar],1,A,'s');
lmiterm([2 2 2 Hvar],1,C2,'s');
lmiterm([2 2 3 Yvar],1,B1);
lmiterm([2 2 3 Hvar],1,D21);
lmiterm([2 2 4 0],C1.');
lmiterm([2 2 4 -Gvar],C2.',D12.');

lmiterm([2 3 3 gammavar],-1,1);
lmiterm([2 3 4 0],D11.');
lmiterm([2 3 4 -Gvar],D21.',D12.');

lmiterm([2 4 4 gammavar],-1,1);

lmisys = getlmis;

% lmiinfo(lmisys);

options = [1e-12, 1e5, 0, 0, 0];

[copt, xopt] = mincx(lmisys, [zeros(1, decnbr(lmisys)-1), 1], options);

argmin_Fcvx_mincx = struct;

argmin_Fcvx_mincx.X = dec2mat(lmisys, xopt, Xvar);
argmin_Fcvx_mincx.Y = dec2mat(lmisys, xopt, Yvar);
argmin_Fcvx_mincx.M = dec2mat(lmisys, xopt, Mvar);
argmin_Fcvx_mincx.H = dec2mat(lmisys, xopt, Hvar);
argmin_Fcvx_mincx.F = dec2mat(lmisys, xopt, Fvar);
argmin_Fcvx_mincx.G = dec2mat(lmisys, xopt, Gvar);
argmin_Fcvx_mincx.gamma = dec2mat(lmisys, xopt, gammavar);

end