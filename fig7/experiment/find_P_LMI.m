function [P, tmin] = find_P_LMI(K,hinf_norm,A,B2,C2,B1,C1,D11,D12,D21)

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
[row, col] = size(Dcl);

setlmis([]);

Pvar = lmivar(1, [2*nx, 1]);

lmiterm([1 1 1 Pvar], Acl.', 1, 's');
lmiterm([1 1 2 Pvar], 1, Bcl);
lmiterm([1 1 3 0], Ccl.');
lmiterm([1 2 2 0], -(hinf_norm)*eye(col));
lmiterm([1 2 3 0], Dcl.');
lmiterm([1 3 3 0], -(hinf_norm)*eye(row));

lmiterm([-2 1 1 Pvar],1,1);

lmisys = getlmis;

[tmin, xfeas] = feasp(lmisys, [0 0 -1 0 1]);
P = dec2mat(lmisys, xfeas, Pvar);

end