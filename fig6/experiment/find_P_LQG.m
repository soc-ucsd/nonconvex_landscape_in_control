function [P, Gamma, tmin] = find_P_LQG(AK,BK,CK,gamma,A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);

nw = nx + ny;
nz = nx + nu;

Acl = [A, B2*CK;
    BK*C2, AK];
Bcl = [B1; BK*D21];
Ccl = [C1, D12*CK];


setlmis([]);
Pvar = lmivar(1, [2*nx, 1]);
Gammavar = lmivar(1, [nx+nu, 1]);
lmiterm([1 1 1 Pvar], Acl.', 1, 's');
lmiterm([1 1 2 Pvar], 1, Bcl);
lmiterm([1 2 2 0], -gamma);

lmiterm([-2 1 1 Pvar], 1,1);
lmiterm([-2 1 2 0], Ccl.');
lmiterm([-2 2 2 Gammavar], 1, 1);

mat = eye(nx+nu);
for i = 1:nx+nu
    vec = mat(i,:);
    lmiterm([3 1 1 Gammavar],vec,vec.');
end
lmiterm([-3 1 1 0],gamma);

lmiterm([-4 1 1 Pvar],1,1);

lmisys = getlmis;

[tmin, xfeas] = feasp(lmisys);%,[0 0 -1 0 1];
P = dec2mat(lmisys, xfeas, Pvar);
Gamma = dec2mat(lmisys, xfeas, Gammavar);

end