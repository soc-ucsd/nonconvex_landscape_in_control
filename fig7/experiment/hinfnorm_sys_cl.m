function J_hinf = hinfnorm_sys_cl(AK,BK,CK,DK,A,B2,C2,B1,C1,D11,D12,D21)
Acl = [A + B2*DK*C2, B2*CK;
    BK*C2, AK];
Bcl = [B1 + B2*DK*D21; BK*D21];
Ccl = [C1 + D12*DK*C2, D12*CK];
Dcl = D11 + D12*DK*D21;

sys_cl = ss(Acl, Bcl, Ccl, Dcl);
eps = 1e-9;
J_hinf = hinfnorm(sys_cl, eps);
end