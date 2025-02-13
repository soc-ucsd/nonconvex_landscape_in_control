function [P] = findP_Riccati(A, B, C, D, gamma)
[ny, nu] = size(D);

B_scl = B/sqrt(gamma);
C_scl = C/sqrt(gamma);
D_scl = D/gamma;

if norm(D_scl'*D_scl) >= 1
    P = [];
    return;
end

AA = A + (B_scl / (eye(nu) - D_scl'*D_scl)) * D_scl'*C_scl;
GG = (B_scl / (eye(nu) - D_scl'*D_scl)) * B_scl';
GG = (GG + GG')/2;
QQ = C_scl'*C_scl + C_scl'* (D_scl/(eye(nu) - D_scl'*D_scl)) * D_scl' * C_scl;
QQ = (QQ + QQ')/2;

P = icare(AA, [], QQ, [], [], [], GG);
end