function tf_Fcvx = is_in_Fcvx_hinf(element_Fcvx,A,B2,C2,B1,C1,D11,D12,D21)

[nx, nu] = size(B2);
[ny, ~] = size(C2);
nw = nx + ny;
nz = nx + nu;

X = element_Fcvx.X;
Y = element_Fcvx.Y;
M = element_Fcvx.M;
H = element_Fcvx.H;
F = element_Fcvx.F;
G = element_Fcvx.G;
gamma = element_Fcvx.gamma; 

AA = [X eye(nx,nx) ; eye(nx,nx) Y];
BB = [A*X+B2*F+(A*X+B2*F).' M.'+A+B2*G*C2 B1+B2*G*D21 (C1*X+D12*F).';...
    M+(A+B2*G*C2).' Y*A+H*C2+(Y*A+H*C2).' Y*B1+H*D21 (C1+D12*G*C2).';...
    (B1+B2*G*D21).' (Y*B1+H*D21).' -gamma*eye(nx+ny,nx+ny) (D11+D12*G*D21).';...
    C1*X+D12*F C1+D12*G*C2 D11+D12*G*D21 -gamma*eye(nz,nz)];

% check if AA > 0 and BB <= 0
% if min(eig(AA)) <= -1e-7 || max(eig(BB)) >= 1e-6
%     tf_Fcvx = false;
% else
%     tf_Fcvx = true;
% end
if min(eig(AA)) > 0 && max(eig(BB)) <= 1e-6
    tf_Fcvx = true;
else
    tf_Fcvx = false;
end