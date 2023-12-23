function [element_Fcvx]=Phi_Hinf(K,P,gamma,A,B1,B2,C1,C2,D11,D12,D21)

[nx, nu] = size(B2);
P12 = P(1:nx,nx+1:2*nx);
P11 = P(1:nx,1:nx);
Pinv = inv(P);
Pinv21 = Pinv(nx+1:2*nx,1:nx);

AK = K.AK; BK = K.BK;
CK = K.CK; DK = K.DK;

X = Pinv(1:nx,1:nx);
Y = P(1:nx,1:nx);
M = P12*BK*C2*X+P11*B2*CK*Pinv21+P11*(A+B2*DK*C2)*X+P12*AK*Pinv21;
H = P11*B2*DK+P12*BK;
F = DK*C2*X+CK*Pinv21;
G = DK;

element_Fcvx = struct;
element_Fcvx.X = X;
element_Fcvx.Y = Y;
element_Fcvx.M = M;
element_Fcvx.H = H;
element_Fcvx.F = F;
element_Fcvx.G = G;
element_Fcvx.gamma = gamma;

end