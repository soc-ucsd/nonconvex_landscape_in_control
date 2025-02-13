
% Fig 1(b) in "Benign Nonconvex Landscapes in Optimal and Robust Control,
%           Part II: Extended Convex Lifting"
%       By Yang Zheng, Chih-Fan Pai, Yujie Tang


clear all;close all;clc;
%%
A=[0 0;0 0];
[nx,~]=size(A);
B=eye(2);
Bw=eye(2);
Omega=Bw*Bw.';%[1 1;1 1];
Q=eye(2);
R=eye(2);

PKstar = icare(A,[],Q,[],[],[],-B*R^(-1)*B.');
Kstar = -R^(-1)*B.'*PKstar;
JKstar = trace(PKstar*Omega);

%% gradient
XKstar=lyap(A+B*Kstar,Omega); %-mat(inv(kron(eye(nx),A+B*Kstar) + kron(A+B*Kstar,eye(nx)))*vec(Omega));
gradstar=2*(R*Kstar+B.'*PKstar)*XKstar;
PKstar1 = lyap((A+B*Kstar).',Q+Kstar.'*R*Kstar); %-mat(inv(kron(eye(nx),(A+B*Kstar).') + kron((A+B*Kstar).',eye(nx)))*vec(Q+Kstar.'*R*Kstar));
gradstar1=2*(R*Kstar+B.'*PKstar1)*XKstar;
%%
n1 = 500; n2 = 500; r = 15;
% k1_range = linspace(Kstar(1)-r, Kstar(1)+1, n1);
% k2_range = linspace(Kstar(2)-1, Kstar(2)+r, n2);
k1_range = linspace(-4,4, n1);
k2_range = linspace(-4,4, n2);

Jlqr = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        K = [Kstar(1,1), k1_range(i);k2_range(j),Kstar(2,2)];

        if max(real(eig(A+B*K)))>=-1e-4
            Jlqr(i,j) = NaN;
            continue;
        end

        PK=lyap((A+B*K).',Q+K.'*R*K);
        if isempty(PK) == 1 % ARE might not have solution since hinf constraint not satisfied
            Jlqr(i,j) = NaN;
        else
            Jlqr(i,j) = trace(Omega*PK);
        end
        % Restrict Jlqr values between 0 and 20
        Jlqr(i, j) = min(max(Jlqr(i, j), 0), 25);
    end
end
[k1,k2] = meshgrid(k1_range,k2_range);
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [7, 6]; %[-k,-k] achieve min
FontSize = 11;

surfc(k1', k2', Jlqr, 'EdgeColor', 'none');
hold on; %contour(k1', k2', Jlqr, 5); hold on;
plot3(0,0,JKstar,'o','Color','r','MarkerSize',6,'MarkerFaceColor','r')
zlim([-5,25])
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');

xlabel('$k_1$','Interpreter','latex','FontSize',FontSize);
ylabel('$k_2$','Interpreter','latex','FontSize',FontSize);
zlabel('$J_{\mathrm{LQR}}(K)$','Interpreter','latex','FontSize',FontSize);
% xticks(-5.5:0.2:0-5)
% yticks(6:0.2:6.5)
clim([0 25]); 
%colorbar; 
