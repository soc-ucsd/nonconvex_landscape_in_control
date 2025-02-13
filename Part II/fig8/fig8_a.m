% Fig 8(a) in "Benign Nonconvex Landscapes in Optimal and Robust Control,
%           Part II: Extended Convex Lifting"
%       By Yang Zheng, Chih-Fan Pai, Yujie Tang

clear all;close all;clc;
%%
A = -1;
[nx,~] = size(A);
B = 1;
[~,nu] = size(B);
Bw = 1;
W = Bw*Bw.';
Q = 0.1; R = 1;
C = [sqrtm(Q);zeros(nu,nx)];
D = [zeros(nx,nu);sqrtm(R)];
%%
% solve optimal state-feedback Hinf (obtain infimum of hinf norm) using LMI method
X = sdpvar(nx,nx);
Y = sdpvar(nu,nx,'full');
gamma = sdpvar(1,1);
obj = gamma;
eps = 1e-5;
constraints = [X >= eps*eye(nx,nx),...
    [A*X+X*A.'+B*Y+Y.'*B.', Bw, X*sqrtm(Q), Y.'*sqrtm(R);...
    Bw.', -gamma*eye(nx), zeros(nx,nx), zeros(nx,nu);...
    sqrtm(Q)*X, zeros(nx,nx), -gamma*eye(nx), zeros(nx,nu);...
    sqrtm(R)*Y, zeros(nu,nx), zeros(nu,nx), -gamma*eye(nu)] <= 0];
ops = sdpsettings('solver','mosek','verbose',1,'debug',1);
optimize(constraints,obj,ops);
JinfLMI= value(obj);
opt_X = value(X);
opt_Y = value(Y);
KinfLMI = opt_Y/opt_X;

% compute Hinf norm of the LMI policy
G = ss(A+B*KinfLMI,eye(nx),sqrtm(Q+KinfLMI.'*R*KinfLMI),0);
Jinf_KinfLMI = norm(G,'inf');
%% plot
n=50;
k_range = linspace(-4, 0.5, n);
Jinf = zeros(n,1);
s = tf('s');
for i = 1:n
    K = k_range(i);
    if max(real(eig(A+B*K)))>=0
        Jinf(i,1)=NaN;
        continue;
    else
        G = ss(A + B*K,Bw,sqrtm(Q+K.'*R*K),0);
        Jinf(i,1) = norm(G,'inf');
    end
end
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [6.5, 6];
FontSize = 11;

plot(k_range, Jinf,'m','linewidth',1.2); hold on;grid on;
plot(KinfLMI,JinfLMI,'o','Color','r','MarkerSize',4,'MarkerFaceColor','r')
xlabel('$k$','Interpreter','latex','FontSize',FontSize); 
ylabel('$J_{\infty}(k)$','Interpreter','latex','FontSize',FontSize);

set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');
%print(gcf,'HinfSF_StationaryPoint_K_1by1_SimpleExample.eps','-depsc2','-r300');