% Fig 8(b) in "Benign Nonconvex Landscapes in Optimal and Robust Control,
%           Part II: Extended Convex Lifting"
%       By Yang Zheng, Chih-Fan Pai, Yujie Tang

clear all;close all;clc;
%%
A = -eye(2);
[nx,~] = size(A);
B = eye(2);
[~,nu] = size(B);
Bw = eye(nx);
W = Bw*Bw.';
Q = eye(nx);
R = eye(nu);
C = [sqrtm(Q);zeros(nu,nx)];
D = [zeros(nx,nu);sqrtm(R)];
%% can be commented for the plot
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
n1=100; n2=100; r=1;
k1_range = linspace(KinfLMI(1,1)-1.1,KinfLMI(1,1)+0.8,n1);
k2_range = linspace(KinfLMI(2,2)-1.1,KinfLMI(2,2)+0.8,n2);
Jinf = zeros(n1,n2);
k_line = linspace(-2,-0.2,n1);
Jinf_line = zeros(n1,1);
for i = 1:n1
    K = k_line(i)*eye(2);
    if max(real(eig(A+B*K))) >= 0
        Jinf_line(i,1) = NaN;
        continue;
    else
        G = ss(A + B*K,Bw,sqrtm(Q+K.'*R*K),0);
        Jinf_line(i,1) = norm(G,'inf');
    end
end
for i = 1:n1
    for j=1:n2
        K = [k1_range(i) KinfLMI(1,2);KinfLMI(2,1) k2_range(j)];
        % K = [k1_range(i) k2_range(j);KinfLMI(2,1) KinfLMI(2,2)];
        if max(real(eig(A+B*K))) >= 0
            Jinf(i,j) = NaN;
            continue;
        else
            G = ss(A + B*K,Bw,sqrtm(Q+K.'*R*K),0);
            Jinf(i,j) = norm(G,'inf');
        end
    end
end

[k1,k2] = meshgrid(k1_range,k2_range);
%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [7, 6]; %[-k,-k] achieve min
FontSize = 11;

surfc(k1', k2', Jinf,'EdgeColor', 'none','FaceAlpha',0.9); 
hold on;

plot3(k_line, k_line, Jinf_line,'-','Color','r','linewidth',1.2)
hold on;

plot3(k_line, k_line, repmat(0.69,n1),'--','Color','r','linewidth',1.2)
hold on;

plot3(KinfLMI(1,1),KinfLMI(2,2),Jinf_KinfLMI,'o','Color','r','MarkerSize',5,'MarkerFaceColor','r')
hold on;

plot3(KinfLMI(1,1),KinfLMI(2,2),0.69,'o','Color','r','MarkerSize',5,'MarkerFaceColor','r')
hold on;

zlim([0.69,0.8]);

xlabel('$k_1$','Interpreter','latex','FontSize',FontSize);
ylabel('$k_2$','Interpreter','latex','FontSize',FontSize);
zlabel('$J_{\infty}(K)$','Interpreter','latex','FontSize',FontSize);

set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex');
clim([0.67 0.83]); 
view(-80,25)
%print(gcf,'HinfSF_StationaryPoint_K_1by2_SimpleExample.eps','-depsc2','-r300');