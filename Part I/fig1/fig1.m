
% Figure 1: measure zero for degenerate policies

% Housekeeping
clear
close all

% Make figure
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [10, 6.8];



FontSize = 11;
%----------------------------------------------

% LQG data

load data_LQG.mat;

pos_x = [0.05 0.38 0.71 0.05 0.38 0.71];
pos_y = [0.58 0.58 0.58 0.08 0.08 0.08];

Ck = [-1, -0.8, -0.4, 0, 0.2, ...
    0.4, 0.6, 0.75, 0.9];
cc = 1;
for kk = 4:9 
    axes1 = axes('Parent',fig,'Position',[pos_x(cc) pos_y(cc) 0.26 0.3]);
    hold(axes1,'on');
    surf(Ak_grid(:, :, kk), Bk_grid(:, :, kk), log(abs(P12vals(:, :, kk))), ...
        'Parent',axes1, 'EdgeColor', 'none');
    title(strcat('$C_{\mathsf{K}}= \;$',string(Ck(kk))),'FontSize',FontSize,...
        'Interpreter','latex');
    grid(axes1,'on');
    hold(axes1,'off');
    set(axes1,'CLim',[-7 2],'FontSize',FontSize,'TickLabelInterpreter','latex',...
        'XTick',[-2 -1 0 1 2]);
    xlim([-1,1]);
    ylim([-2,1]);
    
    cc = cc + 1;
end

print(gcf,'heatmap_LQG_yz.eps','-depsc2','-r300');

%% Hinf data

% Make figure
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [10, 6.8];


load data_hinf.mat;

pos_x = [0.05 0.38 0.71 0.05 0.38 0.71];
pos_y = [0.58 0.58 0.58 0.08 0.08 0.08];

Dk = [-1, -0.8, -0.4, 0, 0.2, ...
    0.4, 0.6, 0.75, 0.9];
tmp = [11 15 23 31 35 39 43 46 49];
cc = 1;
for kk = tmp(4:9) 
    axes1 = axes('Parent',fig,'Position',[pos_x(cc) pos_y(cc) 0.26 0.3]);
    hold(axes1,'on');
    surf(Ak_grid(:, :, kk), Bk_grid(:, :, kk), log(abs(P12vals(:, :, kk))), ...
        'Parent',axes1, 'EdgeColor', 'none');
    title(strcat('$D_{\mathsf{K}}= \;$',string(Dk(cc+3))),'FontSize',FontSize,'Interpreter','latex');
    grid(axes1,'on');
    hold(axes1,'off');
    set(axes1,'CLim',[-7 2],'FontSize',FontSize,'TickLabelInterpreter','latex',...
        'XTick',[-2 -1 0 1]);
    xlim([-2,1]);
    ylim([-4,4]);
    
    cc = cc + 1;
end

print(gcf,'heatmap_Hinf_yz.eps','-depsc2','-r300');
