clear all; close all; clc;
%%
load("dim1.mat"); dim1 = values;
load("dim2.mat"); dim2 = values;
load("dim3.mat"); dim3 = values;
load("dim1_gd.mat");
load("dim2_gd.mat");
load("dim3_gd.mat");
load("optimal_value.mat")

lwidth = 1.5;
FontSize = 12;

%%
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [13, 6];

axes1 = axes('Parent',fig,'Position',[0.15,0.3,0.32,0.6]);
hold(axes1,'on');
semilogy(1:length(dim1),(dim1-Jopt1)/Jopt1,'-','LineWidth',lwidth); hold on;
semilogy(1:length(dim2),(dim2-Jopt2)/Jopt2,'-','LineWidth',lwidth);
semilogy(1:length(dim3),(dim3-Jopt3)/Jopt3,'-','LineWidth',lwidth);
ylabel('Suboptimality','FontSize',FontSize,'Interpreter','latex');
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex',...
    'XTick',[0:40:160],'YTick',[10^(-6),10^(-4),10^(-2),10^0]);
xlabel('Iterations','FontSize',FontSize-1,'Interpreter','latex')
xlim([0,140]);ylim([10^(-6),10]);
grid off;set(gca,'YScale','log');
title('\texttt{HIFOO}','FontSize',FontSize-1,'Interpreter','latex')
box on

axes2 = axes('Parent',fig,'Position',[0.64,0.3,0.32,0.6]);
hold(axes2,'on');
h1 = semilogy(1:length(dim1_gd),dim1_gd,'-','LineWidth',lwidth); hold on
h2 = semilogy(1:500,dim2_gd(1:500),'-','LineWidth',lwidth);
h3 = semilogy(1:500,dim3_gd(1:500),'-','LineWidth',lwidth);
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex',...
    'XTick',[0:200:500],'YTick',[10^(-6),10^(-4),10^(-2),10^0]);
xlabel('Iterations','FontSize',FontSize-1,'Interpreter','latex')
xlim([0,400]);ylim([10^(-6),10]);
title('Gradient Descent','FontSize',FontSize-1,'Interpreter','latex')
grid off; set(gca,'YScale','log');
box on

h = legend([h1,h2,h3],'Instance 1','Instance 2','Instance 3','Interpreter','latex','FontSize',FontSize);
set(h,'box','off')
set(h, 'Position', [0.2, 0.01, .7, .1],'Orientation','horizontal');

print(gcf,'LQG_hifoo_gd.eps','-depsc2','-r300');
