clear all; close all; clc;
%%
load("dim1.mat"); dim1 = values;
load("dim2.mat"); dim2 = values;
load("dim3.mat"); dim3 = values;
load("AC8.mat"); AC8 = values;
load("HE1.mat"); HE1 = values;
load("REA2.mat"); REA2 = values;

BestValues = [dim1(end),dim2(end),5.0829,1.6165,0.0736,REA2(end)];

%%
lwidth = 1.5;
FontSize = 12;
fig = figure();
fig.WindowStyle = 'normal';
fig.Units = 'centimeters';
fig.Position([3,4]) = [11, 6.5];

semilogy(1:length(dim1),(dim1-BestValues(1))./(1+BestValues(1)),'-','LineWidth',lwidth);
hold on;
semilogy(1:length(dim2),(dim2-BestValues(2))./(1+BestValues(2)),'-','LineWidth',lwidth);
semilogy(1:length(dim3),(dim3-BestValues(3))./(1+BestValues(3)),'-','LineWidth',lwidth);
semilogy(1:length(AC8),(AC8-BestValues(4))./(1+BestValues(4)),'-','LineWidth',lwidth);
semilogy(1:400,(HE1(1:400)-BestValues(5))./(1+BestValues(5)),'-','LineWidth',lwidth);
semilogy(1:length(REA2),(REA2-BestValues(6))./(1+BestValues(6)),'-','LineWidth',lwidth);

h = legend('1-dim','2-dim','3-dim','AC8','HE1','REA2','Interpreter','latex','FontSize',FontSize);
set(h,'box','off')

xlabel('Iterations','FontSize',FontSize,'Interpreter','latex')
ylabel('Suboptimality','FontSize',FontSize,'Interpreter','latex');
set(gca,'FontSize',FontSize,'TickLabelInterpreter','latex',...
    'XTick',[0:100:600],'YTick',[10^(-5),10^(-4),10^(-2),10^0]);
grid off; xlim([0,600]); ylim([0.7*10^(-5),10])
print(gcf,'Hinf_hifoo.eps','-depsc2','-r300');
set(gca,'FontSize',FontSize) 