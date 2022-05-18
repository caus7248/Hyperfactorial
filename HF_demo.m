clear;
clc;
close all;
x   = -4:0.025:4;
y   = -4:0.025:4;
[X,Y] = meshgrid(x,y);
Z   = X+1i*Y;
for k =1:50
	m = -5+10*k/50;
    tic;
    [H,L] = HF(m,Z);
%      save(['HL_data' num2str(k+1000) '.mat'],"m","X","Y","L","H",'-mat');
set(gcf,'color','white','Position',[50   50   1400   700])
subplot(1,2,1);
dcolor(X,Y,L);  %phase plot of L = log(H)
xlabel('\Re{z}');
ylabel('\Im{z}');
title(['ln H(' num2str(m) ',z)']);
colorbar;

subplot(1,2,2);
dcolor(X,Y,H);  %phase plot of H
xlabel('\Re{z}');
ylabel('\Im{z}');
title(['H(' num2str(m) ',z)']);
colorbar;

drawnow;
% saveas(gca, ['HL_plot' num2str(k+1000) '.png']);
    toc;
end