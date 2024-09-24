N =200;
Q = 1;
r = 3;
J = 20;
PlotWavelets(N,Q,r,5,5,'radix2');
xlabel('采样时间')
ylabel('分解层数')
title('Q=1 ,r=3 ')
set(gca,'xtick',[],'xticklabel',[]) 
set(gca,'ytick',[],'yticklabel',[]) 
ylabel('')

PlotFreqResps(Q, r, J)