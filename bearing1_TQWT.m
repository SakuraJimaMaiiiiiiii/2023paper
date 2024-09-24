%% Slip_Simulated_Signal_Analysis

clear
close all
clc


%% chosen signal
% original signal
M = csvread('60.csv',1,0);
sig0=M(:,1);
figure,plot(sig0,'color',0.5*[1 1 1])
xlim([0,10000])
ylim([-9,9])
fs=26500;  % sampling freq.

% 20000length signal
x = sig0;
x = x(1:20000);
figure,plot(x)
xlabel('时间 [s]'),ylabel('幅值')
xlim([1*10^4,1.5*10^4])
xticklabels({'1','1.05','1.1','1.15','1.2','1.25','1.3','1.35','1.4','1.45','1.5'})



%% TQWT denoising
xzero=zeros(size(x));
soft = @(x, T) max(1 - T./abs(x), 0) .* x;
N=length(x);
Q = 3; r = 3; J = 15; % High Q-factor wavelet transform
%  J=min(15,floor(jisuandezuida));
[wlets, now] = ComputeWavelets(N,Q,r,J,'radix2');
w = tqwt(x,Q,r,J);              % TQWT
wzero = tqwt(xzero,Q,r,J);      %follow sig

norm_w=w{1}/now(1);      %guiyihua
for i=2:length(now);
    norm_w=[norm_w(:);w{i}(:)/now(i)];
end

thres=mean(norm_w)+std(norm_w)*2.5;

for i=1:length(now)
    wzero{i}=soft(w{i},thres*now(i));
end 
y = itqwt(wzero,Q,r,length(x));  
                                  %figure,plot(y)      reconstrust sig
                                  %figure,plot(sig0(1:length(y)),'b')
figure,plot(y,'color',0.5*[1 1 1])
xlabel('时间[s]')
ylabel('幅值')
xlim([0,3000])
set(gca,'ytick',[])  %隐去y轴坐标值
set(gca,'xtick',[])  %隐去x轴坐标值

%%
blp=abs(fft(abs(hilbert(y))))/length(y)*2;
                                  %figure,plot(blp)
blp(1)=0;
pl=(0:length(y)-1)/length(y)*fs;
figure,plot(pl(1:round(length(y)/2)),blp(1:round(length(y)/2)),'color',0.5*[1 1 1])
xlim([0,450])
xlabel('频率[Hz]')
ylabel('幅值')
set(gca,'ytick',[])  %隐去y轴坐标值
%  set(gca,'xtick',[])  %隐去x轴坐标值
