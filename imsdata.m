clear;
close all;
clc

datapath = 'E:\files\tqwt_code_paper\故障轴承数据\2nd_test\2004.02.19.04.42.39';
data = load(datapath);
sig0=data(:,1);
figure,plot(sig0)
xlim([0.1*10^4,1*10^4])
xlabel('时间 [s]'),ylabel('幅值')
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
fs = 20000;
% 20000length signal
y = sig0;
y = y(1:20000);
figure,plot(y)
xlim([0.1*10^4,1*10^4])
xlabel('时间 [s]'),ylabel('幅值')
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})


% blp = abs(y);
% figure,plot(blp)
% t=(0:length(y)-1)/fs;
% figure,plot(t,y)
% xlabel('时间 [s]'),ylabel('幅值')

t_end= 1;
tt = length(sig0)/t_end/fs; 
guzhang_signal = sig0(1:tt:length(sig0));
t=(0:length(y)-1)/fs;
figure,plot(t,guzhang_signal(1:length(t)),'b--')
x=guzhang_signal;

%% TQWT denoising
x=x(1:20000);
xzero=zeros(size(x));
N=length(x);
Q = 3; r = 3; J = 15; % High Q-factor wavelet transform
Jtext=computeJmax(N,Q,r);
J=min(15,Jtext);
[wlets, now] = ComputeWavelets(N,Q,r,J,'radix2');
w = tqwt(x,Q,r,J);              % TQWT
wzero = tqwt(xzero,Q,r,J);      %follow sig

norm_w=w{1}/now(1);      %guiyihua
for i=2:length(now)
    norm_w=[norm_w(:);w{i}(:)/now(i)];
end

%% definition of the threshold
thres=mean(norm_w)+std(norm_w)*1.5;
% thres=thselect(norm_w,'minimaxi'); % 'heursure','sqtwolog','minimaxi','rigrsure'
%% wavelet coefficient 
alpha=0.9;
for i=1:J+1
    wzero{i}=compute_soft(w{i},thres,alpha);
end
%% inverse TQWT and obtain the denoised result
y = itqwt(wzero,Q,r,length(x));  
                                  %figure,plot(y)      reconstrust sig
                                  %figure,plot(sig0(1:length(y)),'b')


figure,plot(y)
xlim([0.1*10^4,1*10^4])
xlabel('时间 [s]'),ylabel('幅值')
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
ylabel('幅值')



f = linspace(0,fs/2,N/2);
A1 = abs(y)/(N/2);
figure,plot(f,A1(1:N/2))


blp=abs(fft(abs(hilbert(y))))/length(y)*2;
blp(1)=0;
pl=(0:length(y)-1)/length(y)*fs;
figure,plot(pl(1:round(length(y)/2)),blp(1:round(length(y)/2)))
xlabel('频率 [HZ]')
ylabel('幅值')
xlim([0,1200])


