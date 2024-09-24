%% hard/soft threshold tqwt denoise

clear
close all
clc

%% Two fault modes with common resonant frequency
fr1=3;  % convolution freq.
a1=1;    % amplitude for sig1;
fn1=625; % resonant freq.
zeta=400; % decaying rate
T1=1/16;  % fault frequency for sig1;
fs=10000;  % sampling freq.

% Generating the signal
t1=0:1/fs:T1;
sig_1=(a1*exp(-zeta*t1).*cos(2*pi*fn1*t1));  % single impulse
figure,plot(sig_1)
% xlabel('时间[s]')
% ylabel('幅值')
% xlim([0,400])
% title('冲击单元')
t=0:1/fs:2;
sig0=zeros(size(t));
for i=1:fix(t(end)/T1)
    if i==1 || i==fix(t(end)/T1)
        inde=max(find(t<=T1*(i-1)));
    else
        inde=max(find(t<=T1*(i-1)+(rand(1)-0.5)*2*1/100*T1));
%         inde=max(find(t<=T1*(i-1)+dd(mod(round(rand(1)*10),3)+1)*13/100*T1));
%         inde=max(find(t<=T1*(i-1)+dd(randperm(3,1))*0.75/100*T1));
    end
    sig0(inde:inde+length(sig_1)-1)=sig_1;
end
% sig1=(cos(2*pi*fr1*t)+1).*sig0;  % fault signal 1
sig1=sig0;
sig=awgn(sig1,-8,'measured');
% figure,plot(t,sig,'b')
figure,plot(t,sig,'color',0.5*[1 1 1])
hold on
plot(t,sig1,'b-')
legend('含噪信号','故障信号')
axis tight
xlim([0 2])
xlabel('time [s]'),ylabel('幅值')
% frequency analysis
pp=abs(fft(sig))/length(sig)*2;
pp(1)=0;
ff=(0:length(pp)-1)/length(pp)*fs;
n1=floor(length(pp)/2);
figure,plot(ff(1:n1),pp(1:n1))
%  xlim([0 1000])
xlabel('频率 [Hz]'),ylabel('幅值')
x=sig;

t=t(1:length(t)-1);
x=x(1:length(x)-1);
figure,plot(t,x,'b')
xlabel('time [s]'),ylabel('幅值')


%% denoise

soft=wthresh(x,'s',0.4);
figure,plot(soft)
xlabel('时间 [s]'),ylabel('幅值')'
xlim([0.4*10^4,1*10^4])
xticklabels({'0.4','0.5','0.6','0.7','0.8','0.9','1'})
title('软阈值去噪')
% hard=wthresh(x,'h',0.3);
% figure,plot(hard)

lev=3;
n=5;
wf='db1';

hard=wden(x,'heursure','h','one',lev,wf);
figure,plot(hard)
xlabel('时间 [s]'),ylabel('幅值')
xlim([0.4*10^4,1*10^4])
xticklabels({'0.4','0.5','0.6','0.7','0.8','0.9','1'})
title('硬阈值去噪')



%% Frequency Analysis

blp=abs(fft(abs(hilbert(soft))))/length(soft)*2;
figure,plot(blp)

blp(1)=0;
pl=(0:length(soft)-1)/length(soft)*fs;
figure,plot(pl(1:round(length(soft)/2)),blp(1:round(length(soft)/2)))
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,250])
title('软阈值去噪信号包络谱')

blp=abs(fft(abs(hilbert(hard))))/length(hard)*2;
figure,plot(blp)

blp(1)=0;
pl=(0:length(hard)-1)/length(hard)*fs;
figure,plot(pl(1:round(length(hard)/2)),blp(1:round(length(hard)/2)))
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,250])
title('硬阈值去噪信号包络谱')



