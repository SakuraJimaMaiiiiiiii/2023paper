%% Slip_Simulated_Signal_Analysis

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
% figure,plot(sig_1)
% xlabel('时间[s]')
% ylabel('幅值')
% xlim([0,300])

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

%% TQWT test
Q = 3; r = 3; J = 15;     % High Q-factor wavelet transform
N=200;
w = tqwt(x,Q,r,J);              % TQWT
y = itqwt(w,Q,r,length(x));             % Inverse TQWT
figure,plot(x,'b');
hold on
plot(y,'r--');
plot(y-x,'g')

J1 = 1; J2 = J;
figure, clf
PlotWavelets(N,Q,r,J1,J2);
figure,
PlotFreqResps(Q, r, J)


%% TQWT denoising
load x;
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
for i=2:length(now);
    norm_w=[norm_w(:);w{i}(:)/now(i)];
end

%% definition of the threshold
thres=mean(norm_w)+std(norm_w)*2.5;
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
xlabel('时间[s]'),ylabel('幅值')
xlim([0.4*10^4,1*10^4])
xticklabels({'0.4','0.5','0.6','0.7','0.8','0.9','1'})

figure,plot(sig1(1:length(y)),'b')
hold on 
plot(y,'r--') 
xlabel('time [s]'),ylabel('幅值')
legend('仿真故障特征信号','重构故障特征信号')

blp=abs(fft(abs(hilbert(y))))/length(y)*2;
figure,plot(blp)

blp(1)=0;
pl=(0:length(y)-1)/length(y)*fs;
figure,plot(pl(1:round(length(y)/2)),blp(1:round(length(y)/2)))
xlabel('频率 [Hz]'),ylabel('幅值')
xlim([0,250])