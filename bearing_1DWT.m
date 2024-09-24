%% Slip_Simulated_Signal_Analysis

clear
close all
clc

%% chosen signal
% original signal
M = csvread('60.csv',1,0);
sig0=M(:,1);
figure,plot(sig0,'y--')
fs=26500;  % sampling freq.

% 20000length signal
y = sig0;
y = y(1:20000);
figure,plot(y,'r--')
xlabel('time [s]'),ylabel('Amplitude')
axis tight
blp = abs(y);
figure,plot(blp)
t=(0:length(y)-1)/fs;
figure,plot(t,y,'b')
xlabel('time [s]'),ylabel('Amplitude')

% guzhang signal
t_end= 1.28;
tt = length(sig0)/t_end/fs; 
guzhang_signal = sig0(1:tt:length(sig0));
t=(0:length(y)-1)/fs;
figure,plot(t,guzhang_signal(1:length(t)),'b--');

%% DB9 离散小波变换
N=4;         %分解层数
n=20000;
[C,L]=wavedec(y,N,'db6');
A=zeros(N+1,n);       %储存重构信号，按行，由低频到高频，逼近到细节信号
figure();
% title('DB6 DWT')
for i=3
    if i==1
        A(i,:)=wrcoef('a',C,L,'db6',N);
        subplot(N+1,1 ,i),plot(t,A(i,:))
        ylim([min(A(i,:)-0.005) max(A(i,:))+0.005])
        set(gca,'xtick',[],'xticklabel',[])
       
    else
        A(i,:)=wrcoef('d',C,L,'db6',N+2-i);
        subplot(N+1,1 ,i),plot(t,A(i,:))
        ylim([min(A(i,:)-0.005) max(A(i,:))+0.005])
        set(gca,'xtick',[],'xticklabel',[])
%     set(gca,'ytick',[],'yticklabel',[])    
%     axis off   %去掉横纵坐标轴
    end
end
figure,plot(A(i,:),'color',0.5*[1 1 1])
y_DWT_blp = abs(fft(abs(hilbert(A(i,:)))))/length(A(i,:))*2;
y_DWT_blp(1) = 0;
bl_fr = (0:length(A(i,:))-1)/length(A(i,:))*fs;
figure,plot(bl_fr(1:fix(length(bl_fr)/2)),y_DWT_blp(1:fix(length(bl_fr)/2)),'color',0.5*[1 1 1])
xlabel('频率 [HZ]')
ylabel('幅值')
xlim([0,450])
ylim([-1.5,1.5])
%  
% 
% xlim([0,3000])
%  xlabel('时间[s]')
%  ylabel('幅值')

% set(gca,'ytick',[])  %隐去y轴坐标值
%  set(gca,'xtick',[])  %隐去x轴坐标值