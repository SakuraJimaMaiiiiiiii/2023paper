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

%% fast Kurtogram
nlevel=4;
c = Fast_Kurtogram_1(y,nlevel,fs);
xlim([0,400])
xlabel('时间（s）')
xlim([0,0.4])
ylim([-2,2])
ylabel('幅值')
set(gca,'ytick',[],'yticklabel',[]) 
% set(gca,'ytick',[])  %隐去y轴坐标值
%  set(gca,'xtick',[])  %隐去x轴坐标值