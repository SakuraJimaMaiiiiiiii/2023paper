% % 定义信号信息
% fs=2^6;    %采样频率
% dt=1/fs;    %时间精度
% timestart=-8;
% timeend=8;
% t=(0:(timeend-timestart)/dt-1)*dt+timestart;
% L=length(t);
% 
% z=4*sin(2*pi*linspace(6,12,L).*t);

%% chosen signal
% original signal
M = csvread('106.csv',1,0);
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
z=y;
%matlab自带的小波变换
%新版本
figure(1)
[wt,f,coi] = cwt(z,'amor',fs);
pcolor(t,f,abs(wt));
xlim([0.1,0.15])

shading interp
