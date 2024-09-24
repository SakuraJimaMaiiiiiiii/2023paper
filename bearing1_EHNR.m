%% Slip_Simulated_Signal_Analysis

clear
close all
clc

%% chosen signal
% original signal
M = csvread('99.csv',1,0);
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


%% Q,R inter
Q_range = 1:1:10;               
r_range = 3:1:5;
ka_range = 0.2:0.1:0.9;
EHNR = zeros();
enp_sig_sqrt =zeros();
et=zeros();
delt_p=1.5;
gamma = 0.85;
T1 = 1/107.9;           %外圈故障频率
for i = 1 :length(Q_range)
    Q= Q_range(i);
    for j = 1:length(r_range)
        r = r_range(j);
        beta=2/(Q+1);
        alpha=1-beta/r;
        Jmax=floor(log(beta*length(y)/8)/log(1/alpha));
        J = min(10,Jmax);
        for k=1:length(ka_range)
            ka = ka_range(k);
            [x,~]=TQWT_SR_GMC_penalty_fun(y',Q,r,J,gamma,ka,0);
            y_GMC = itqwt(x,Q,r,length(y));
            enp_sig_trans = abs(hilbert(y_GMC));
            [et_corr,lags] = xcorr(enp_sig_trans);         %自相关 第一个系数为值 第二个为时滞
            [peaks0,loc_pks0]=findpeaks(et_corr,'minpeakdistance',round(fs*T1*4/5));
            [peaks1,loc_pks1]=findpeaks(et_corr(length(y_GMC)+1:end),'minpeakdistance',round(fs*T1*4/5));
            EHNR(i,j,k) = max(peaks1)/(max(peaks0)-max(peaks1));     
        end
    end
end

%% choose Q,r,J
[i,j]=max(EHNR);
i=squeeze(i);
j=squeeze(j);
[i_,j_]=max(i);
[i__,j__]=max(i_);
max_k=j__;
max_j=j_(max_k);
max_i=j(max_j,max_k);
Q = Q_range(max_i);
r = r_range(max_j);
ka = ka_range(max_k);
J = min(10,computeJmax(length(y),Q,r));

%% TQWT sparse representation
[x,v]=TQWT_SR_GMC_penalty_fun(y',Q,r,J,gamma,ka,0);
y_GMC = itqwt(x,Q,r,length(y));
figure,plot(t,guzhang_signal(1:20000),'r--',t,y_GMC,'b')
y_GMC_blp = abs(fft(abs(hilbert(y_GMC))))/length(y_GMC)*2;
y_GMC_blp(1) = 0;
bl_fr = (0:length(y_GMC)-1)/length(y_GMC)*fs;
figure,plot(bl_fr(1:fix(length(bl_fr)/2)),y_GMC_blp(1:fix(length(bl_fr)/2)))