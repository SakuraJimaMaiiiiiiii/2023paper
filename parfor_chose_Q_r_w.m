%% import signal
clear
clc
close all
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
sig=awgn(sig1,-5,'measured');
% figure,plot(t,sig,'b')
figure,plot(t,sig,'color',0.5*[1 1 1])
hold on
plot(t,sig1,'b-')
axis tight
xlim([0 2])
xlabel('time [s]'),ylabel('Amplitude')
% frequency analysis
pp=abs(fft(sig))/length(sig)*2;
pp(1)=0;
ff=(0:length(pp)-1)/length(pp)*fs;
n1=floor(length(pp)/2);
%figure,plot(ff(1:n1),pp(1:n1),'b')
% xlim([400 1000])
xlabel('Frequency [Hz]'),ylabel('Magnitude')
y=sig;
y=y(1:20000);
soft = @(x, T) max(1 - T./abs(x), 0) .* x;

%% signal parameters
fs=10000;  % sampling freq.
t=(0:length(y)-1)/fs;
figure,plot(t,y,'b')

%% TQWT parameters
% High Q-factor wavelet transform
gamma=0.85;
%%% w = tqwt(x,Q,r,J);              % TQWT
%%%y = itqwt(w,Q,r,length(x));             % Inverse TQWT




%% Q,R inter
Q_range = 1:1:10;               
r_range = 3:0.5:10;
ka_range = 0.3:0.1:1;
N=200;
snr =zeros();
num_r = length(r_range);
num_ka = length(ka_range);
parfor i = 1 :length(Q_range)
    Q= Q_range(i);
    for j = 1:num_r
        r = r_range(j);
        J = min(10,computeJmax(N,Q,r));
        for k=1:num_ka
            ka = ka_range(k);
            [x,~]=TQWT_SR_GMC_penalty_fun(y,Q,r,J,gamma,ka,0);
            y_GMC = itqwt(x,Q,r,length(y));
            snr(i,j,k)=SNR_singlech(y_GMC,sig1(1:length(y)));
        end
    end
end

%% choose Q,r,J
[i,j]=max(snr);
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
J = min(10,computeJmax(N,Q,r));

%% TQWT sparse representation
[x,v]=TQWT_SR_GMC_penalty_fun(y,Q,r,J,gamma,ka,0);
y_GMC = itqwt(x,Q,r,length(y));
figure,plot(t,sig1(1:length(y)),'r--',t,y_GMC,'b')

%% SIGNAL-NOISE RATIO
snr=SNR_singlech(y_GMC,sig1(1:length(y)));



%% kurtosis
ku=kurtosis(y_GMC);
