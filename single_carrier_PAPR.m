clear,close,clc all

function [s,time] = modulation(x,Ts,Nos,Fc)
% modulation(X,1,32,1)
% Ts : Sampling period
% Nos: Oversampling factor
% Fc : Carrier frequency

Nx = length(x); % 4
offset = 0; 
if nargin < 5
  scale = 1;
  T = Ts/Nos; % Scale and Oversampling period for Baseband
else
  scale = sqrt(2);
  T=1/Fc/2/Nos; % Scale and Oversampling period for Passband 
end
t_Ts = [0:T:Ts-T];
time = [0:T:Nx*Ts-T]; % One sampling interval and whole interval
tmp = 2*pi*Fc*t_Ts+offset;
len_Ts = length(t_Ts); % 8
cos_wct = cos(tmp)*scale;
sin_wct = sin(tmp)*scale;
for n = 1:Nx
  s((n-1)*len_Ts+1:n*len_Ts) = real(x(n))*cos_wct-imag(x(n))*sin_wct;
end
endfunction

function [PAPR_dB, AvgP_dB, PeakP_dB] = PAPR(x)
Nx = length(x);
xI = real(x);
xQ = imag(x);
Power = xI.*xI + xQ.*xQ;
PeakP = max(Power);
PeakP_dB = 10*log10(PeakP);
AvgP = sum(Power)/Nx;
AvgP_dB = 10*log10(AvgP);
PAPR_dB = 10*log10(PeakP/AvgP);
endfunction

%single_carrier_PAPR.m
figure(1); clf
Ts=1;   % Sampling period
L=8;    % Oversampling factor
Fc=1;   % Carrier frequency
b=2;
M=2^b;  % Modulation order or Alphabet size

XX = randint(1,4,M); % 1 3 2 2
X = pskmod(XX,M,pi/4,'gray'); % M-PSK/QAM symbol for [0:M-1]
L_=L*4; % Oversampling factor to make it look like continuous-time
[xt_pass_,time_] = modulation(X,Ts,L_,Fc); % Continuous-time
[xt_pass,time] = modulation(X,Ts,L,Fc); % L times oversampling

for i_s=1:M
  xt_base(L*(i_s-1)+1:L*i_s) = X(i_s)*ones(1,L);
end

PAPR_dB_base = PAPR(xt_base);
subplot(311); stem(time,real(xt_base),'b'); hold on; ylabel('S_{I}');
subplot(312); stem(time,imag(xt_base),'b'); hold on; ylabel('S_{Q}');
subplot(313); stem(time,abs(xt_base).^2,'b'); hold on;
title(['PAPR = ' num2str(round(PAPR_dB_base*100)/100) 'dB']);
xlabel('samples'); ylabel('|S_{I}(n)|^{2}+|S_{Q}(n)|^{2}');

figure(2); clf
PAPR_dB_pass = PAPR(xt_pass);
subplot(211); stem(time,xt_pass,'r'); hold on;
plot(time_,xt_pass_,'k'); ylabel('S(n)');

subplot(212); stem(time,xt_pass.*xt_pass,'r'); hold on;
plot(time_,xt_pass_.*xt_pass_,'k');
title(['PAPR = ' num2str(round(PAPR_dB_pass*100)/100) 'dB']);
xlabel('samples'); ylabel('|S(n)|^{2}');
% PAPRs of baseband/passband signals
PAPRs_of_baseband_passband_signals=[PAPR_dB_base; PAPR_dB_pass]
