% PAPR_of_Chu.m
clear,close,clc all

function [xt, time] = IFFT_oversampling(X,N,L)
if nargin < 3
  L = 1;
end
NL = N*L;
T = 1/NL;
time = [0:T:1-T];
X = X(:).';
xt = L*ifft([X(1:N/2) zeros(1,NL-N) X(N/2+1:end)], NL);
endfunction

function [PAPR_dB, AvgP_dB, PeakP_dB] = PAPR(x)
% PAPR_dB  : PAPR[dB]
% AvgP_dB  : Average power[dB]
% PeakP_dB : Maximum power[dB]

Nx=length(x);
xI=real(x);
xQ=imag(x);
Power = xI.*xI + xQ.*xQ;
PeakP = max(Power);
PeakP_dB = 10*log10(PeakP);
AvgP = sum(Power)/Nx;
AvgP_dB = 10*log10(AvgP);
PAPR_dB = 10*log10(PeakP/AvgP);
endfunction

N = 16; % 16 point IFFT
L = 4;
i = [0:N-1]; 
k = 3; % gcd(k,N) = 1
X = exp(j*k*pi/N*(i.*i));
[x,time] = IFFT_oversampling(X,N);
PAPRdB = PAPR(x);
[x_os,time_os] = IFFT_oversampling(X,N,L);
PAPRdB_os = PAPR(x_os);

subplot(121)
plot(x,'ro');
hold on;
plot(x_os,'k*');
axis([-0.4 0.4 -0.4 0.4]);
%axis('equal');

plot(0.25*exp(j*pi/180*[0:359])) % circle with radius 0.25 

subplot(122)
plot(time,abs(x),'ro', time_os,abs(x_os),'k:*');
title('IFFT(X_i(k)), k=3, N=16, L=1, 4');
ylabel('|IFFT(X_i(k))|');
xlabel('Time (normalized by symbol duration)');
legend('L = 1','L = 4');
PAPRdB_without_and_with_oversampling=[PAPRdB  PAPRdB_os]
