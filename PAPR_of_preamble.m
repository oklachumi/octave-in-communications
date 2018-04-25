% PAPR_of_preamble.m
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

N = 1024;
L = 4;
Npreamble = 114;
n = 0:Npreamble-1;
for i = 1:Npreamble
  X = load(['.\\Wibro-Preamble\\Preamble_sym' num2str(i-1) '.dat']);
  X = X(:,1);
  X = sign(X);
  X = fftshift(X);
  x = IFFT_oversampling(X,N);
  PAPRdB(i) = PAPR(x);
  x_os = IFFT_oversampling(X,N,L);
  PAPRdB_os(i) = PAPR(x_os);
end
plot(n,PAPRdB,'-ro', n,PAPRdB_os,':*');
title('PAPR of IEEE 802.16e preamble without and with oversampling');
ylabel('|IFFT(X_i(k))|');
xlabel('Preamble index [0-113]');
legend('L = 1','L = 4', 'location','southeast');
