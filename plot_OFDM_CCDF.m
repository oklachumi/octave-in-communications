clear,close,clc all
% CCDF of OFDM Signal

%function [mod_object] = mapper(b,N)
%% If N is given, it generates a block of N random 2^b-PSK/QAM modulated symbols.
%% Otherwise, it generates a block of 2^b-PSK/QAM modulated symbols for [0:2^b-1].
%
%M=2^b; % Modulation order or Alphabet (Symbol) size
%if b==1
%  Mod='BPSK';
%  A=1;
%  mod_object=pskmod(N,2);
%elseif b==2
%  Mod='QPSK';
%  A=1;
%  mod_object=pskmod(N,4,pi/4,'gray');
%else
%  Mod=[num2str(2^b) 'QAM'];
%  Es=1;
%  A=sqrt(3/2/(M-1)*Es);
%  mod_object=qammod(N,M);
%end
%endfunction

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

Ns = 2.^[6:10];
b=2;
M=2^b;
Nblk = 1e3;

zdBs = [4:0.1:10]; % x-axis
N_zdBs = length(zdBs);

CCDF_formula=inline('1-((1-exp(-z.^2)).^N)','N','z');
for n = 1:length(Ns) % n = 1:5    
  N=Ns(n); % N = 64, 128,..., 1024
  x=zeros(Nblk,N); % 1000*N
  sqN=sqrt(N);
  for k = 1:Nblk
    XX = randint(1,N,M);
    X = pskmod(XX,M,pi/4,'gray');
    x(k,:) = ifft(X,N)*sqN;
    CFx(k) = PAPR(x(k,:)); % Cumulative Distribution Function
  end
  %s2 = mean(mean(abs(x)))^2/(pi/2);
  CCDF_theoretical=CCDF_formula(N,10.^(zdBs/20)); % Complementary CDF
  for i = 1:N_zdBs
    CCDF_simulated(i) = sum(CFx>zdBs(i))/Nblk;
  end
  semilogy(zdBs,CCDF_theoretical,'k-');
  hold on;
  grid on;
  semilogy(zdBs(1:3:end),CCDF_simulated(1:3:end),'k:*');
end
axis([zdBs([1 end]) 1e-2 1]);
title('OFDM system with N-point FFT');
xlabel('PAPR0 [dB]');
ylabel('CCDF=Probability(PAPR>PAPR0)');
legend('Theoretical','Simulated');
