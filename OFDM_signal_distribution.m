% OFDM_signal_distribution.m
clear,close,clc all

N=8;
b=2;
M=2^b;
Nos=16; % oversampling factor = 16
NNos=N*Nos;
T=1/NNos;
time=[0:T:1-T];
XX=randint(1,N,M);
X=pskmod(XX,M,pi/4,'gray');
X(1) = 0+j*0; % A block of 16 QPSK symbols with no DC-subcarrier 
for i = 1:N
  if i <= N/2
    x = ifft([zeros(1,i-1) X(i) zeros(1,NNos-i+1)],NNos);
  else
    x = ifft([zeros(1,NNos-N+i-1) X(i) zeros(1,N-i)],NNos);
  end
  xI(i,:) = real(x);
  xQ(i,:) = imag(x);
end
sum_xI = sum(xI);
sum_xQ = sum(xQ);

figure(1); clf
subplot(311)
plot(time,xI,'k:','linewidth',1); hold on;
plot(time,sum_xI,'b','linewidth',2);
title(['N=' num2str(N)]);
ylabel('x_{I}(t)');
axis([0 1 min(sum_xI) max(sum_xI)]);

subplot(312)
plot(time,xQ,'k:','linewidth',1); hold on;
plot(time,sum_xQ,'b','linewidth',2);
ylabel('x_{Q}(t)');
axis([0 1 min(sum_xQ) max(sum_xQ)]);

subplot(313)
plot(time,abs(sum_xI+j*sum_xQ),'b','linewidth',2); hold on;
ylabel('|x(t)|'); xlabel('t');

clear('xI');
clear('xQ');

N=2^4; % N = 16
NNos=N*Nos; % NNos = 16*16 = 256
T=1/NNos;
time=[0:T:1-T];
Nhist=1e3;
for k = 1:Nhist
  XX=randint(1,N,M);
  X=pskmod(XX,M,pi/4,'gray'); 
  X(1)=0+j*0; % A block of 16 QPSK symbols with no DC-subcarrier 
  for i = 1:N
    if i<= N/2
      x = ifft([zeros(1,i-1) X(i) zeros(1,NNos-i+1)],NNos);
    else
      x = ifft([zeros(1,NNos-N/2+i-N/2-1) X(i) zeros(1,N-i)],NNos);
    end
    xI(i,:) = real(x);
    xQ(i,:) = imag(x);
  end
  HistI(NNos*(k-1)+1:NNos*k) = sum(xI);
  HistQ(NNos*(k-1)+1:NNos*k) = sum(xQ);
end

N_bin = 30;
figure(2); clf

subplot(311)
% [counts(個數),centers(在這範圍內的均值)] = hist(___)
% bar(centers,counts)
[xI_dist,bins] = hist(HistI,N_bin);
bar(bins,xI_dist/sum(xI_dist),'k');
title(['N=' num2str(N)]);
ylabel('pdf of x_{I}(t)');

subplot(312)
[xQ_dist,bins] = hist(HistQ,N_bin);
bar(bins,xQ_dist/sum(xQ_dist),'k');
ylabel('pdf of x_{Q}(t)');

subplot(313)
[xabs_dist,bins] = hist(abs(HistI+j*HistQ),N_bin);
bar(bins,xabs_dist/sum(xabs_dist),'k');
ylabel('pdf of |x(t)|');
xlabel('x_{0}');
