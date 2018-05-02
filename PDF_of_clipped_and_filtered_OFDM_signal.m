% PDF_of_clipped_and_filtered_OFDM_signal.m
% QPSK/OFDM system for analyzing the performance of clipping and filtering technique
clear,close,clc all

function [x_clipped,sigma] = clipping(x,CL,sigma)
% CL   : Clipping Level
% sigma: sqrt(variance of x)
if nargin < 3
  x_mean = mean(x);
  x_dev = x-x_mean;
  sigma = sqrt(x_dev*x_dev'/length(x));
end
CL = CL*sigma;
x_clipped = x;  
ind = find(abs(x)>CL); % Indices to clip
x_clipped(ind) = x(ind)./abs(x(ind))*CL;
endfunction

function [xt,time] = IFFT_oversampling(X,N,L)
if nargin < 3
  L = 1;
end
NL = N*L;
T = 1/NL;
time = [0:T:1-T];
X = X(:).';
xt = L*ifft([X(1:N/2) zeros(1,NL-N) X(N/2+1:end)], NL);
endfunction

function y = add_CP(x,Ncp)
% Add cyclic prefix
y = [x(:,end-Ncp+1:end) x];
endfunction

CR = 1.2; % Clipping Ratio
b = 2; % QPSK 2^b bits
N = 128; % FFT size
Ncp = 32; % CP size
fs = 1e6; % sampling freq
L = 8; % over sampling factor
Tsym = 1/(fs/N); % sampling freq
Ts = 1/(fs*L); % sampling period 
fc = 2e6;
wc = 2*pi*fc; % Carrier frequency
t = [0:Ts:2*Tsym-Ts]/Tsym; % time vector
t0 = t((N/2-Ncp)*L); % t(256)
f = [0:fs/(N*2):L*fs-fs/(N*2)]-L*fs/2;
Fs = 8;
Norder = 104;
dens = 20; % Sampling frequency, Order, and Density factor of filter
FF = [0 1.4 1.5 2.5 2.6 Fs/2]; % Stopband/Passband/Stopband frequency edge vector
WW = [6 1 6]; % Stopband/Passband/Stopband weight vector

h = remez(Norder,FF/(Fs/2),[0 0 1 1 0 0],WW,"bandpass",dens); % BPF coefficients

XX = randint(1,N,4);
X = pskmod(XX,4,pi/4,'gray');
% X = mapper(b,N);
X(1) = 0; % QPSK modulation

x = IFFT_oversampling(X,N,L); % IFFT and oversampling
x_b = add_CP(x,Ncp*L); % Add CP
x_b_os = [zeros(1,(N/2-Ncp)*L), x_b, zeros(1,N*L/2)]; % Oversampling 256(0) 256(CP) 1024(data) 512(0)
x_p = sqrt(2)*real(x_b_os.*exp(j*2*wc*t)); % From baseband to passband

x_p_c = clipping(x_p,CR);
X_p_c_f= fft(filter(h,1,x_p_c));

x_p_c_f = ifft(X_p_c_f);
x_b_c_f = sqrt(2)*x_p_c_f.*exp(-j*2*wc*t); % From passband to baseband

figure(1); clf
nn = (N/2-Ncp)*L + [1:N*L]; % 257~1280
nn1 = N/2*L + [-Ncp*L+1:0]; % 257~512
nn2 = N/2*L + [0:N*L]; % 512~1536

subplot(221)
plot(t(nn1)-t0, abs(x_b_os(nn1)),'b:'); % CP
hold on;
plot(t(nn2)-t0, abs(x_b_os(nn2)),'k-'); % BB signal
axis([t([nn1(1) nn2(end)])-t0 0 max(abs(x_b_os))]);
title(['Baseband signal, with CP']);
xlabel('t (normalized by symbol duration)');
ylabel('abs(x''[m])');

subplot(223) % baseband
XdB_p_os = 20*log10(abs(fft(x_b_os)));
plot(f,fftshift(XdB_p_os)-max(XdB_p_os),'k');
xlabel('frequency [Hz]');
ylabel('PSD [dB]');
axis([f([1 end]) -100 0]);

subplot(222)
% [counts(個數),centers(在這範圍內的均值)] = hist(___)
% bar(centers,counts)
[pdf_x_p,bin] = hist(x_p(nn),50);
bar(bin,pdf_x_p/sum(pdf_x_p),'k');
xlabel('x');
ylabel('pdf');
title(['Unclipped passband signal']);

subplot(224) % passband
XdB_p = 20*log10(abs(fft(x_p)));
plot(f,fftshift(XdB_p)-max(XdB_p),'k');
xlabel('frequency [Hz]');
ylabel('PSD [dB]');
axis([f([1 end]) -100 0]);

figure(2); clf
subplot(221)
[pdf_x_p_c,bin] = hist(x_p_c(nn),50);
bar(bin,pdf_x_p_c/sum(pdf_x_p_c),'k');
title(['Clipped passband signal, CR=' num2str(CR)]);
xlabel('x');
ylabel('pdf');

subplot(223)
XdB_p_c = 20*log10(abs(fft(x_p_c)));
plot(f,fftshift(XdB_p_c)-max(XdB_p_c),'k');
xlabel('frequency [Hz]');
ylabel('PSD [dB]');
axis([f([1 end]) -100 0]);

subplot(222)
[pdf_x_p_c_f,bin] = hist(real(x_p_c_f(nn)),50); 
bar(bin,pdf_x_p_c_f/sum(pdf_x_p_c_f),'k');
title(['Passband signal after clipping and filtering, CR=' num2str(CR)]);
xlabel('x');
ylabel('pdf');
axis([min(bin) max(bin) min(pdf_x_p_c_f/sum(pdf_x_p_c_f)) max(pdf_x_p_c_f/sum(pdf_x_p_c_f))]);

subplot(224)
XdB_p_c_f = 20*log10(abs(X_p_c_f));
plot(f,fftshift(XdB_p_c_f)-max(XdB_p_c_f),'k');
xlabel('frequency [Hz]'); ylabel('PSD [dB]');
axis([f([1 end]) -100 0]);

figure(3); clf
subplot(221)
stem(h,'k');
xlabel('tap');
ylabel('Filter coefficient h[n]');
axis([1, length(h), min(h), max(h)]);

subplot(222)
HdB = 20*log10(abs(fft(h,length(X_p_c_f))));
plot(f,fftshift(HdB),'k');
xlabel('frequency [Hz]');
ylabel('Filter freq response H [dB]');
axis([f([1 end]) -100 0]);

subplot(223)
[pdf_x_p_c_f,bin] = hist(abs(x_b_c_f(nn)),50);
bar(bin,pdf_x_p_c_f/sum(pdf_x_p_c_f),'k');
title(['Baseband signal after clipping and filtering, CR=' num2str(CR)]);
xlabel('|x|');
ylabel('pdf');
axis([min(bin) max(bin) min(pdf_x_p_c_f/sum(pdf_x_p_c_f)) max(pdf_x_p_c_f/sum(pdf_x_p_c_f))]);

subplot(224)
XdB_b_c_f = 20*log10(abs(fft(x_b_c_f)));
plot(f,fftshift(XdB_b_c_f)-max(XdB_b_c_f),'k');
xlabel('frequency [Hz]');
ylabel('PSD [dB]');
axis([f([1 end]) -100 0]);
