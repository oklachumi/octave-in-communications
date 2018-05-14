%channel_estimation.m
% for LS/DFT Channel Estimation with linear/spline interpolation
clear,close,clc all;

function H_LS = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt)
% LS channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       N         = FFT size
%       Nps       = Pilot spacing
%       int_opt   = 'linear' or 'spline'
% output:
%       H_LS      = LS channel etimate

Np = Nfft/Nps; % # of pilot
k = 1:Np;
LS_est(k) = Y(pilot_loc(k))./Xp(k); % LS channel estimation
if lower(int_opt(1)) == 'l',
    method = 'linear'; 
else
    method = 'spline';  
end
H_LS = interpolate(LS_est,pilot_loc,Nfft,method); % Linear/Spline interpolation
endfunction

function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR)
% function H_MMSE = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,ts,SNR)
% MMSE channel estimation function
% Inputs:
%       Y         = Frequency-domain received signal
%       Xp        = Pilot signal
%       pilot_loc = Pilot location
%       Nfft      = FFT size
%       Nps       = Pilot spacing
%       h         = Channel impulse response
%       ts        = Sampling time
%       SNR       = Signal-to-Noise Ratio[dB]
% output:
%       H_MMSE     = MMSE channel estimate

%H = fft(h,N);
snr = 10^(SNR*0.1);
Np = Nfft/Nps; % # of pilot
k = 1:Np;   
H_tilde = Y(1,pilot_loc(k))./Xp(k); % LS estimate
k = 0:length(h)-1; % k_ts = k*ts; 
hh = h*h';
tmp = h.*conj(h).*k; % tmp = h.*conj(h).*k_ts; 
r = sum(tmp)/hh;
r2 = tmp*k.'/hh; % r2 = tmp*k_ts.'/hh; A.' = transpose(A)
tau_rms = sqrt(r2-r^2); % rms delay
df = 1/Nfft; % 1/(ts*Nfft);
j2pi_tau_df = j*2*pi*tau_rms*df;
K1 = repmat([0:Nfft-1].',1,Np); % K1: Nfft*Np, each row: 0:Nfft-1
K2 = repmat([0:Np-1],Nfft,1); % K2: Nfft*Np, each row: 0:Np-1
rf = 1./(1+j2pi_tau_df*(K1-K2*Nps)); % rf[k]
K3 = repmat([0:Np-1].',1,Np); % K3: Np*Np, each column: 0:Np-1
K4 = repmat([0:Np-1],Np,1); % K4: Np*Np, each column: 0:Np-1
rf2 = 1./(1+j2pi_tau_df*Nps*(K3-K4)); % rf[k]
Rhp = rf;
Rpp = rf2 + eye(length(H_tilde),length(H_tilde))/snr;
H_MMSE = transpose(Rhp*inv(Rpp)*H_tilde.'); % MMSE channel estimate
endfunction

function H_interpolated = interpolate(H_est,pilot_loc,Nfft,method)
% Input:        H_est    = Channel estimate using pilot sequence
%           pilot_loc    = location of pilot sequence
%                Nfft    = FFT size
%              method    = 'linear'/'spline'
% Output: H_interpolated = interpolated channel

if pilot_loc(1) > 1
    slope = (H_est(2)-H_est(1))/(pilot_loc(2)-pilot_loc(1));
    H_est = [H_est(1)-slope*(pilot_loc(1)-1) H_est];
    pilot_loc = [1 pilot_loc];
end
if pilot_loc(end) < Nfft
    slope = (H_est(end)-H_est(end-1))/(pilot_loc(end)-pilot_loc(end-1));  
    H_est = [H_est H_est(end)+slope*(Nfft-pilot_loc(end))];
    pilot_loc = [pilot_loc Nfft];
end
if lower(method(1)) == 'l'
    H_interpolated = interp1(pilot_loc,H_est,[1:Nfft]);   
else
    H_interpolated = interp1(pilot_loc,H_est,[1:Nfft],'spline');
end
endfunction

figure(1); clf;
figure(2); clf;
Nfft = 32;
Ng = Nfft/8; % Ng = Add CP = 4
Nofdm = Nfft+Ng;
Nsym = 100;
Nps = 4; % Pilot spacing
Np = Nfft/Nps; % Numbers of pilots per OFDM symbol
Nd = Nfft-Np;  % Numbers of datas per OFDM symbol
Nbps = 4;
M = 2^Nbps; % Number of bits per (modulated) symbol

Es = 1;
A = sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor
%fs = 10e6;  ts = 1/fs;  % Sampling frequency and Sampling period
SNRs = [0:3:30];
sq2 = sqrt(2);
for i = 1:length(SNRs)
  SNR = SNRs(i); 
  rand('seed',1);
  randn('seed',1);
  MSE = zeros(1,6);
  nose = 0; % Number_of_symbol_errors
  for nsym = 1:Nsym
    Xp = 2*(randn(1,Np)>0)-1;    % Pilot sequence generation: randn -1 and 1
    %Data = ((2*(randn(1,Nd)>0)-1) + j*(2*(randn(1,Nd)>0)-1))/sq2; % QPSK modulation
    msgint = randint(1,Nfft-Np,M);    % bit generation
    Data = qammod(msgint,M)*A;
    %Data = modulate(mod_object, msgint); Data = modnorm(Data,'avpow',1)*Data;   % normalization
    ip = 0;
    pilot_loc = [];
    for k = 1:Nfft % 在頻域的特定位置加入導頻和數據  
      if mod(k,Nps) == 1
        X(k) = Xp(floor(k/Nps)+1);
        pilot_loc = [pilot_loc k];
        ip = ip+1;
      else
        X(k) = Data(k-ip); % ip指示了當前OFDM符號中已經加入的導頻的數量
      end
    end
    x = ifft(X,Nfft);                            % IFFT
    xt = [x(Nfft-Ng+1:Nfft) x];                  % Add CP
    h = [(randn+j*randn) (randn+j*randn)/2];     % generates a (2-tap) channel
    H = fft(h,Nfft);
    channel_length = length(h);                  % True channel and its time-domain length
    H_power_dB = 10*log10(abs(H.*conj(H)));      % True channel power in dB
    y_channel = conv(xt, h);                     % Channel path (convolution)
    sig_pow = mean(y_channel.*conj(y_channel));
    %y_aw(1,1:Nofdm) = y(1,1:Nofdm) + ...
    %   sqrt((10.^(-SNR/10))*sig_pow/2)*(randn(1,Nofdm)+j*randn(1,Nofdm)); % Add noise(AWGN)
    yt = awgn(y_channel,SNR,'measured');  
    y = yt(Ng+1:Nofdm);                   % Remove CP
    Y = fft(y);                           % FFT
    for m = 1:3
      if m == 1
        H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'linear');
        method = 'LS-linear'; % LS estimation with linear interpolation
      elseif m == 2
        H_est = LS_CE(Y,Xp,pilot_loc,Nfft,Nps,'spline');
        method = 'LS-spline'; % LS estimation with spline interpolation
      else
        H_est = MMSE_CE(Y,Xp,pilot_loc,Nfft,Nps,h,SNR);
        method='MMSE'; % MMSE estimation
      end
      H_est_power_dB = 10*log10(abs(H_est.*conj(H_est)));
      h_est = ifft(H_est);
      h_DFT = h_est(1:channel_length); 
      H_DFT = fft(h_DFT,Nfft); % DFT-based channel estimation
      H_DFT_power_dB = 10*log10(abs(H_DFT.*conj(H_DFT)));
      if nsym == 1
        figure(1);
        subplot(319+2*m); plot(H_power_dB,'b','linewidth',1); grid on; hold on;
        plot(H_est_power_dB,'r:+','Markersize',4,'linewidth',1);
        axis([0 32 -20 10]); title(method);
        xlabel('Subcarrier Index'); ylabel('Power [dB]');
        legend('True Channel',method,4); set(gca,'fontsize',10);
        subplot(320+2*m); plot(H_power_dB,'b','linewidth',1); grid on; hold on;
        plot(H_DFT_power_dB,'r:+','Markersize',4,'linewidth',1);
        axis([0 32 -20 10]); title([method ' with DFT']);
        xlabel('Subcarrier Index'); ylabel('Power [dB]');
        legend('True Channel',[method ' with DFT'],4); set(gca,'fontsize',10);
      end
      MSE(m) = MSE(m) + (H-H_est)*(H-H_est)';
      MSE(m+3) = MSE(m+3) + (H-H_DFT)*(H-H_DFT)';
    end
    Y_eq = Y./H_est;
    if nsym >= Nsym-10
      figure(2);
      subplot(121);
      plot(Y,'.','Markersize',10);
      title(['Before channel compensation']);
      %axis([-3 3 -3 3]); axis('equal'); set(gca,'fontsize',10); 
      hold on;
      subplot(122);
      plot(Y_eq,'.','Markersize',10);
      title(['After channel compensation']);
      %axis([-3 3 -3 3]); axis('equal'); set(gca,'fontsize',10);
      hold on;      
    end
    ip = 0;
    for k = 1:Nfft
      if mod(k,Nps) == 1
        ip = ip+1;
      else
        Data_extracted(k-ip) = Y_eq(k);
      end
    end
    msg_detected = qamdemod(Data_extracted/A,M);
    nose = nose + sum(msg_detected~=msgint);
  end   
  MSEs(i,:) = MSE/(Nfft*Nsym);
end   
Number_of_symbol_errors = nose
figure(3); clf;
semilogy(SNRs',MSEs(:,1),'-x', SNRs',MSEs(:,3),'-o');
xlabel('SNR [dB]'); ylabel('BER');
legend('LS-linear','MMSE');
fprintf('MSE of LS-linear/LS-spline/MMSE Channel Estimation = %6.4e/%6.4e/%6.4e\n',MSEs(end,1:3));
fprintf('MSE of LS-linear/LS-spline/MMSE Channel Estimation with DFT = %6.4e/%6.4e/%6.4e\n',MSEs(end,4:6));
