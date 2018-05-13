% plot_Ray_Ric_channel.m
clear,close,clc all

function H = Ray_model(L)
% Rayleigh Channel Model
% Input : L  : # of channel realization
% Output: H  : Channel vector

H = (randn(1,L)+j*randn(1,L))/sqrt(2);
endfunction

function H=Ric_model(K_dB,L)
% Rician Channel Model
%   Input:
%       K_dB   : K factor [dB]
%       L      : # of channel realization
%   Output:
%       h      : channel vector

K=10^(K_dB/10);
H = sqrt(K/(K+1)) + sqrt(1/(K+1))*Ray_model(L);
endfunction

N = 200000;
level = 30;
K_dB = [-40 15 30];
Rayleigh_ch = zeros(1,N);
Rician_ch = zeros(2,N);
color = ['k'];
line = ['-'];
marker = ['s','o','^','*'];

% Rayleigh model
% [counts,centers] = hist(___)
Rayleigh_ch = Ray_model(N); 
[temp,x] = hist(abs(Rayleigh_ch(1,:)),level);   
plot(x,temp,['k-' marker(1)]);
hold on

% Rician model
for i = 1:length(K_dB);
  Rician_ch(i,:) = Ric_model(K_dB(i),N);
  [temp x] = hist(abs(Rician_ch(i,:)),level);   
  plot(x,temp,['k-' marker(i+1)]);
end
xlabel('x');
ylabel('Occurance');
legend('Rayleigh','Rician, K = -40 dB','Rician, K = 15 dB','Rician, K = 30 dB');
