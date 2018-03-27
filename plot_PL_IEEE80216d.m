clear,close,clc all
function PL = PL_IEEE80216d(fc,d,type,htx,hrx,corr_fact,mod)
% IEEE 802.16d model
% Input - fc       : carrier frequency
%         d        : between base and terminal
%         type     : selects 'A', 'B', or 'C' 
%         htx      : height of transmitter
%         hrx      : height of receiver
%         corr_fact: if shadowing exists, set to 'ATnT' or 'Okumura'. Otherwise, 'NO'
%         mod      : set to 'mod' to get the modified IEEE 802.16d model
% output - PL      : path loss [dB]

Mod = 'UNMOD';
if nargin > 6
  Mod = upper(mod);
end
if nargin == 6 && corr_fact(1) == 'm'
  Mod = 'MOD';
  corr_fact = 'NO';  
elseif nargin < 6
  corr_fact = 'NO';
  if nargin == 5 && hrx(1) == 'm'
    Mod = 'MOD';
    hrx = 2;
  elseif nargin < 5
    hrx = 2;
    if nargin == 4 && htx(1) == 'm'
      Mod = 'MOD';
      htx = 30;
    elseif nargin < 4
      htx = 30;
      if nargin == 3 && type(1) == 'm'
        Mod = 'MOD';
        type = 'A';         
      elseif nargin < 3
        type = 'A';   
      end 
    end
  end    
end
d0 = 100;
Type = upper(type);
if Type ~= 'A' && Type ~= 'B' && Type ~= 'C'
  disp('Error: The selected type is not supported');
  return;
end

switch upper(corr_fact)
  case 'ATNT'
    Cf = 6*log10(fc/2e9);
    C_Rx = -10.8*log10(hrx/2);
  case 'OKUMURA'
    Cf = 6*log10(fc/2e9);
    if hrx <= 3
      C_Rx = -10*log10(hrx/3);  
    else
      C_Rx = -20*log10(hrx/3); 
    end
  case 'NO'
    Cf = 0;
    C_Rx = 0;
end

if Type == 'A'
  a = 4.6;
  b = 0.0075;
  c = 12.6;
elseif Type == 'B'
  a = 4;
  b = 0.0065;
  c = 17.1;
else
  a = 3.6;
  b = 0.005;
  c = 20;
end

lamda = 299792458/fc;
gamma = a-b*htx+c/htx;
d0_pr = d0;
if Mod(1) == 'M'
  d0_pr = d0*10^-((Cf+C_Rx)/(10*gamma));
end
A = 20*log10(4*pi*d0_pr/lamda)+Cf+C_Rx;

for k = 1:length(d)
  if d(k) > d0_pr
    PL(k) = A + 10*gamma*log10(d(k)/d0);
  else
    PL(k) = -10*log10((lamda/(4*pi*d(k)))^2);
  end
end
end

% plot_PL_IEEE80216d.m
fc = 2e9;
htx = [30 30];
hrx = [2 10];
distance = [1:1000];
for k = 1:2    
  y_IEEE16d(k,:) = PL_IEEE80216d(fc,distance,'A',htx(k),hrx(k),'atnt');
  y_MIEEE16d(k,:) = PL_IEEE80216d(fc,distance,'A',htx(k),hrx(k),'atnt','mod');
end
subplot(121)
semilogx(distance,y_IEEE16d(1,:),'k:','linewidth',1.5), hold on
semilogx(distance,y_IEEE16d(2,:),'k-','linewidth',1.5), grid on 
title(['IEEE 802.16d Path loss Models, f_c=',num2str(fc/1e6),'MHz'])
axis([1 1000 10 150]), xlabel('Distance[m]'), ylabel('Pathloss[dB]')
legend('h_{Tx}=30m, h_{Rx}=2m','h_{Tx}=30m, h_{Rx}=10m',2)
subplot(122)
semilogx(distance,y_MIEEE16d(1,:),'k:','linewidth',1.5), hold on
semilogx(distance,y_MIEEE16d(2,:),'k-','linewidth',1.5), grid on 
title(['Modified IEEE 802.16d Path loss Models, f_c=', num2str(fc/1e6), 'MHz'])
axis([1 1000 10 150]), xlabel('Distance[m]'), ylabel('Pathloss[dB]')
legend('h_{Tx}=30m, h_{Rx}=2m','h_{Tx}=30m, h_{Rx}=10m',2)

