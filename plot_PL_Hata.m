clear,close,clc all

function PL = PL_Hata(fc,d,htx,hrx,Etype)
% Hata Model
% Input
%       fc    : carrier frequency [Hz]
%       d     : between base station and mobile station [m]
%       htx   : height of transmitter [m]
%       hrx   : height of receiver [m]
%       Etype : Environment Type('urban','suburban','open')
% output
%       PL    : path loss [dB]

if nargin < 5
    Etype = 'URBAN';
end
fc = fc/(1e6);
if fc >= 150 && fc <= 200
    C_Rx = 8.29*(log10(1.54*hrx))^2 - 1.1;
elseif fc > 200
    C_Rx = 3.2*(log10(11.75*hrx))^2 - 4.97;
else
    C_Rx = 0.8+(1.1*log10(fc)-0.7)*hrx-1.56*log10(fc);
end
PL = 69.55+26.16*log10(fc)-13.82*log10(htx)-C_Rx...
    +(44.9-6.55*log10(htx))*log10(d/1000);
EType = upper(Etype);
if EType(1) == 'S'
    PL = PL-2*(log10(fc/28))^2-5.4;
elseif EType(1) == 'O' 
    PL = PL+(18.33-4.78*log10(fc))*log10(fc)-40.97;
end
end

% plot_PL_Hata.m
fc = 1.5e9;
htx = 30;
hrx = 2;
distance = [1:2:31].^2; 
y_urban = PL_Hata(fc,distance,htx,hrx,'urban');
y_suburban = PL_Hata(fc,distance,htx,hrx,'suburban');
y_open = PL_Hata(fc,distance,htx,hrx,'open');
semilogx(distance,y_urban,'k-s', distance,y_suburban,'k-o', distance,y_open,'k-^');
grid on
axis([1 1000 40 110]);
title(['Hata PL model, f_c=', num2str(fc/1e6), 'MHz']);
xlabel('Distance [m]');
ylabel('Path loss [dB]');
legend('urban','suburban','open area','location','northwest');
