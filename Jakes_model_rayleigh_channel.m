clear,clc,close all

function [h]=rayleigh(fd,t)
%該程式利用改進的jakes模型來產生單徑的平坦型瑞利衰落信道
%IEEE Commu letters, Vol.6, NO.6, JUNE 2002
%輸入變數說明：
%  fd:信道的最大多普勒頻移 單位Hz     
%   t:信號的抽樣時間序列 抽樣間隔單位s  
%   h:為輸出的瑞利信道函數 一個時間函數複序列 

    %假設的入射波數目
    N=40; 
    wm=2*pi*fd;
    %每象限的入射波數目即振盪器數目
    N0=N/4;
    %信道函數的實部
    Tc=zeros(1,length(t));
    %信道函數的虛部
    Ts=zeros(1,length(t));
    %歸一化功率係數
    P_nor=sqrt(1/N0);
    %區別個條路徑的均勻分佈隨機相位
    theta=2*pi*rand(1,1)-pi;
    for ii=1:N0
            %第i條入射波的入射角 
            alfa(ii)=(2*pi*ii-pi+theta)/N;
            %對每個子載波而言在(-pi,pi)之間均勻分佈的隨機相位
            fi_tc=2*pi*rand(1,1)-pi;
            fi_ts=2*pi*rand(1,1)-pi;
            %計算衝激響應函數
            Tc=Tc+cos(cos(alfa(ii))*wm*t+fi_tc);
            Ts=Ts+cos(sin(alfa(ii))*wm*t+fi_ts);
    end;
%乘歸一化功率係數得到傳輸函數
h=P_nor*(Tc+j*Ts);
endfunction

fd=10;              %Doppler shift 10
ts=1/1000;          %信道sampling time
t=0:ts:1;           
h1=rayleigh(fd,t);  

fd=20;              
h2=rayleigh(fd,t);

subplot(2,1,1),plot(20*log10(abs(h1(1:1/ts))));
title('fd =10Hz power vs t');
xlabel('time');ylabel('power');
subplot(2,1,2),plot(20*log10(abs(h2(1:1/ts))));
title('fd=20Hz power vs t');
xlabel('time');ylabel('power');

