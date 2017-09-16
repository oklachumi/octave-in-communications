clear,clc,close all
function [h]=rayleigh(fd,t)
%該程式利用改進的jakes模型來產生單徑的平坦型Rayleigh衰落信道
%IEEE Commu letters, Vol.6, NO.6, JUNE 2002
%輸入變數說明：
%  fd:信道的最大多普勒頻移 單位Hz     
%   t:信號的抽樣時間序列 抽樣間隔單位s  
%   h:為輸出的Rayleigh信道函數 一個時間函數複序列 

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
  for i=1:N0
    %第i條入射波的入射角 
    alfa(i)=(2*pi*i-pi+theta)/N;
    %對每個子載波而言在(-pi,pi)之間均勻分佈的隨機相位
    fi_tc=2*pi*rand(1,1)-pi;
    fi_ts=2*pi*rand(1,1)-pi;
    %計算衝激響應函數
    Tc=Tc+cos(cos(alfa(i))*wm*t+fi_tc);
    Ts=Ts+cos(sin(alfa(i))*wm*t+fi_ts);
  end;
%乘歸一化功率係數得到傳輸函數
h=P_nor*(Tc+j*Ts);
endfunction

function [DATA]=intdump(IN,num)
  outidx=1;
  for z=1:num:length(IN)
    DATA(outidx)=sum(IN(z:z+num-1))/num;
    outidx=outidx+1;
  end
% return DATA
endfunction

nsamp = 8;              %矩形脈衝取樣點數
numsymb = 10000;        %每種SNR下的傳輸的符號數
ts=1/(numsymb*nsamp);
t=(0:numsymb*nsamp-1)*ts;

M=4;                    
SNR=-3:3;               
grayencod=[0 1 3 2];    
for i=1:length(SNR)
  msg=randsrc(1,numsymb,[0:3]);                   %產生發送符號
  %msg_gr=grayencod(msg+1);                       %進行Gray編碼映射
  msg_tx=pskmod(msg,M,pi/M,'gray');               %QPSK調製
  msg_tx=msg_tx(ones(nsamp,1),:)(:).';            %rectpulse shaping
  h=rayleigh(10,t);                               %生成Rayleigh衰落
  msg_tx_rale=h.*msg_tx;                          %信號通過Rayleigh衰落信道
  msg_rx=awgn(msg_tx,SNR(i),'measured');          %通過AWGN信道
  msg_rx_rale=awgn(msg_tx_rale,SNR(i),'measured');%Rayleigh + AWGN
  msg_rx=intdump(msg_rx,nsamp);                   %rectpulse shaping integral&dump
  msg_rx_rale=intdump(msg_rx_rale,nsamp);         %rectpulse shaping integral&dump  
  msg_demod=pskdemod(msg_rx,M,pi/M,'gray');                 %QPSK解調
  msg_demod_rale=pskdemod(msg_rx_rale,M,pi/M,'gray');       %QPSK解調
  [errorBit BER(i)]          =biterr(msg,msg_demod,log2(M));%計算BER
  [errorBit_rale BER_rale(i)]=biterr(msg,msg_demod_rale,log2(M));
  [errorSym SER(i)]          =symerr(msg,msg_demod);        %計算SER
  [errorSym_rale SER_rale(i)]=symerr(msg,msg_demod_rale); 
end

scatterplot(msg_tx);                              
title('AWGN Tx constellation');
xlabel('I');
ylabel('Q');
scatterplot(msg_rx);                              
title('AWGN Rx constellation');
xlabel('I');
ylabel('Q');
scatterplot(msg_tx_rale);                              
title('AWGN+Rayleigh Tx constellation');
xlabel('I');
ylabel('Q');
scatterplot(msg_rx_rale);                              
title('AWGN+Rayleigh Rx constellation');
xlabel('I');
ylabel('Q');

figure
semilogy(SNR,BER,'-ro',SNR,SER,'-r*',SNR,BER_rale,'-b+',SNR,SER_rale,'-b^');
legend('AWGN\_BER','AWGN\_SER','Rayleigh+AWGN\_BER','Rayleigh+AWGN\_SER','location','southwest');
title('QPSK in AWGN+Rayleigh fading channel capability');
xlabel('SNR(dB)');
ylabel('BER and SER');
