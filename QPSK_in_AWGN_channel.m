clear,close,clc all
function [DATA]=intdump(IN,num)
  outidx=1;
    for z=1:num:length(IN)
      DATA(outidx)=sum(IN(z:z+num-1))/num;
      outidx=outidx+1;
    end
% return DATA
end
M=4;                                              %QPSK的符號類型
nsamp=8;
numsymb=1e5;                                      %每種SNR下的傳輸的符號數
SNR=-3:3;
grayencod=[0 1 3 2];                              %Gray編碼格式
for i=1:length(SNR)
  msg=randsrc(1,numsymb,[0:3]);                   %產生發送符號
  %msg_gr=grayencod(msg+1);                       %進行Gray編碼映射
  msg_tx=pskmod(msg,M,pi/M,'gray');               %QPSK調製
  msg_tx=msg_tx(ones(nsamp,1),:)(:).';            %rectpulse shaping
  msg_rx=awgn(msg_tx,SNR(i),'measured');          %通過AWGN信道
  msg_rx=intdump(msg_rx,nsamp);                   %rectpulse shaping integral&dump
  msg_demod=pskdemod(msg_rx,M,pi/M,'gray');       %QPSK解調
  %[dummy graydecod]=sort(grayencod);             %Gray編碼逆映射
  %graydecod=graydecod-1;                         %Gray編碼逆映射
  %msg_demod=graydecod(msg_demod+1);              %Gray編碼逆映射
  [errorBit BER(i)]=biterr(msg,msg_demod,log2(M));%計算BER
  [errorSym SER(i)]=symerr(msg,msg_demod);        %計算SER
endfor

scatterplot(msg_tx);                              %畫出發射信號的星座圖
title('Tx constellation');
xlabel('I');
ylabel('Q');
scatterplot(msg_rx);                              %畫出接收信號的星座圖
title('Rx constellation');
xlabel('I');
ylabel('Q');
figure
semilogy(SNR,BER,'-ro',SNR,SER,'-r*');            %畫出BER和SNR隨SNR變化的曲線
legend('BER','SER');
title('QPSK in AWGN channel');
xlabel('SNR(dB)');
ylabel('BER and SER');