clear,close,clc all
snr=10;      %信噪比，單位為dB 
fftl=128;    %FFT的長度
N=6;         %一個幀結構中OFDM信號的個數
para=128;    %平行傳輸的子載波個數
gsl=32;      %保護時隙的長度
%**************OFDM信號的產生***************
signal=randint(1,para*N*2);
%產生0、1隨機序列，符號個數為para*N*2（子通道數*調製水平*每個子通道中符號數）
for i=1:para
    for j=1:N*2
        sigpara(i,j)=signal(i*j);
        %串並變換，將隨機產生的二進位矩陣變換為列數為para，行數為N*2的矩陣
    end
end
%*******以下進行QPSK調製，將資料分為I、Q兩路********
for j=1:N;
    ich(:,j)=sigpara(:,2*j-1);
    qch(:,j)=sigpara(:,2*j);
end
kmod=1./sqrt(2);
ich1=ich.*kmod;
qch1=qch.*kmod;
x=ich1+qch1.*I;            %產生complex信號
y=ifft(x);                 %IFFT將頻域信號轉換為時域信號
ich2=real(y);              %I通道取時域信號的實部
qch2=imag(y);              %Q通道取時域信號的虛部
%**********以下插入保護時隙**********
ich3=[ich2(fftl-gsl+1:fftl,:);ich2];
qch3=[qch2(fftl-gsl+1:fftl,:);qch2];
%**********以下進行並串變換**********
ich4=reshape(ich3,1,(fftl+gsl)*N);
qch4=reshape(qch3,1,(fftl+gsl)*N);
%以下為系統發送端形成的信號
Tdata=ich4+qch4.*I;
%**********以下為系統接收端進行解調的過程**********
Rdata=awgn(Tdata,snr,'measured');   %對接收到的信號加入AWGN
%**********以下為接收端移去保護時隙**********
idata=real(Rdata);
qdata=imag(Rdata);
idata1=reshape(idata,fftl+gsl,N);
qdata1=reshape(qdata,fftl+gsl,N);
idata2=idata1(gsl+1:gsl+fftl,:);
qdata2=qdata1(gsl+1:gsl+fftl,:);
%**********以下為系統接收端進行FFT**********
Rx=idata2+qdata2*I;
ry=fft(Rx);
Rich=real(ry);
Rqch=imag(ry);
Rich=Rich/kmod;
Rqch=Rqch/kmod;
%**********以下為接收端進行QPSK解調**********
for j=1:N;
    Rpara(:,2*j-1)=Rich(:,j);
    Rpara(:,2*j)=Rqch(:,j);
end
Rsig=reshape(Rpara,1,para*N*2);
Rsig=Rsig > 0.5;   %抽樣判決

[errorBit BER]=biterr(Rsig(1:128),signal(1:128))      %計算BER

figure(1)
subplot(2,1,1)
stem(Rsig(1:128))
ylabel('demodulation data');
grid on;
subplot(2,1,2)
stem(signal(1:128))
ylabel('origional data');
grid on;	 
