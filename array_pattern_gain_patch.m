%ArrayPatGain - This program is an updated and combined version of the programs "ANTGAIN
% and "ARRAYPATGAIN" in the 4th edition of Kraus's Electromagnetics.  This program
% computes and plots the field pattern of a uniform linear array of sources.  
clc, clear;
timesrun=0;
while timesrun<1000,
    if timesrun==0,
        SP=1;
	      SP=0.49/sqrt(2.27);
        PH=0;
        N=2;
        MF=1;
    else,
        SP=input('Enter element spacing in wavelengths: ');
        PH=input('Enter phase difference between elements in degrees: ');
        N=input('Enter number of elements: ');
        MF=input('Enter pattern multiplication factor: ')';
    end;
   
% This version modified to include isotropic AND short dipole sources
% both of which are displayed on the same polar plot.  This version
% does not do phase plots.


% Compute fields for plotting
    A=0.01:0.01:6.27;                     % Angle over 2pi radians
    B=0.01:0.01:3.14;                     % B is theta -- the angle from 0 degrees
    
    U=(2*pi*SP*cos(A)+(pi*PH/180))/2;     % compute argument
    FP=(sin(N*U)./sin(U));                % compute isotropic electric field
    FP_n=FP./max(FP);			  % nomalized the electric field
    %FPD=(sin(N*U)./sin(U)).*cos(A);      % computer dipole electric field
    FPD=(cos((2*pi*SP*cos(B))/2.)-cos(2*pi*SP/2.))./sin(B);
    FPD_n=FPD./max(FPD);		  % normalized the electric field
    R=MF.*abs(FP_n);                      % multiply as appropriate
    Rdipole=MF.*abs(FPD_n);               % compute for array of short dipoles    
    pause(1);

% For isotropic source    
% compute the beam area in theta

                                          %    0 <theta (B) < 180 degrees
    W=(2*pi*SP*cos(B)+(pi*PH/180))/2;     % compute psi -- phase shift
    PP=(sin(N*W)./sin(W)).^2;             % compute unnormalized power
    Z=0.01*sin(B).*PP;                    % differential power is PP*sin(theta)*d(theta)
    SUM=sum(Z);                           % integrate the elements over 180 degrees

    DR=(2*(N^2))/SUM;                     % numerical directivity is Pmax/Pavg.
                                          % N^2 is Pmax 
                                          % D = 4*pi/omegaA                                                                                    
    DBI=10*log10(DR);                     % express in dBi

% For dipole source
% compute the beam area in theta
%   B=0.01:0.01:3.14;                     % B is theta -- the angle from 0 degrees
                                          %    0 <theta (B) < 180 degrees
    W=(2*pi*SP*cos(B)+(pi*PH/180))/2;     % compute psi -- phase shift
    PP=(sin(N*W)./sin(W).*cos(B)).^2;     % compute unnormalized power
    PP_d=((cos((2*pi*SP.*cos(B))/2.)-cos(2*pi*SP/2.))./sin(B)).^2;
    Z=0.01*sin(B).*PP;                    % differential power is PP*sin(theta)*d(theta)
    Z_d=0.01*sin(B).*PP_d;
    SUM=sum(Z);                           % integrate the elements over 180 degrees
    SUM_d=sum(Z_d);
    
    DRdipole=2./SUM_d;                    %(2*(N^2))/SUM; % numerical directivity is Pmax/Pavg.
                                          % N^2 is Pmax 
                                          % D = 4*pi/omegaA                                          
    DBIdipole=10*log10(DRdipole);         % express in dBi
    theta=(0:2*pi/626:2*pi);
    theta_d=(0:pi/313:pi);             

% plot antenna pattern in polar coordinates
% plot isotropic pattern and then plot dipole pattern in same fig.
    polar(theta,(round(R*100))/100,'r'); 
    set (gca, 'rtick', 0:0.2:1, 'ttick', 0:10:350);
    hold on;
    polar(theta_d,(round(Rdipole*100))/100,'b');
    set (gca, 'rtick', 0:0.2:1, 'ttick', 0:10:350);
    title(['Field pattern of ',num2str(N),' sources spaced by ',num2str(SP),' lambda,  ','phase delta= ',num2str(PH),' deg']);
    text(max(R),max(R),['Dtheta = ',num2str((round(DBI*100))/100),' dBi']);
    text(max(R),-max(R),['Dphi = ',num2str((round(DBIdipole*100))/100),' dBi']);
    hold off;

    sep=max(R)/20;
    xcent=sep*((-(N-1)/2)-1);
    for elcount=1:N
        xcent=xcent+sep;
        ycent=0;
        radius=sep/5;
        xp=xcent+[-radius:radius/10:radius];
        yp=ycent+real(sqrt((radius^2)-(xp-xcent).^2));
        patch(xp,yp,'k');
        patch(xp,-yp,'k');
    end;

    rp=input('Enter 1 for rectangular plot: ');
    if rp==1
        figure(2);
        plot(theta*(180/pi),R);
        grid on;
        %zoom on;
        title(['Field Pattern of ',num2str(N),' sources spaced by ',num2str(SP),' lambda,  ','phase delta= ',num2str(PH),' deg']);
        xlabel('Angle (deg)');
        ylabel('Amplitude (linear)');
    end;

    logp=input('Enter 1 for dB plot: ');
    if logp==1
        figure(3);
        plot(theta*(180/pi),10*log10(abs(R)));
        grid on;
        title(['Field Pattern of ',num2str(N),' sources spaced by ',num2str(SP),' lambda,  ','phase delta= ',num2str(PH),' deg']);
        xlabel('Angle (deg)');
        ylabel('Amplitude (dB)');
    end;

    powerp=input('Enter 1 for power pattern: ');
    if powerp==1
        figure(4);
        plot(theta*(180/pi),20*log10(abs(R)));
        grid on;
        title(['Power Pattern of ',num2str(N),' sources spaced by ',num2str(SP),' lambda,  ','phase delta= ',num2str(PH),' deg']);
        xlabel('Angle (deg)');
        ylabel('Amplitude (dB)');
    end;

    timesrun=timesrun+1;
    done=input('Enter 1 to modify parameters: ');
    if done~=1,
        timesrun=1000;
    end;

end;
