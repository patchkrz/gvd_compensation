%
% PROPAGATION METHOD FOR THE SOLUTION OF NLSE
%
% i Fz - beta2/2 Ftt + gamma*|F|^2 F=0
%

clc
clear all
close all

% MATERIAL PROPERTIES

% beta21=-1;
% beta22=-beta21;

beta2 = -1;
beta2_cpn = -1*beta2;
% beta2 = 1;
gamma=0; % making this we are sure that we focused on self-phase modulation effect is discarded.

% TEMPORAL COORDINATE
t0=-150;
t1=150;
nt=5000;
t=linspace(t0,t1,nt);
deltat=(t1-t0)/(nt-1);

% SPATIAL COORDINATE
z0=0;
z1=5;
nz=1000;
z=linspace(z0,z1,nz);
deltaz=(z1-z0)/(nz-1); % deltaz = h
h=deltaz;

% FREQUENCY COORDINATE 
indfreq=-nt/2:1:nt/2-1;
omega=(pi./t1).*indfreq;

% INPUT ENVELOPE
% tau0=1;
% amp=1;
% FIN=amp*exp(-t.^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE
% C=1;
% FING=amp*exp(-(1+i*C)/2*t.^2/tau0^2);  % CHIRPED PULSES
% FIN=FING;

% C=-2;
% tau0=1;
% FIN=sech(t/tau0).*exp(-i*C*t.^2/(2*tau0^2)); % HYPERBOLIC CHIRPED SECH PULSES

% C=0;
% m=4;
% FIN=exp(-(1+i*C)/2*(t/tau0).^(2*m));

% 5 gaussian pulses 1 0 1 1 1 0 0 1
tau0=1; 
amp=1;

% tshift0=-40;
% FIN0=amp*exp(-(t+tshift0).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 0
% 
% tshift1=-10;
% FIN1=amp*exp(-(t+tshift1).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 1

tshift2=0;
FIN2=1*amp*exp(-(t+tshift2).^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE 2

% tshift3=10; 
% FIN3=amp*exp(-(t+tshift3).^2/(2*tau0^2)); % UNCHIRPEDGAUSSIAN PULSE 3
% 
% tshift4=30; 
% FIN4=amp*exp(-(t+tshift4).^2/(2*tau0^2)); % UNCHIRPEDGAUSSIAN PULSE 4

FIN=FIN2; % +FIN0+FIN1+FIN3+FIN4
% FIN=amp*exp(-t.^2/(2*tau0^2)); % UNCHIRPED GAUSSIAN PULSE

%FOURIER TRANSFORM
%FFING=deltat*fftshift(fft(FING));
FFIN=deltat*fftshift(fft(FIN));

% figure
% plot(omega,abs(FFIN),omega,abs(FFING),'r')
% xlabel('\omega')
% ylabel('|Fourier Transform F|')
% grid
% set(findall(gcf,'type','text'),'FontSize',30)
% set(gca,'FontSize',30)
% stopprof

%%%%%%%%%%%%%%%%%%%%%%
% CORE OF THE PROGRAM
%%%%%%%%%%%%%%%%%%%%%%

F=zeros(ceil(nt),floor(nz));
FF=zeros(ceil(nt),floor(nz));

F(:,1)=FIN;
FF(:,1)=FFIN;

q=FIN;

for loop_step=2:1:nz
    
    % LINEAR DISPERSIVE STEP
    qs=deltat*fftshift(fft(ifftshift(q))); % FFT
    
    qs_old=qs;
    
    prop=beta2/2*omega.^2;

    if loop_step < (nz/2)-1  
    fact = 1j*beta2/2*omega.^2*h;
    qs   = qs_old.*exp(fact);
    q    = (1/deltat)*fftshift(ifft(ifftshift(qs))); 
    % calc of the envelope under dispersive effect
    else
    fact = 1j*beta2_cpn/2*omega.^2*h;
    qs   = qs_old.*exp(fact);
    q    = (1/deltat)*fftshift(ifft(ifftshift(qs))); 
    end

    % fact=j*prop*h;
    % qs=qs_old.*exp(fact); % calculation of the propagation in the frequency domain
    % q=(1/deltat)*fftshift(ifft(ifftshift(qs))); % coming back in the time domain
    
    %NONLINEAR CHI3 EFFECT STEP
    q_old=q;
    q=q_old.*exp(j*gamma*abs(q_old).^2*h);
    
    %SAVE DATA EVERY NZ
    F(:,loop_step)=q;
    FF(:,loop_step)=deltat*fftshift(fft(ifftshift(q)));
        
end

save dati

deltaOmega_in = gradient(unwrap(angle(F(:,1))),t); 
deltaOmega_out = gradient(unwrap(angle(F(:,end))),t); 

plot(t,deltaOmega_in,t,deltaOmega_out)
xlabel('t')
ylabel('\delta \omega')
title('Induced Chirp For Single Pulse After Dispersion Compensation')
legend('IN','OUT')
set(gca,'FontSize',12)
grid

% subplot(2,2,1)
% plot(t,abs(FIN),t,abs(F(:,end)))
% xlabel('t')
% ylabel('|F|')
% title('Input-Output Envelope')
% legend('IN','OUT')
% set(gca,'FontSize',12)
% grid

% subplot(2,2,2)
% phase = [angle(F(1:length(nt)/2,end)') angle(F((length(nt)/2)+1:end,end)')];
% plot(t,angle(F(:,1)),t,unwrap(angle(F(:,end))))
% xlabel('t')
% ylabel('angle F')
% title('Input-Output Phase Comparision')
% legend('IN','OUT')
% grid

% mesh(t,z,unwrap(angle(F)')); 
% xlabel('t')
% ylabel('z')
% zlabel('angle F')
% title('Dispersion Compensation Effect On Phase ')
% xlim([-9 9]); 
% set(gca,'FontSize',12)
% grid on

% subplot(2,2,3)
% plot(omega,abs(FFIN),omega,abs(FF(:,end)))
% xlabel('\omega')
% ylabel('|FF|')
% legend('IN','OUT')
% title('Frequency Content Comparision')
% grid
 
% subplot(2,2,4)
% plot(omega,angle(FIN),omega,unwrap(angle(FF(:,end)')))
% xlabel('\omega')
% ylabel('angle(FF)')
% legend('IN','OUT')
% title('Phase Change Comparision')
% %xlim([-100 100]); ylim([0 100]);
% set(gca,'FontSize',12)
% grid

% subplot(1,2,1)
% mesh(t,z,abs(F)')
% xlabel('t')
% ylabel('z')
% zlabel('|F|')
% title('Pulse Propagation')
% xlim([-150 150]); ylim([0 5]);
% set(gca,'FontSize',12)
% grid
% 
% subplot(1,2,2)
% mesh(t,z,abs(F)')
% view(2)
% xlim([-150 150]);
% xlabel('t')
% ylabel('z')
% zlabel('|F|')
% title('Pulse Propagation')
% set(gca,'FontSize',12)
% grid
 
% subplot(1,2,2)
% mesh(omega,z,abs(FF)')
% xlabel('\omega')
% ylabel('z') 
% zlabel('|FFT F|')
% title('Frequency Content Through Propagation')
% set(gca,'FontSize',12)
% grid

% plot(t,abs(F(:,1)),t,abs(F(:,end)),'r:',t,abs(F(:,50)))
% legend('IN','OUT','MIDDLE')
% subplot(2,1,1)
% plot(t,abs(F(:,1)))
% title('Envelope Of Input Field')
% xlabel('t')
% ylabel('|F|')
% grid
 
% subplot(2,1,2)
% plot(t,abs(F(:,end)))
% title('Envelope Of Output Field')
% xlabel('t')
% ylabel('|F|')
% grid
% set(findall(gcf,'type','text'),'FontSize',12)
% set(gca,'FontSize',12)

% figure
% mesh(t,z,unwrap(angle(F))')
% xlabel('t')
% ylabel('z')
% zlabel('angle F')
 
% figure
% plot(t,angle(F(:,1)),t,unwrap(angle(F(:,end))))
% xlabel('t')
% ylabel('angle F')
% legend('IN','OUT')
% xlim([-100 100])
% set(gca,'FontSize',12)
% title('Input-Output Phase Comparision')
% grid

% figure
% mesh(omega,z,abs(FF)')
% xlabel('\omega')
% ylabel('z') 
% zlabel('|FFT F|')
% set(findall(gcf,'type','text'),'FontSize',12)
% set(gca,'FontSize',12) 

% subplot(2,1,1)
% plot(omega,abs(FF(:,1)),omega,abs(FF(:,end)))
% xlabel('\omega')
% ylabel('|FFT F|')
% legend('IN','OUT')
% title('Frequency Content For 8-bits Case Through Propagation')
% xlim([-4 4])
% set(findall(gcf,'type','text'),'FontSize',12)
% set(gca,'FontSize',12)
% grid
 
% subplot(2,1,2)
% plot(omega,abs(FF(:,length(FF)/2)),'r')
% xlabel('\omega')
% ylabel('|FFT F|')
% legend('MIDDLE')
% xlim([-4 4])
% set(findall(gcf,'type','text'),'FontSize',12)
% set(gca,'FontSize',12)
% grid
