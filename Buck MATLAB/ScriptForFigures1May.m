%% Script for MS Thesis Figures

presentJournal(0)
clear
close all
clc
%%   
% Color code for each undersampling factor
plotcolor.ULA       = [1 0 0];          % ULA: red
plotcolor.Marray    = [0 .25 .9];       % Marray: blue
plotcolor.Narray    = [0 .75 0];      % Narray: green
% plotcolor.CSACBF    = [1 .6 0];         % CSACBF: orange
plotcolor.CSACBF    = [0 0.5 0.5];         % CSAcbf: orange
plotcolor.CSAPP     = [0 0 0];          % CSAPP: purple
plotcolor.sim       = [.25 .25 .25];    % sim: almost black

uNfft = 2048;
ur = linspace(0,1,uNfft/2);
ul = linspace(-1,0-1/uNfft,uNfft/2);
u = [ul ur];

%% Chapter 2.1 ULA
% ULA Beampattern
L = 30;
v = ones(L,1);
wULA = v/L;
WULA = fftshift(fft(wULA,uNfft));

f1 = figure(1);
set(f1,'DefaultAxesFontSize',24)
clf
plot(u,20*log10(abs(WULA)),'color',plotcolor.ULA)
ylim([-35 0])
title(['ULA Beampattern, L = ' num2str(L)])
xlabel('u = sin($\theta$)'); ylabel('[dB]')

%% ULA Coarray
Iula = v';
kappaULA = conv(Iula,fliplr(Iula));
gamma = 1-L:L-1;

f2 = figure(2);
set(f2,'DefaultAxesFontSize',24)
clf
stem(gamma,kappaULA,'MarkerSize',10,'color',plotcolor.ULA)
title('ULA Coarray, $\kappa_{ula}[\gamma]$')
xlabel('$\gamma$'); ylabel('Redundancies')

%% M-array Coarray
M = 2;
I_M = zeros(1,L);
I_M(1:M:end) = 1;
kappaM = conv(I_M,fliplr(I_M));

f3 = figure(3);
set(f3,'DefaultAxesFontSize',24)
clf
stem(gamma,kappaM,'MarkerSize',10,'color',plotcolor.Marray)
title('$M$-array Coarray, $\kappa_{M}[\gamma]$')
xlabel('$\gamma$'); ylabel('Redundancies')

% M-array Beampattern
wM = v/L*M.*I_M.';
WM = fftshift(fft(wM,uNfft));

%% Increased Undersampling Beampatterns
N = 3;
I_N = zeros(1,L);
I_N(1:N:end) = 1;
wN = v/L*N.*I_N.';
WN = fftshift(fft(wN,uNfft));

USF = 5;
I_USF = zeros(1,L);
I_USF(1:USF:end) = 1;
wUSF = v/L*USF.*I_USF.';
WUSF = fftshift(fft(wUSF,uNfft));

f4 = figure(4);
set(f4,'DefaultAxesFontSize',24)
clf
subplot(311)
plot(u,20*log10(abs(WM)),'color',plotcolor.Marray)
ylim([-35 0])
title('Effects of Increased Undersampling')

subplot(312)
plot(u,20*log10(abs(WN)),'color',plotcolor.Narray)
ylim([-35 0])
ylabel('[dB]')

subplot(313)
plot(u,20*log10(abs(WUSF)),'m')
ylim([-35 0])
xlabel('u')

%% Chapter 2.2 CSA
% Basic CSA 
L = 6;
v = ones(L,1);
wULA = v/L;
WULA = fftshift(fft(wULA,uNfft));

% Basic undersampled array sampling functions
M = 2;
I_M = zeros(1,L);
I_M(1:M:end) = 1;

N = 3;
I_N = zeros(1,L);
I_N(1:N:end) = 1;
I_Lambda = double(or(I_M,I_N));

% Basic undersampled beampattern
wM = v/L*M.*I_M.';
WM = fftshift(fft(wM,uNfft));
wN = v/L*N.*I_N.';
WN = fftshift(fft(wN,uNfft));

kappaMN = conv(I_M,fliplr(I_N));
KMN = fftshift(fft(kappaMN/L,uNfft));

f5 = figure(5);
set(f5,'DefaultAxesFontSize',24)
clf
subplot(211)
hold on
plot(u,20*log10(abs(WM)),'color',plotcolor.Marray)
plot(u,20*log10(abs(WN)),'color',plotcolor.Narray)
ylim([-25 0])
title(['Full ULA vs. CSA$_{\mathrm{pp}}$: L = ' num2str(L) ', M = ' num2str(M) ', N = ' num2str(N)])
legend('M','N','location','southeast')
ylabel('[dB]')

subplot(212)
hold on
plot(u,20*log10(abs(WULA)),'color',plotcolor.ULA)
plot(u,10*log10(abs(KMN)),'color',plotcolor.CSAPP)
ylim([-25 0])
xlabel('u'); ylabel('[dB]')
legend('ULA','CSA$_{\mathrm{pp}}$','location','southeast')


%% Extended aperture analysis
P = 3;
Le = M*N*P;
gamma = 1-Le:Le-1;
v = ones(Le,1);
kappaULA = conv(v,flipud(v));
wULAe = v/sqrt(Le)^2;
WULAe = fftshift(fft(wULAe,uNfft));

I_Me = repmat(I_M,[1 P]);
I_Ne = repmat(I_N,[1 P]);
I_Lambdae = repmat(I_Lambda,[1 P]);
Lambda = sum(I_Lambdae);

wMe = v/Le*M.*I_Me.';
WMe = fftshift(fft(wMe,uNfft));
wNe = v/Le*N.*I_Ne.';
WNe = fftshift(fft(wNe,uNfft));

kappaMNe = conv(I_Me,fliplr(I_Ne));
KMNe = fftshift(fft(kappaMNe/P/Le,uNfft));

f6 = figure(6);
set(f6,'DefaultAxesFontSize',24)
clf
subplot(211)
hold on
plot(u,20*log10(abs(WMe)),'color',plotcolor.Marray)
plot(u,20*log10(abs(WNe)),'color',plotcolor.Narray)
ylim([-35 0])
title(['Full ULA vs. CSA$_{\mathrm{pp}}$: L = ' num2str(Le) ', M = ' num2str(M) ', N = ' num2str(N) ', $\varepsilon$ = ' num2str(P)])
legend('M','N','location','southeast')
ylabel('[dB]')

subplot(212)
hold on
plot(u,20*log10(abs(WULAe)),'color',plotcolor.ULA)
plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSAPP)
ylim([-35 10*log10(Le)])
xlabel('u'); ylabel('[dB]')
legend('ULA','CSA$_{\mathrm{pp}}$','location','southeast')

f7 = figure(7);
set(f7,'DefaultAxesFontSize',24)
clf
stem(gamma,kappaMNe,'x','MarkerSize',10,'color',plotcolor.CSAPP)
% title('CSA$_{\mathrm{pp}}$ coarray, $\kappa^{pp}_{csa}[\gamma]$')
ylim([0 5.01])
xlabel('$\gamma$'); ylabel('Redundancies')

%% CSA_CBF
kappaLambda = conv(I_Lambdae,fliplr(I_Lambdae));
KLambda = fftshift(fft(kappaLambda/Lambda^2,uNfft));

f8 = figure(8);
set(f8,'DefaultAxesFontSize',50)
clf

% stem(gamma,kappaULA,'^','MarkerSize',12,'color',plotcolor.ULA)
% hold on
stem(gamma,kappaLambda,'+','MarkerSize',30,'color',plotcolor.CSACBF)
hold on
stem(gamma,kappaMNe,'x','MarkerSize',30,'color',plotcolor.CSAPP)
% title('CSA$_{\mathrm{cbf}}$ coarray, $\kappa^{cbf}_{csa}[\gamma]$')
% ylim([0 Lambda+.01])
xlabel('$\gamma$'); ylabel('Redundancies')
ylim([0 13])
legend('$\kappa_{csa}^{cbf}[\gamma]$','$\kappa_{csa}^{pp}[\gamma]$')

f9 = figure(9);
set(f9,'DefaultAxesFontSize',40)
clf
hold on
plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSAPP)
plot(u,10*log10(abs(KLambda)),'color',plotcolor.CSACBF)
ylim([-35 10*log10(Le)])
title(['CSA$_{\mathrm{pp}}$ vs. CSA$_{\mathrm{cbf}}$ Beampatterns: L = ' num2str(Le) ', M = 2, N = 3, $\varepsilon$ = 5'])
xlabel('u'); ylabel('[dB]')
legend('Product','CBF','location','southeast')

f10 = figure(10);
set(f10,'DefaultAxesFontSize',50)
clf
hold on
% quiver(-2/M,0,0,N/(M+N-1),'color',plotcolor.Marray)
% quiver(2/N,0,0,M/(M+N-1),'color',plotcolor.Narray)
% quiver(-2/N,0,0,-1/(M+N-1),'m')
% quiver(-2/M/N,0,0,1/(M+N-1),'color',plotcolor.CSACBF)
% 
% quiver(2/M,0,0,N/(M+N-1),'color',plotcolor.Marray)
% quiver(-2/N,0,0,M/(M+N-1),'color',plotcolor.Narray)
% 
% quiver(-2/N,0,0,(M-1)/(M+N-1),'color',plotcolor.CSACBF)
% quiver(2/N,0,0,(M-1)/(M+N-1),'color',plotcolor.CSACBF)
% quiver(2/N,0,0,-1/(M+N-1),'m')
% 
% quiver(2/M,0,0,(N-1)/(M+N-1),'color',plotcolor.CSACBF)
% quiver(2/M,0,0,-1/(M+N-1),'m')
% quiver(-2/M,0,0,(N-1)/(M+N-1),'color',plotcolor.CSACBF)
% quiver(-2/M,0,0,-1/(M+N-1),'m')
% 
% quiver(-2/M/N,0,0,-1/(M+N-1),'m')
% quiver(2/M/N,0,0,-1/(M+N-1),'m')
% quiver(2/M/N,0,0,1/(M+N-1),'color',plotcolor.CSACBF)

stem(-2/M,(N/(M+N-1)),'s','filled','MarkerSize',25,'color',plotcolor.Marray)
stem(2/N,(M/(M+N-1)),'^','filled','MarkerSize',25,'color',plotcolor.Narray)
stem(-2/N,-(1/(M+N-1)),'vm','filled','MarkerSize',25)
stem(-2/M/N,(1/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF)


stem(2/M,(N/(M+N-1)),'s','filled','MarkerSize',25,'color',plotcolor.Marray)
stem(-2/N,(M/(M+N-1)),'^','filled','MarkerSize',25,'color',plotcolor.Narray)

stem(-2/N,((M-1)/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF)
stem(2/N,((M-1)/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF)
stem(2/N,-(1/(M+N-1)),'vm','filled','MarkerSize',25)

stem(2/M,((N-1)/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF)
stem(2/M,-(1/(M+N-1)),'vm','filled','MarkerSize',25)
stem(-2/M,((N-1)/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF)
stem(-2/M,-(1/(M+N-1)),'vm','filled','MarkerSize',25)

stem(-2/M/N,-(1/(M+N-1)),'vm','filled','MarkerSize',25)
stem(2/M/N,-(1/(M+N-1)),'vm','filled','MarkerSize',25)
h = stem(2/M/N,(1/(M+N-1)),'+','MarkerSize',30,'color',plotcolor.CSACBF);
set(h,'ShowBaseLine','off');


% plot(u,abs(KMNe),'color',plotcolor.CSAPP)
line_fewer_markers(u,abs(KMNe),45,'kx','MarkerSize',30,'Spacing','x','LockOnMax',1);
% plot(u,abs(KLambda),'color',plotcolor.CSACBF)
line_fewer_markers(u,abs(KLambda),47,'+','MarkerSize',30,'color',plotcolor.CSACBF,'Spacing','x','LockOnMax',1);
xlim([-1.03 1.03])
ylim([-0.35 1.1])
% title('Finite Aperture CSA$_{\mathrm{cbf}}$ and CSA$_{\mathrm{pp}}$ Beampatterns with Subarray Impulses')
[legh,objh,outh,outm] = legend('$M$','$N$','$P$','$\left|K^{cbf}_{csa}(u)\right|$','location','northeast','Orientation','Horizontal');
xlabel('$u$'); ylabel('Linear')


%%

% f1 = @(x,y)line_fewer_markers(x,y,5,'ro');
% f2 = @(x,y)line_fewer_markers(x,y,8,'-.bs','Spacing','curve','markerfacecolor','g');%,'LegendLine','off');
% [yH,lh1,lh2] = plotyy(x,y1,x,y2,f1,f2)

f11 = figure(11);
set(f11,'DefaultAxesFontSize',50)
clf
hold on
% plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSAPP)
line_fewer_markers(u,10*log10(abs(KMNe)),150,'kx','MarkerSize',30,'Spacing','x','lockonmax',1);
% plot(u,10*log10(abs(KLambda)),'color',plotcolor.CSACBF)
line_fewer_markers(u,10*log10(abs(KLambda)),175,'+','MarkerSize',30,'color',plotcolor.CSACBF,'Spacing','x','lockonmax',1);
% plot(u,20*log10(abs(WULAe)),'color',plotcolor.ULA)
line_fewer_markers(u,20*log10(abs(WULAe)),150,'ro','MarkerSize',20,'MarkerFaceColor','r','Spacing','x','lockonmax',1);
ylim([-40 0])
% title(['ULA, CSA$_{\mathrm{pp}}$, CSA$_{\mathrm{cbf}}$ Beampatterns: L = ' num2str(Le) ', M = 2, N = 3, $\beta$ = 5'])
xlabel('$u$'); ylabel('[dB]')
legend('CSA$_{\mathrm{pp}}$','CSA$_{\mathrm{cbf}}$','ULA','location','northeast')