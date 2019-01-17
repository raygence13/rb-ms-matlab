%% Script for MS Thesis and Journal Figures

clear
close all
clc
presentJournal(0)
%%   
% Color code for each undersampling factor
plotcolor.ULA       = [1 0 0];          % ULA: red
plotcolor.Marray    = [0 .25 .9];       % Marray: blue
plotcolor.Narray    = [0 .75 0];      % Narray: green
% plotcolor.CSAcbf    = [1 .6 0];         % CSACBF: orange
plotcolor.CSAcbf    = [0 0.5 0.5];         % CSAcbf: orange
plotcolor.CSApp     = [0 0 0];          % CSAPP: purple
plotcolor.sim       = [.25 .25 .25];    % sim: almost black

uNfft = 2048;
ur = linspace(0,1,uNfft/2);
ul = linspace(-1,0-1/uNfft,uNfft/2);
u = [ul ur];

%% Chapter 2.1 ULA
% ULA Beampattern
L = 12;
v = ones(L,1);
wULA = v/L;
WULA = fftshift(fft(wULA,uNfft));

f1 = figure(1);
% clf
plot(u,20*log10(abs(WULA)),'color',plotcolor.ULA)
ylim([-35 0])
title(['ULA Beampattern, L = ' num2str(L)])
xlabel('u = sin($\theta$)'); ylabel('[dB]')

%% ULA Coarray
Iula = v';
kappaULA = conv(Iula,fliplr(Iula));
gamma = 1-L:L-1;

f2 = figure(2);
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
xlabel('$u$')

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
plot(u,10*log10(abs(KMN)),'color',plotcolor.CSApp)
ylim([-25 0])
xlabel('$u$'); ylabel('[dB]')
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
plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSApp)
ylim([-35 0])
xlabel('$u$'); ylabel('[dB]')
legend('ULA','CSA$_{\mathrm{pp}}$','location','southeast')

f7 = figure(7);
clf
stem(gamma,kappaMNe,'x','color',plotcolor.CSApp)
title('CSA$_{\mathrm{pp}}$ coarray, $\kappa^{pp}_{csa}[\gamma]$')
ylim([0 5.01])
xlabel('$\gamma$'); ylabel('Redundancies')

%% CSA_CBF
kappaLambda = conv(I_Lambdae,fliplr(I_Lambdae));
KLambda = fftshift(fft(kappaLambda/Lambda^2,uNfft));

presentJournal(0)
f8 = figure(8);
clf
% plot(-11:11,kappaULA,'ro','MarkerSize',25,'MarkerFaceColor',plotcolor.ULA)
stem(-11:11,kappaULA,'ro','MarkerSize',25,'MarkerFaceColor',plotcolor.ULA)
hold on
stem(-17:17,kappaLambda,'+','MarkerSize',30,'color',plotcolor.CSAcbf)
stem(-17:17,kappaMNe,'x','MarkerSize',30,'color',plotcolor.CSApp)
% title('CSA$_{\mathrm{cbf}}$ coarray, $\kappa^{cbf}_{csa}[\gamma]$')
ylim([0 Lambda+1])
xlabel('$\gamma$'); ylabel('Redundancies')
legend('ULA','CSA$_{\textrm{cbf}}$','CSA$_{\textrm{pp}}$')

f9 = figure(9);
clf
hold on
plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSApp)
plot(u,10*log10(abs(KLambda)),'color',plotcolor.CSAcbf)
ylim([-35 0])
title(['CSA$_{\mathrm{pp}}$ vs. CSA$_{\mathrm{cbf}}$ Beampatterns: L = ' num2str(Le) ', $M = 2$, $N = 3$, $P = 5$'])
xlabel('$u$'); ylabel('[dB]')
legend('Product','CBF','location','southeast')

f10 = figure(10);
set(f10,'DefaultAxesFontSize',50)
hold on
Mplot = stem(-2/M,(N/(M+N-1)),'s','filled','color',plotcolor.Marray);
Nplot = stem(2/N,(M/(M+N-1)),'^','filled','color',plotcolor.Narray);
Pplot = stem(-2/N,-(1/(M+N-1)),'vm','filled');
Kplot = stem(-2/M/N,(1/(M+N-1)),'+','color',plotcolor.CSAcbf);


stem(2/M,(N/(M+N-1)),'s','filled','color',plotcolor.Marray)
stem(-2/N,(M/(M+N-1)),'^','filled','color',plotcolor.Narray)

stem(-2/N,((M-1)/(M+N-1)),'+','color',plotcolor.CSAcbf)
stem(2/N,((M-1)/(M+N-1)),'+','color',plotcolor.CSAcbf)
stem(2/N,-(1/(M+N-1)),'vm','filled')

stem(2/M,((N-1)/(M+N-1)),'+','color',plotcolor.CSAcbf)
stem(2/M,-(1/(M+N-1)),'vm','filled')
stem(-2/M,((N-1)/(M+N-1)),'+','color',plotcolor.CSAcbf)
stem(-2/M,-(1/(M+N-1)),'vm','filled')

stem(-2/M/N,-(1/(M+N-1)),'vm','filled')
stem(2/M/N,-(1/(M+N-1)),'vm','filled')
h = stem(2/M/N,(1/(M+N-1)),'+','color',plotcolor.CSAcbf);
set(h,'ShowBaseLine','off');


% plot(u,abs(KMNe),'color',plotcolor.CSApp)
Kppplot = line_fewer_markers(u,abs(KMNe),45,'kx','MarkerSize',30,'Spacing','x','LockOnMax',1);
% plot(u,abs(KLambda),'color',plotcolor.CSAcbf)
line_fewer_markers(u,abs(KLambda),47,'+','MarkerSize',30,'color',plotcolor.CSAcbf,'Spacing','x','LockOnMax',1);
xlim([-1.03 1.03])
ylim([-0.35 1.1])
% title('Finite Aperture CSA$_{\mathrm{cbf}}$ and CSA$_{\mathrm{pp}}$ Beampatterns with Subarray Impulses')
[legh,objh,outh,outm] = legend({'$M$','$N$','$P$','$\left|K^{cbf}_{csa}(u)\right|$'},'location','northeast','Orientation','Horizontal','FontSize',50);
xlabel('$u$'); ylabel('Linear')


%%

% f1 = @(x,y)line_fewer_markers(x,y,5,'ro');
% f2 = @(x,y)line_fewer_markers(x,y,8,'-.bs','Spacing','curve','markerfacecolor','g');%,'LegendLine','off');
% [yH,lh1,lh2] = plotyy(x,y1,x,y2,f1,f2)

f11 = figure(11);
clf
hold on
nmarkers = 150;
% plot(u,10*log10(abs(KMNe)),'color',plotcolor.CSApp)
line_fewer_markers(u,10*log10(abs(KMNe)),nmarkers,'kx','MarkerSize',32,'Spacing','x','lockonmax',1);
% plot(u,10*log10(abs(KLambda)),'color',plotcolor.CSAcbf)
line_fewer_markers(u,10*log10(abs(KLambda)),nmarkers,'+','MarkerSize',32,'color',plotcolor.CSAcbf,'Spacing','x','lockonmax',1);
% plot(u,20*log10(abs(WULAe)),'color',plotcolor.ULA)
line_fewer_markers(u,20*log10(abs(WULAe)),nmarkers,'r:o','MarkerSize',18,'MarkerFaceColor','r','Spacing','x','lockonmax',1);
ylim([-40 0])
% title(['ULA, CSA$_{\mathrm{pp}}$, CSA$_{\mathrm{cbf}}$ Beampatterns: L = ' num2str(Le) ', M = 2, N = 3, $\beta$ = 5'])
xlabel('$u$'); ylabel('[dB]')
legend('CSA$_{\mathrm{pp}}$','CSA$_{\mathrm{cbf}}$','ULA','location','northeast')


%% AG Plots with simulated data


% load Jan302018_CSA_SimResults_M2N5P10_SNRin-0dB_uS-0_lalpha-95_nSamples15000.mat;
% tempD = derived.ULA;
% tempS = simulated.ULA;
% % simulated.ULA = simulated_M2N5B10_denseAG.ULA;
% load Jan302018_CSA_SimResults_M5N6P10_SNRin-0dB_uS-0_lalpha-95_nSamples15000.mat;
% derived.ULA = tempD;
% simulated.ULA = tempS;
% % derived.CSAcbf = derived_M5N6B10_denseAG.CSAcbf;
% % simulated.CSAcbf = simulated_M5N6B10_denseAG.CSAcbf;
% % 
% % derived.CSApp = derived_M5N6B10_denseAG.CSApp;
% % simulated.CSApp = simulated_M5N6B10_denseAG.CSApp;
% % alpha = alpha_M5N6B10_denseAG;
% 
% f12 = figure(12);
% for b = 1:1
%     clf
%     semilogx(1-alpha,10*log10(abs(derived.ULA.AG(b,:))),'color',plotcolor.ULA)
%     hold on % needs to be after first semilogx command
%     semilogx(1-alpha,10*log10(abs(derived.CSAcbf.AG(b,:))),'color',plotcolor.CSAcbf)
%     semilogx(1-alpha,10*log10(abs(derived.CSApp.AG(b,:))),'color',plotcolor.CSApp)
%     semilogx(1-alpha(1:2:end),10*log10(abs(simulated.ULA.AG(b,1:2:end))),':o','MarkerSize',20,'MarkerFaceColor',plotcolor.ULA,'Color',plotcolor.ULA)
%     semilogx(1-alpha(1:2:end),10*log10(abs(simulated.CSAcbf.AG(b,1:2:end))),':+','MarkerSize',30,'color',plotcolor.CSAcbf)
%     semilogx(1-alpha(1:2:end),10*log10(abs(simulated.CSApp.AG(b,1:2:end))),':x','MarkerSize',30,'color',plotcolor.CSApp)
%     grid on
%     ax = gca;
%     xtick = linspace(.1,1,10); set(ax,'XTick',xtick);
%     xticklabel = abs(round(xtick.*100)./100 - 1);
%     set(ax,'Xdir','reverse','XTickLabel',xticklabel)
% %     ax.XTickLabelRotation=45;
% %     xlim([1-alpha(end) 1]);
%     xlim([1-alpha(end) 1]);
%     ylim([0 inf]);
%     xlabel('$\alpha$'); ylabel('[dB]')
%     legend('ULA','','CSA$_{\mathrm{cbf}}$','','CSA$_{\mathrm{pp}}$','','location','northeast')
% end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Noise Only Biased PSD estimate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Jan302018_CSA_SimResults_M2N5P10_SNRin-0dB_uS-0_lalpha-10_nSamples15000.mat;
tempD = derived.ULA;
tempS = simulated.ULA;

load Jan302018_CSA_SimResults_M5N6P10_SNRin-0dB_uS-0_lalpha-10_nSamples15000.mat;
derived.ULA = tempD;
simulated.ULA = tempS;
u0 = 0;
% load SimResultsdenseAG_M2N5B10.mat;
% derived.CSAcbf = derived_M2N5B10_denseAG.CSAcbf;
% simulated.CSAcbf = simulated_M2N5B10_denseAG.CSAcbf;
% 
% derived.CSApp = derived_M2N5B10_denseAG.CSApp;
% simulated.CSApp = simulated_M2N5B10_denseAG.CSApp;
% alpha = alpha_M2N5B10_denseAG;
presentJournal(0)
f2 = figure(2);
clf
set(f2,'DefaultAxesFontSize',50)
plotcolor.truePSD = [0.6 0.6 0.6];
aplot = 8;
clf
hold on
%     box on
%     for plotind = 1:length(zeta)
plot(u,10*log10(abs(derived.truePSDNoise{1,aplot}))','LineWidth',10,'color',plotcolor.truePSD)
xlim([-0.4 0.4])
% ylim([-0.5 0.5])
xlabel('$u$'); ylabel('[dB]')

f3 = figure(3);
clf
set(f3,'DefaultAxesFontSize',50)
line_fewer_markers(u,10*log10(abs(KMNe)),153,'kx','MarkerSize',30,'color',plotcolor.CSApp)
hold on
line_fewer_markers(u,10*log10(abs(KLambda)),150,'+','MarkerSize',20,'color',plotcolor.CSAcbf)
line_fewer_markers(u,20*log10(abs(WULA)),149,'o','MarkerFaceColor',plotcolor.ULA,'MarkerSize',15,'color',plotcolor.ULA)
xlabel('$u$'); ylabel('[dB]')
legend('CSA$_{\textrm{pp}}$','CSA$_{\textrm{cbf}}$','ULA')
ylim([-50 0])
xlim([-0.4 0.4])
% plot(u,10*log10(abs(derived.ULA.NoisePSD{1,aplot})),'color',plotcolor.ULA)
% line_fewer_markers(u,10*log10(abs(simulated.ULA.NoisePSD{1,aplot})),25,':o','MarkerFaceColor','r','MarkerSize',20,'color',plotcolor.ULA)
% 
% plot(u,10*log10(abs(derived.CSApp.NoisePSD{1,aplot})),'color',plotcolor.CSApp)
% line_fewer_markers(u,10*log10(abs(simulated.CSApp.NoisePSD{1,aplot})),25,'k:x','MarkerSize',30,'color',plotcolor.CSApp)
% 
% plot(u,10*log10(abs(derived.CSAcbf.NoisePSD{1,aplot})),'color',plotcolor.CSAcbf)
% line_fewer_markers(u,10*log10(abs(simulated.CSAcbf.NoisePSD{1,aplot})),23,':+','MarkerSize',30,'color',plotcolor.CSAcbf)
% legend({'true PSD','ULA','','CSA$_{\mathrm{pp}}$','','CSA$_{\mathrm{cbf}}$',''},'location','northeast')

% xlabel('$u$'); ylabel('[dB]')
% ylim([-6 6])

%%
    
f3 = figure(3);
set(f3,'DefaultAxesFontSize',50)
plotcolor.truePSD = [1 0.6 0.78];
aplot = 3;
clf
hold on
%     box on
%     for plotind = 1:length(zeta)
sigline = plot([u0 u0],[-50 50].','-.m','LineWidth',6);
ULAderplot = plot(u,10*log10(abs(derived.ULA.SigNoisePSD{1,aplot})),'LineWidth',30,'color',plotcolor.truePSD);

CSAppderplot = plot(u,10*log10(abs(derived.CSApp.SigNoisePSD{1,aplot})),'color',plotcolor.CSApp);
CSAppsimplot = line_fewer_markers(u,10*log10(abs(simulated.CSApp.SigNoisePSD{1,aplot})),25,'k:x','MarkerSize',30,'color',plotcolor.CSApp);

CSAcbfderplot = plot(u,10*log10(abs(derived.CSAcbf.SigNoisePSD{1,aplot})),'color',plotcolor.CSAcbf);
CSAcbfsimplot = line_fewer_markers(u,10*log10(abs(simulated.CSAcbf.SigNoisePSD{1,aplot})),23,':+','MarkerSize',30,'color',plotcolor.CSAcbf);
legend([sigline ULAderplot CSAppderplot CSAppsimplot CSAcbfderplot CSAcbfsimplot],...
    {'Signal','ULA','CSA$_{\mathrm{pp}}$','','CSA$_{\mathrm{cbf}}$',''},'location','northeast')
xlabel('$u$'); ylabel('[dB]')
ylim([-10 25])
% %% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%    rho vs alpha plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load  rhoanalysis.mat
% f13 = figure(13);
% clf
% semilogx(1-alpha,rhoNoiseder,'r','LineWidth',4)
% hold on
% semilogx(1-alpha(1:2:end),abs(rhoNoisesim(1:2:end)),':dr','MarkerSize',9)
% %     semilogx(1-alpha,(rhoSigNoiseder(b,:)),'LineWidth',4,'color',[0 0 1])
% %     semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim(1:2:end)),':d','MarkerSize',9,'color',[0 0 1])
% semilogx(1-alpha,(rhoSigNoiseder_20dB),'m','LineWidth',4)%,'color',[1 .6 0])
% semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_20dB(1:2:end)),'m:d','MarkerSize',9)%,'color',[1 .6 0])
% semilogx(1-alpha,(rhoSigNoiseder_10dB),'b','LineWidth',4)
% semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_10dB(1:2:end)),'b:d','MarkerSize',9)
% semilogx(1-alpha,(rhoSigNoiseder_0dB),'LineWidth',4,'color',[0 .6 .1])
% semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim_0dB(1:2:end)),':d','MarkerSize',9,'color',[0 .6 .1])
% grid on
% ax = gca;
% xtick = linspace(.1,1,10);
% set(ax,'XTick',xtick);
% xticklabel = abs(round(xtick.*100)./100 - 1);
% set(ax,'XTickLabel',xticklabel);
% set(ax,'Xdir','reverse')
% %     ax.XTickLabelRotation=45;
% xlim([0.01 1]); ylim([0 1.1]);
% xlabel('$\alpha$')
% ylabel('$\rho$')
% legend('$\rho_{0}$','','$\rho_{1}$','SNR$_{in} = -20$dB','$\rho_{1}$','SNR$_{in} = -10$dB','$\rho_{1}$','SNR$_{in} = 0$dB','Location','SouthEast')
% title(['Correlation Coefficient $\rho$ vs. $\alpha$, $M$ = ' num2str(M) ', $N$ = ' num2str(N)])
%     

%% ROC plots
% load SimResultsROC_M2N5B10.mat
% derived.ULA = derived_M2N5B10_wROCs.ULA;
% simulated.ULA = simulated_M2N5B10_wROCs.ULA;
% 
% load SimResultsROC_M5N6B10.mat
% derived.CSAcbf = derived_M5N6B10_wROCs.CSAcbf;
% simulated.CSAcbf = simulated_M5N6B10_wROCs.CSAcbf;
% derived.CSApp = derived_M5N6B10_wROCs.CSApp;
% simulated.CSApp = simulated_M5N6B10_wROCs.CSApp;
% alpha = alpha_M5N6B10_wROCs;
f14 = figure(14);
presentJournal(0)
        ROC1 = 4;
        ROC2 = 4;
        ROC3 = 8;
        ROC4 = 10;
        skip1 = 3;
        endskip = 8;
        clf
%         plot(log10(derived.ULA.Pfa{1,ROC1}(1:5:end)),log10(derived.ULA.Pd{1,ROC1}(1:5:end)),':o','MarkerSize',30,'color',plotcolor.ULA)
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         plot(log10(derived.CSAcbf.Pfa{1,ROC1}(1:30:end)),log10(derived.CSAcbf.Pd{1,ROC1}(1:30:end)),':o','MarkerSize',30,'color',plotcolor.CSAcbf)
%         plot(log10(derived.CSApp.Pfa{1,ROC1}(1:40:end)),log10(derived.CSApp.Pd{1,ROC1}(1:40:end)),':o','MarkerSize',30,'color',plotcolor.CSApp)
% %         
%         loglog((derived.ULA.Pfa{1,ROC1}),(derived.ULA.Pd{1,ROC1}),'color',plotcolor.ULA)
%         hold on
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         loglog((derived.CSAcbf.Pfa{1,ROC1}),(derived.CSAcbf.Pd{1,ROC1}),'color',plotcolor.CSAcbf)
%         loglog((derived.CSApp.Pfa{1,ROC1}),(derived.CSApp.Pd{1,ROC1}),'color',plotcolor.CSApp)
%         
%         loglog((simulated.ULA.Pfa{1,ROC1}(1:5:end)),(simulated.ULA.Pd{1,ROC1}(1:5:end)),':^','MarkerSize',20,'MarkerFaceColor',plotcolor.ULA,'color',plotcolor.ULA)
%         hold on
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         loglog((simulated.CSAcbf.Pfa{1,ROC1}(1:5:end)),(simulated.CSAcbf.Pd{1,ROC1}(1:5:end)),':+','MarkerSize',30,'color',plotcolor.CSAcbf)
%         loglog((simulated.CSApp.Pfa{1,ROC1}(1:5:end)),(simulated.CSApp.Pd{1,ROC1}(1:5:end)),':x','MarkerSize',30,'color',plotcolor.CSApp)
% %         
%         
%         semilogx((derived.ULA.Pfa{1,ROC1}),(derived.ULA.Pd{1,ROC1}),'color',plotcolor.ULA)
%         hold on
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         semilogx((derived.CSAcbf.Pfa{1,ROC1}),(derived.CSAcbf.Pd{1,ROC1}),'color',plotcolor.CSAcbf)
%         semilogx((derived.CSApp.Pfa{1,ROC1}),(derived.CSApp.Pd{1,ROC1}),'color',plotcolor.CSApp)
%         
%         semilogx((simulated.ULA.Pfa{1,ROC1}(1:5:end)),(simulated.ULA.Pd{1,ROC1}(1:5:end)),':o','MarkerSize',20,'MarkerFaceColor',plotcolor.ULA,'color',plotcolor.ULA)
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         semilogx((simulated.CSAcbf.Pfa{1,ROC1}(1:5:end)),(simulated.CSAcbf.Pd{1,ROC1}(1:5:end)),':+','MarkerSize',30,'color',plotcolor.CSAcbf)
%         semilogx((simulated.CSApp.Pfa{1,ROC1}(1:5:end)),(simulated.CSApp.Pd{1,ROC1}(1:5:end)),':x','MarkerSize',30,'color',plotcolor.CSApp)
%         
%         
        semilogx((derived.ULA.Pfa{1,ROC2}),(derived.ULA.Pd{1,ROC2}),'color',plotcolor.ULA)
        hold on
%         plot(log10(ULAPfaDerived{1,ROC2}(1:5:end)),log10(ULAPdDerived{1,ROC2}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((derived.CSAcbf.Pfa{1,ROC2}),(derived.CSAcbf.Pd{1,ROC2}),'color',plotcolor.CSAcbf)
        semilogx((derived.CSApp.Pfa{1,ROC2}),(derived.CSApp.Pd{1,ROC2}),'color',plotcolor.CSApp)
        
        semilogx((simulated.ULA.Pfa{1,ROC2}(1:5:end)),(simulated.ULA.Pd{1,ROC2}(1:5:end)),':p','MarkerSize',20,'MarkerFaceColor',plotcolor.ULA,'color',plotcolor.ULA)
%         plot(log10(ULAPfaDerived{1,ROC2}(1:5:end)),log10(ULAPdDerived{1,ROC2}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((simulated.CSAcbf.Pfa{1,ROC2}(1:5:end)),(simulated.CSAcbf.Pd{1,ROC2}(1:5:end)),':p','MarkerSize',20,'color',plotcolor.CSAcbf)
        semilogx((simulated.CSApp.Pfa{1,ROC2}(1:5:end)),(simulated.CSApp.Pd{1,ROC2}(1:5:end)),':p','MarkerSize',20,'color',plotcolor.CSApp)
        
        
        semilogx((derived.ULA.Pfa{1,ROC3}),(derived.ULA.Pd{1,ROC3}),'color',plotcolor.ULA)
        hold on
%         plot(log10(ULAPfaDerived{1,ROC3}(1:5:end)),log10(ULAPdDerived{1,ROC3}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((derived.CSAcbf.Pfa{1,ROC3}),(derived.CSAcbf.Pd{1,ROC3}),'color',plotcolor.CSAcbf)
        semilogx((derived.CSApp.Pfa{1,ROC3}),(derived.CSApp.Pd{1,ROC3}),'color',plotcolor.CSApp)
        
        semilogx((simulated.ULA.Pfa{1,ROC3}(1:5:end)),(simulated.ULA.Pd{1,ROC3}(1:5:end)),':*','MarkerSize',30,'color',plotcolor.ULA)
%         plot(log10(ULAPfaDerived{1,ROC3}(1:5:end)),log10(ULAPdDerived{1,ROC3}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((simulated.CSAcbf.Pfa{1,ROC3}(1:5:end)),(simulated.CSAcbf.Pd{1,ROC3}(1:5:end)),':*','MarkerSize',30,'color',plotcolor.CSAcbf)
        semilogx((simulated.CSApp.Pfa{1,ROC3}(1:5:end)),(simulated.CSApp.Pd{1,ROC3}(1:5:end)),':*','MarkerSize',30,'color',plotcolor.CSApp)
        
        
        semilogx((derived.ULA.Pfa{1,ROC4}),(derived.ULA.Pd{1,ROC4}),'color',plotcolor.ULA)
        hold on
%         plot(log10(ULAPfaDerived{1,ROC4}(1:5:end)),log10(ULAPdDerived{1,ROC4}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((derived.CSAcbf.Pfa{1,ROC4}),(derived.CSAcbf.Pd{1,ROC4}),'color',plotcolor.CSAcbf)
        semilogx((derived.CSApp.Pfa{1,ROC4}),(derived.CSApp.Pd{1,ROC4}),'color',plotcolor.CSApp)
        
        semilogx((simulated.ULA.Pfa{1,ROC4}(1:2:end)),(simulated.ULA.Pd{1,ROC4}(1:2:end)),':d','MarkerSize',20,'color',plotcolor.ULA)
%         plot(log10(ULAPfaDerived{1,ROC4}(1:5:end)),log10(ULAPdDerived{1,ROC4}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
        semilogx((simulated.CSAcbf.Pfa{1,ROC4}(1:2:end)),(simulated.CSAcbf.Pd{1,ROC4}(1:2:end)),':d','MarkerSize',20,'color',plotcolor.CSAcbf)
        semilogx((simulated.CSApp.Pfa{1,ROC4}(1:2:end)),(simulated.CSApp.Pd{1,ROC4}(1:2:end)),':d','MarkerSize',20,'color',plotcolor.CSApp)
% 
% 
% %         plot(log10(derived.ULA.Pfa{1,ROC2}(1:8:end)),log10(derived.ULA.Pd{1,ROC2}(1:8:end)),':d','MarkerSize',10,'color',plotcolor.ULA)
% %         plot(log10(derived.CSAcbf.Pfa{1,ROC2}(1:15:end)),log10(derived.CSAcbf.Pd{1,ROC2}(1:15:end)),':d','MarkerSize',10,'color',plotcolor.CSAcbf)
% %         plot(log10(derived.CSApp.Pfa{1,ROC2}(1:60:end)),log10(derived.CSApp.Pd{1,ROC2}(1:60:end)),':d','MarkerSize',10,'color',plotcolor.CSApp)
% %         
% %         plot(log10(derived.ULA.Pfa{1,ROC3}(1:20:end)),log10(derived.ULA.Pd{1,ROC3}(1:20:end)),':s','MarkerSize',10,'color',plotcolor.ULA)
%         plot(log10(derived.CSAcbf.Pfa{1,ROC3}(1:50:end)),log10(derived.CSAcbf.Pd{1,ROC3}(1:50:end)),':s','MarkerSize',10,'color',plotcolor.CSAcbf)
%         plot(log10(derived.CSApp.Pfa{1,ROC3}(1:100:end)),log10(derived.CSApp.Pd{1,ROC3}(1:100:end)),':s','MarkerSize',10,'color',plotcolor.CSApp)
        
%         legend('',['$\alpha$ = ' num2str(alpha(ROC1))],'','',['$\alpha$ = ' num2str(alpha(ROC2))],'','',['$\alpha$ = ' num2str(alpha(ROC3))],'','location','southeast')
        ylabel('$\mathcal{P}_1$')
        xlabel('$\mathcal{P}_0$')
%         title({'ROCs: ULA (red), CSA$_{\mathrm{cbf}}$ (orange), and CSA$_{\mathrm{pp}}$ (black)',['SNR$_{in}=$ ' num2str(10*log10(varS/varW)) 'dB, L = ' num2str(L) ', $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ]})
%         xlim([log10(min(simulated.ULA.Pfa{1,ROC1}))-.1 0])
%         title({['SNR$_{in}=$ ' num2str(10*log10(varS/varW)) 'dB ROCs: ULA (red), CSA$_{\mathrm{cbf}}$ (orange), and CSA$_{\mathrm{pp}}$ (black)'],['ULA Aperture = 50$\lambda_o$, CSA Aperture = 150$\lambda_o$, $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ]})

%         ylim([log10(min(simulated.CSAcbf.Pd{1,ROC1}))-.01 0])
        ylim([.6 1])
        xlim([10^(-4) 1])
%         axis tight