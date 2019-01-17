%%
plotcolor.ULA       = [1 0 0];          % ULA: red
plotcolor.M    = [0 .25 .9];       % M: blue
plotcolor.N    = [0 0.6 0.2];      % N: green
plotcolor.CSAcbf    = [1 .6 0];         % CSAcbf: orange
plotcolor.CSApp     = [0 0 0];          % CSApp: purple
plotcolor.sim       = [.25 .25 .25];    % sim: almost black
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Array Gain plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
present(0)
load SimResultsdenseAG_M2N3B25.mat;
load SimResultsdenseAG_M2N5B10.mat;
load SimResultsdenseAG_M5N6B10.mat;
alpha_AG = alpha_M2N3B25_denseAG;
ULA_AG_Derived = derived_M2N5B10_denseAG.ULA.AG;
ULA_AG_Simulated = simulated_M2N5B10_denseAG.ULA.AG;

CSAcbf_Aperture_AG_derived = derived_M2N5B10_denseAG.CSAcbf.AG;
CSAcbf_Aperture_AG_simulated = simulated_M2N5B10_denseAG.CSAcbf.AG;
CSAcbf_Sensor_AG_derived1 = derived_M5N6B10_denseAG.CSAcbf.AG;
CSAcbf_Sensor_AG_simulated1 = simulated_M5N6B10_denseAG.CSAcbf.AG;
CSAcbf_Sensor_AG_derived2 = derived_M2N3B25_denseAG.CSAcbf.AG;
CSAcbf_Sensor_AG_simulated2 = simulated_M2N3B25_denseAG.CSAcbf.AG;

CSApp_Aperture_AG_derived = derived_M2N5B10_denseAG.CSApp.AG;
CSApp_Aperture_AG_simulated = simulated_M2N5B10_denseAG.CSApp.AG;
CSApp_Sensor_AG_derived1 = derived_M5N6B10_denseAG.CSApp.AG;
CSApp_Sensor_AG_simulated1 = simulated_M5N6B10_denseAG.CSApp.AG;
CSApp_Sensor_AG_derived2 = derived_M2N3B25_denseAG.CSApp.AG;
CSApp_Sensor_AG_simulated2 = simulated_M2N3B25_denseAG.CSApp.AG;


f1 = figure(1);
set(f1,'DefaultAxesFontSize',24)
for b = 1:1
    clf
    semilogx(1-alpha_AG,10*log10(abs(ULA_AG_Derived)),'color',plotcolor.ULA)
    hold on % needs to be after first semilogx command
    semilogx(1-alpha_AG,10*log10(abs(CSAcbf_Sensor_AG_derived1)),'color',plotcolor.CSAcbf)
    semilogx(1-alpha_AG,10*log10(abs(CSApp_Sensor_AG_derived1)),'color',plotcolor.CSApp)
    semilogx(1-alpha_AG,10*log10(abs(ULA_AG_Simulated)),':^','MarkerSize',12,'Color',plotcolor.ULA)
    semilogx(1-alpha_AG,10*log10(abs(CSAcbf_Sensor_AG_simulated1)),':+','MarkerSize',15,'color',plotcolor.CSAcbf)
    semilogx(1-alpha_AG,10*log10(abs(CSApp_Sensor_AG_simulated1)),':x','MarkerSize',15,'color',plotcolor.CSApp)
%     semilogx(1-alpha_AG,10*log10(abs(CSAcbf_Aperture_AG_derived)),'color',plotcolor.CSAcbf)
%     semilogx(1-alpha_AG,10*log10(abs(CSAcbf_Aperture_AG_simulated)),':^','color',plotcolor.CSAcbf)
%     semilogx(1-alpha_AG,10*log10(abs(CSApp_Aperture_AG_derived)),'color',plotcolor.CSApp)
%     semilogx(1-alpha_AG,10*log10(abs(CSApp_Aperture_AG_simulated)),':^','color',plotcolor.CSApp)
%     
    grid on
    ax = gca;
    xtick = linspace(.1,1,10); set(ax,'XTick',xtick);
    xticklabel = abs(round(xtick.*100)./100 - 1);
%     xticklabel = ['.1';'.2';'.3';'.4';'.5';'.6';'.7';'.8';'.9'];
    set(ax,'Xdir','reverse','XTickLabel',xticklabel)
%     ax.XTickLabelRotation=45;
%     xlim([1-alpha(end) 1]);
    xlim([1-alpha_AG(end) 1]);
    ylim([0 inf]);
    xlabel('$\alpha$'); ylabel('[dB]')
    legend('','Derived','','','Monte-Carlo','','location','northeast')
%     legend('ULA','','M','','N','','CSA$_{\mathrm{cbf}}$','','CSA$_{\mathrm{pp}}$','','location','northeast')
%     set(legend,'FontSize',20)
%     set(findall(gcf,'-property','FontSize',18))
%     title(['Array Gain for ULA of L = '...
%         num2str(L) ', and CSA Parameters $\beta$ = ' num2str(beta(b)) ', M = ' num2str(M) ', N = ' num2str(N) ', $u_s$ = ' num2str(round(u0*100)/100)])
%     title({'Fixed Aperture Constraint:',['ULA of L = ' num2str(L) ', $\Lambda$ = ' num2str(Lambda) ', M = ' num2str(M) ', N = ' num2str(N) ', $\beta$ = ' num2str(beta)]})
%     title({'Array Gain vs. $\alpha$ - Fixed No. Sensors Constraint (100)',['ULA Aperture $\approx$ 50$\lambda_o$, CSA Aperture $\approx$ 150$\lambda_o$, M = 5, N = 6, $\beta$ = 10']})
end

%%
load SimResultsROC_M5N6B10.mat;
load SimResultsROC_M2N5B10.mat;

ULA_Pfa_Derived = derived_M2N5B10_wROCs.ULA.Pfa;
ULA_Pd_Derived = derived_M2N5B10_wROCs.ULA.Pd;
ULA_Pfa_Simulated = simulated_M2N5B10_wROCs.ULA.Pfa;
ULA_Pd_Simulated = simulated_M2N5B10_wROCs.ULA.Pd;

CSAcbf_Pfa_Derived = derived_M5N6B10_wROCs.CSAcbf.Pfa;
CSAcbf_Pd_Derived = derived_M5N6B10_wROCs.CSAcbf.Pd;
CSAcbf_Pfa_Simulated = simulated_M5N6B10_wROCs.CSAcbf.Pfa;
CSAcbf_Pd_Simulated = simulated_M5N6B10_wROCs.CSAcbf.Pd;

CSApp_Pfa_Derived = derived_M5N6B10_wROCs.CSApp.Pfa;
CSApp_Pd_Derived = derived_M5N6B10_wROCs.CSApp.Pd;
CSApp_Pfa_Simulated = simulated_M5N6B10_wROCs.CSApp.Pfa;
CSApp_Pd_Simulated = simulated_M5N6B10_wROCs.CSApp.Pd;

alpha_ROC = alpha_M2N5B10_wROCs;
% ULA_Pd_Derived = derived_M2N5B10_wROCs.ULA.Pd;
% ULA_Pd_Derived = derived_M2N5B10_wROCs.ULA.Pd;
f12 = figure(10);
set(f12,'DefaultAxesFontSize',24);
    for dummy = 1
        ROC1 = 4;
        
        ROC2 = 7;
        ROC3 = 10;
        
        skip = 2;
%         clf
        hold on
        
        plot(log10(ULA_Pfa_Simulated{1,ROC1}(1:skip:end)),log10(ULA_Pd_Simulated{1,ROC1}(1:skip:end)),':^','MarkerSize',12,'color',plotcolor.ULA)
        plot(log10(CSAcbf_Pfa_Simulated{1,ROC1}(1:skip:end)),log10(CSAcbf_Pd_Simulated{1,ROC1}(1:skip:end)),':+','MarkerSize',15,'color',plotcolor.CSAcbf)
        plot(log10(CSApp_Pfa_Simulated{1,ROC1}(1:skip:end)),log10(CSApp_Pd_Simulated{1,ROC1}(1:skip:end)),':x','MarkerSize',18,'color',plotcolor.CSApp)
        
        plot(log10(ULA_Pfa_Simulated{1,ROC2}(1:skip:end)),log10(ULA_Pd_Simulated{1,ROC2}(1:skip:end)),':s','MarkerSize',10,'color',plotcolor.ULA)
        plot(log10(CSAcbf_Pfa_Simulated{1,ROC2}(1:skip:end)),log10(CSAcbf_Pd_Simulated{1,ROC2}(1:skip:end)),':s','MarkerSize',10,'color',plotcolor.CSAcbf)
        plot(log10(CSApp_Pfa_Simulated{1,ROC2}(1:skip:end)),log10(CSApp_Pd_Simulated{1,ROC2}(1:skip:end)),':s','MarkerSize',10,'color',plotcolor.CSApp)
        
        plot(log10(ULA_Pfa_Simulated{1,ROC3}),log10(ULA_Pd_Simulated{1,ROC3}),':d','MarkerSize',10,'color',plotcolor.ULA)
        plot(log10(CSAcbf_Pfa_Simulated{1,ROC3}),log10(CSAcbf_Pd_Simulated{1,ROC3}),':d','MarkerSize',10,'color',plotcolor.CSAcbf)
        plot(log10(CSApp_Pfa_Simulated{1,ROC3}),log10(CSApp_Pd_Simulated{1,ROC3}),':d','MarkerSize',10,'color',plotcolor.CSApp)
% %         
        plot(log10(ULA_Pfa_Derived{1,ROC1}),log10(ULA_Pd_Derived{1,ROC1}),'color',plotcolor.ULA)
        plot(log10(CSAcbf_Pfa_Derived{1,ROC1}),log10(CSAcbf_Pd_Derived{1,ROC1}),'color',plotcolor.CSAcbf)
        plot(log10(CSApp_Pfa_Derived{1,ROC1}),log10(CSApp_Pd_Derived{1,ROC1}),'color',plotcolor.CSApp)
%         
        plot(log10(ULA_Pfa_Derived{1,ROC2}),log10(ULA_Pd_Derived{1,ROC2}),'color',plotcolor.ULA)
        plot(log10(CSApp_Pfa_Derived{1,ROC2}),log10(CSApp_Pd_Derived{1,ROC2}),'color',plotcolor.CSApp)
        plot(log10(CSAcbf_Pfa_Derived{1,ROC2}),log10(CSAcbf_Pd_Derived{1,ROC2}),'color',plotcolor.CSAcbf)
        
        plot(log10(ULA_Pfa_Derived{1,ROC3}),log10(ULA_Pd_Derived{1,ROC3}),'color',plotcolor.ULA)
        plot(log10(CSAcbf_Pfa_Derived{1,ROC3}),log10(CSAcbf_Pd_Derived{1,ROC3}),'color',plotcolor.CSAcbf)
        plot(log10(CSApp_Pfa_Derived{1,ROC3}),log10(CSApp_Pd_Derived{1,ROC3}),'color',plotcolor.CSApp)
        
        ylim([-.07 0])
        xlim([-3.5 0])
%         legend('','Monte-Carlo','','','Derived','','','location','southeast')
        legend('',['$\alpha = $ ' num2str(alpha_ROC(ROC1))],'','',['$\alpha = $ ' num2str(alpha_ROC(ROC2))],'','',['$\alpha = $ ' num2str(alpha_ROC(ROC3))],'','location','northwest')
        ylabel('log$_{10}(P_d)$')
        xlabel('log$_{10}(P_{fa})$')
%         title({['ROCs for  $\alpha = $ ' num2str(alpha_ROC(ROC1)) ' - Fixed No. Sensors Constraint (100)'],...
%             'ULA Aperture $\approx$ 50$\lambda_o$, CSA Aperture $\approx$ 150$\lambda_o$, M = 5, N = 6, $\beta$ = 10'})
%         title({['ROCs for  increasing $\alpha$ - Fixed No. Sensors Constraint (100)'],...
%             'ULA Aperture $\approx$ 50$\lambda_o$, CSA Aperture $\approx$ 150$\lambda_o$, M = 5, N = 6, $\beta$ = 10'})

%         axis tight
    end
    
    %%
ULA_xSigNoisePDF_Derived = derived_M2N5B10_wROCs.CSAcbf.xSigNoisePDF;
ULA_SigNoisePDF_Derived = derived_M2N5B10_wROCs.CSAcbf.SigNoisePDF;
ULA_xSigNoisePDF_Simulated = simulated_M2N5B10_wROCs.CSAcbf.xSigNoisePDF;
ULA_SigNoisePDF_Simulated = simulated_M2N5B10_wROCs.CSAcbf.SigNoisePDF;

CSAcbf_xSigNoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.xSigNoisePDF;
CSAcbf_SigNoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.SigNoisePDF;
CSAcbf_xSigNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.xSigNoisePDF;
CSAcbf_SigNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.SigNoisePDF;

CSApp_xSigNoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.xSigNoisePDF;
CSApp_SigNoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.SigNoisePDF;
CSApp_xSigNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.xSigNoisePDF;
CSApp_SigNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.SigNoisePDF;




ULA_xNoisePDF_Derived = derived_M2N5B10_wROCs.ULA.xNoisePDF;
ULA_NoisePDF_Derived = derived_M2N5B10_wROCs.ULA.NoisePDF;
ULA_xNoisePDF_Simulated = simulated_M2N5B10_wROCs.ULA.xNoisePDF;
ULA_NoisePDF_Simulated = simulated_M2N5B10_wROCs.ULA.NoisePDF;

CSAcbf_xNoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.xNoisePDF;
CSAcbf_NoisePDF_Derived = derived_M5N6B10_wROCs.CSAcbf.NoisePDF;
CSAcbf_xNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.xNoisePDF;
CSAcbf_NoisePDF_Simulated = simulated_M5N6B10_wROCs.CSAcbf.NoisePDF;

CSApp_xNoisePDF_Derived = derived_M5N6B10_wROCs.CSApp.xNoisePDF;
CSApp_NoisePDF_Derived = derived_M5N6B10_wROCs.CSApp.NoisePDF;
CSApp_xNoisePDF_Simulated = simulated_M5N6B10_wROCs.CSApp.xNoisePDF;
CSApp_NoisePDF_Simulated = simulated_M5N6B10_wROCs.CSApp.NoisePDF;

f5 = figure(5);
set(f5,'DefaultAxesFontSize',24)
clf
aplot = 1
skip2 = 2;
subplot(2,4,[1:2 5:6])
hold on
plot(CSAcbf_xSigNoisePDF_Derived{1,aplot},CSAcbf_SigNoisePDF_Derived{1,aplot}...
    ,'LineWidth',4,'color',plotcolor.CSAcbf);
plot(CSAcbf_xSigNoisePDF_Simulated{1,aplot}(1:end-80),CSAcbf_SigNoisePDF_Simulated{1,aplot}(1:end-80)...
    ,':d','color',plotcolor.CSAcbf);
plot(CSAcbf_xSigNoisePDF_Simulated{1,aplot}(end-80:skip2:end),CSAcbf_SigNoisePDF_Simulated{1,aplot}(end-80:skip2:end)...
    ,':d','color',plotcolor.CSAcbf);


skip2 = 5;
hold on
plot(ULA_xSigNoisePDF_Derived{1,aplot},ULA_SigNoisePDF_Derived{1,aplot}...
    ,'LineWidth',4,'color',plotcolor.ULA);
%     /trapz(derived.ULA.xSigNoisePDF{1,aplot},derived.ULA.SigNoisePDF{1,aplot})...
    
plot(ULA_xSigNoisePDF_Simulated{1,aplot}(1:end-80),ULA_SigNoisePDF_Simulated{1,aplot}(1:end-80)...
    ,':d','color',plotcolor.ULA);
plot(ULA_xSigNoisePDF_Simulated{1,aplot}(end-80:skip2:end),ULA_SigNoisePDF_Simulated{1,aplot}(end-80:skip2:end)...
    ,':d','color',plotcolor.ULA);
ylabel('Fixed No. Sensors')
xlabel('x')
% ULAxmax = derived.ULA.xSigNoisePDF{1,aplot}(end);
% ULAymax = round(max(derived.ULA.SigNoisePDF{1,aplot})*100)/100;
% axis([0 2 0 ULAymax+.2*ULAymax])

plot(CSApp_xSigNoisePDF_Derived{1,aplot},abs(CSApp_SigNoisePDF_Derived{1,aplot})...
    ,'LineWidth',4,'color',plotcolor.CSApp);
plot(CSApp_xSigNoisePDF_Simulated{1,aplot}(1:end-80),CSApp_SigNoisePDF_Simulated{1,aplot}(1:end-80)...
    ,':d','color',plotcolor.CSApp);
plot(CSApp_xSigNoisePDF_Simulated{1,aplot}(end-80:skip2:end),CSApp_SigNoisePDF_Simulated{1,aplot}(end-80:skip2:end)...
    ,':d','color',plotcolor.CSApp);




subplot(2,4,[3:4 7:8])
hold on

        hold on
%         bar(simulated.CSAcbf.xNoisePDF{1,aplot},simulated.CSAcbf.NoisePDF{1,aplot}...
%             ,'FaceColor',plotcolor.sim);
        plot(derived.CSAcbf.xNoisePDF{1,aplot},derived.CSAcbf.NoisePDF{1,aplot}...
            ,'LineWidth',4,'color',plotcolor.CSAcbf);
        plot(simulated.CSAcbf.xNoisePDF{1,aplot}(1:end-80),simulated.CSAcbf.NoisePDF{1,aplot}(1:end-80)...
            ,':d','color',plotcolor.CSAcbf);
        plot(simulated.CSAcbf.xNoisePDF{1,aplot}(end-80:skip2:end),simulated.CSAcbf.NoisePDF{1,aplot}(end-80:skip2:end)...
            ,':d','color',plotcolor.CSAcbf);

        
        skip1 = 1;
        skip2 = 5;
        hold on
        plot(derived.ULA.xNoisePDF{1,aplot},derived.ULA.NoisePDF{1,aplot}...
            /trapz(derived.ULA.xNoisePDF{1,aplot},derived.ULA.NoisePDF{1,aplot})...
            ,'LineWidth',4,'color',plotcolor.ULA);
        plot(simulated.ULA.xNoisePDF{1,aplot}(1:end-80),simulated.ULA.NoisePDF{1,aplot}(1:end-80)...
            ,':d','color',plotcolor.ULA);
        plot(simulated.ULA.xNoisePDF{1,aplot}(end-80:skip2:end),simulated.ULA.NoisePDF{1,aplot}(end-80:skip2:end)...
            ,':d','color',plotcolor.ULA);
        xlabel(['ULA: $\chi^2$-Distributed'])
        ULAxmax = derived.ULA.xNoisePDF{1,aplot}(end);
        ULAymax = round(max(derived.ULA.NoisePDF{1,aplot})*100)/100;
        axis([0 ULAxmax 0 ULAymax+.2*ULAymax])

        
        hold on
%         bar(simulated.CSApp.xNoisePDF{1,aplot},simulated.CSApp.NoisePDF{1,aplot}...
%             ,'FaceColor',plotcolor.sim);
        plot(derived.CSApp.xNoisePDF{1,aplot},abs(derived.CSApp.NoisePDF{1,aplot})...
            ,'LineWidth',4,'color',plotcolor.CSApp);
        plot(simulated.CSApp.xNoisePDF{1,aplot}(1:end-80),simulated.CSApp.NoisePDF{1,aplot}(1:end-80)...
            ,':d','color',plotcolor.CSApp);
        plot(simulated.CSApp.xNoisePDF{1,aplot}(end-80:skip2:end),simulated.CSApp.NoisePDF{1,aplot}(end-80:skip2:end)...
            ,':d','color',plotcolor.CSApp);
%         xlabel(['CSA$_{\mathrm{pp}}$: $\mathrm{R}_M \cdot \mathrm{R}_N$ Distributed'])
%         CSAprodxmax = derived.CSApp.xNoiseCDF{1,aplot}(end);
%         CSAprodymax = round(max(derived.CSApp.NoisePDF{1,aplot})*100)/100;
%         axis([0 CSAprodxmax 0 CSAprodymax+.2*CSAprodymax])
NoisePDFanimation = getframe(gcf);
if doMovie
    writeVideo(NoisePDFmovie,NoisePDFanimation);
end