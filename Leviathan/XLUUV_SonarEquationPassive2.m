% XLUUV_SonarEquationPassive.m
% Sonar equation analysis for XLUUV
% This function computes Passive Signal Excess
% SE == Signal Excess = SL - TL - NL + DI - DT
% TL == Transmission Loss (defined as positive), function of range  and
%       frequency; uses Thorp model for the attenuation coeficient
% NL == Noise Level, function of frequency
% DI == Directivity Index, function of size of receiver, steering direction
%       and Frequency; currently assume uniformly weighted planar array
% DT == Detection Threshold, a function of Probabilities of False Alarm and
%       Detection and Time-Bandwidth product

% Fyzodeen Khan, NUWCNPT Code 1513; 11 Sept 2017
% Ray Bautista, NUWCNPT Code 1511; 11 Oct 2017

present(0)
close all
clear all

islocal = 1;
%-------------------SETUP PARAMETERS---------------------------------------
% Source Level model
fc = 1000;     % corner frequency for SL change
SL1k = 140;    % SL below fc Hz (fixed from f(1) : 1000)

% Frequency band
f0 = 10000;     % maximum frequency (receiver design), Hz
f = 20:f0;      % frequency range, Hz

% Array spacing and DI
fd = 6000;      % desing frequency of receiver, Hz
c = 4900;       % sound speed, feet per sec

N = 10;         % Number of elements is one plane of panel receiver
M = 5;          % Number of elements is one plane of panel receiver
phi_s = 0;      % steering direction for DI calculation; 0 deg == broadside

% Range interval for transmission Loss computation
R = [10000 20000 30000 40000]; % Range, yards

% Broadband processor parameters
%W = 1000;       % passive processing bandwidth, Hz
%T = 1;          % passive integration time, sec
TW = 1000;      % Time-Bandwidth product
NumLines = 1;   % number of lines integrated
SPLoss = 0;     % expected loss in signal processing, dB
s = 1; s2 = s*s;	% H0 & H1 noise power
mu0 = 0;            % noise mean under H0
Pfa = 1e-4;         % False alarm probability
Pd1 = 0.5;          % Detection probability
det_indx = [5 : 1 : 45];	% Specify different detection index

% Ambient noise background parameters
SS = 3;         % sea-state level
ShipLevel = 5;  % shipping level
RainLevel = 0;  % Rain level
OT_flag = 'n';  % plot flag: 'y'=yes, 'n' = no
TL_flag = 's';  % Transmission Loss flag: 's'=spherical, 'c'=cylindrical
plot_flag = 'n';% plot flag: 'y'=yes, 'n' = no

% Additional levels for ambient and Source level
NL0 = 0;    % additional noise level => background may be louder
SL0 = 0;	% additional source level 

%-------------COMPUTE SONAR EQUATION PARAMETERS----------------------------

% Generate model for vessel Source Level
idx1k = find(f == fc);
SL(1: idx1k-1) = SL1k;  % put corner frequency at fc
SL(idx1k : length(f)) = SL1k - 20*log10(f(idx1k : length(f))/fc);
SL = SL0 + SL'; % This is -20 dB per decade starting at f0_sl

%%
% Compute the Detection Threshold
if islocal
    DT = -9.295764367342025; % if no norminv function
    % DT is obtained for following parameters
    % Need to use full MATLAB to calculate new DTs
%     m = sqrt(det_indx) * s;         % mean level from detection index equation
%     tau = norminv(1-Pfa, mu0, s);   % detection threshold under H0 fro Pfa
%     Pd = 1 - normcdf(tau, m, s);    % Pd for varying mean levels
%     mu1 = interp1(Pd, m, Pd1);      % find noise mean under H1 for desired Pd1
%     det_indx2 = (mu1/s)^2;          % compute detection index for desired Pfa and Pd1
%     d_dB = 5*log10(det_indx2);
%     %DT = d_dB - 5*log10(T*W) - 5*log10(NumLines);   % Detection threshold
%     DT = d_dB - 5*log10(TW);   % Detection threshold
else
    m = sqrt(det_indx) * s;         % mean level from detection index equation
    tau = norminv(1-Pfa, mu0, s);   % detection threshold under H0 fro Pfa
    Pd = 1 - normcdf(tau, m, s);    % Pd for varying mean levels
    mu1 = interp1(Pd, m, Pd1);      % find noise mean under H1 for desired Pd1
    det_indx2 = (mu1/s)^2;          % compute detection index for desired Pfa and Pd1
    d_dB = 5*log10(det_indx2);
    %DT = d_dB - 5*log10(T*W) - 5*log10(NumLines);   % Detection threshold
    DT = d_dB - 5*log10(TW);   % Detection threshold
end
%%
% Compute Transmission Loss as function of Range and Frequency
TL = TransmissionLoss(f/1000, R, TL_flag);
TLc = TransmissionLoss(f/1000, R, 'c');

% Compute Ambient Noise Level
NL = NL0 + (OceanNoise(f, SS, ShipLevel, RainLevel, OT_flag, plot_flag))';
% NL = repmat(NL, 1, length(R));

% Compute receiver parameters
lambda = c/fd;      % wavelength, feet
d = lambda/2;       % Element separation in receiver; feet
dx = d; dy = d;     % equal separation in x and y
Lx = (N - 1) * dx;  % array length in x
Ly = (M - 1) * dy;  % array length in y

% Compute array Directivity Index
% [DI,Ns,Ms] = DI_PlanarArray([Lx], Ly, f, fd, c, phi_s,'A');
[DI] = DI_PlanarArray(N, M, f, fd, c, phi_s,'S');

[DI2] = DI_PlanarArray(17, 3, f, fd, c, phi_s,'S');

% SL at receiver input
% SL_r0 = SL - TL - NL;	% SL at receiver hydrophone
SL_r0 = bsxfun(@plus,SL - NL, - TL );	% SL at receiver hydrophone

% Signal Excess at broadband detector output
% SE_r0 = SL - NL + DI - TL - DT; % Passive Signal Excess Equation
SE_r0 = bsxfun(@plus,(SL - NL + DI), - TL) - DT; % Passive Signal Excess Equation
SE_r02 = bsxfun(@plus,(SL - NL + DI2), - TL) - DT; % Passive Signal Excess Equation

% SEc_r0 = SL - TLc - NL + DI - DT;
SEc_r0 = bsxfun(@plus,SL - NL + DI, - TLc) - DT;

%SE_r0 = SL - TL - NL - DT;  % SE - DI
% DI_hat  = DT - SL + TL + NL;
DI_hat  = DT + bsxfun(@plus,-SL + NL, TL);
DI_hat(DI_hat < 0) = 0;

% DIc_hat  = DT - SL + TLc + NL;
DIc_hat  = DT + bsxfun(@plus,-SL + NL, TLc);
DIc_hat(DIc_hat < 0) = 0;

FOM = SL - NL + DI - DT;    % Figure of Merit
df = 100;
f_l = fc;
[ SL_T,NL_T,TW ] = PBD( f,f0,f_l,df,SL(:,1),NL(:,1));



%%
% Plot results
figure
semilogx(f/1000, SL,...
         f/1000, NL,...
         f/1000, TL);
v = axis;
axis([v(1:3) SL1k+10])
legendinfo = cell(1,length(R));
for idx = 1:length(R)
    legendinfo{idx} = ['TL, R = ' num2str(R(idx)/1e3) ' kyds'];
end
legendinfo = [{'SL'},{['NL = SS ' num2str(SS) ' Ship Level ' num2str(ShipLevel)]},legendinfo];
legend(legendinfo,'Location','southwest','FontSize',18)
xlabel('Frequency, KHz','FontSize',18)
ylabel('Level, dB // {\mu}Pa @1 yd','FontSize',18)
title('Source and Noise Spectra, Expected Spherical Transmission Loss','FontSize',18)

%%
figure
semilogx(f/1000, SL_r0)
xlabel('f, KHz','FontSize',18)
ylabel('SL{_r} , dB // {\mu}Pa/Hz @1 yd','FontSize',18)
title(['Signal to Noise Ratio at Hydrophone: SL-TL-NL'],'FontSize',18)
legendinfo = cell(1,length(R));
for idx = 1:length(R)
    legendinfo{idx} = ['R = ' num2str(R(idx)/1e3) ' kyds SS' num2str(SS)];
end
legend(legendinfo,'Location','southwest','FontSize',18)


%%
figure
legendinfo = cell(1,length(R));
% semilogx(f/1000, SE_r0,'-')
% for idx = 1:length(R)
%     legendinfo{idx} = ['R = ' num2str(R(idx)/1e3) ' kyds Spherical Spreading'];
% end
% hold on
semilogx(f/1000, SE_r02,'-');
for idx = 1:length(R)
    legendinfo{idx} = ['R = ' num2str(R(idx)/1e3) ' kyds Spherical Spreading'];
end
legend(legendinfo,'location','Northwest','FontSize',24)
hold on
semilogx(f/1000, zeros(1,length(f)),'k')
v = axis;
axis([v(1:3) v(4)+10])
xlabel('f, KHz','Fontsize',24)
ylabel('SE, dB // ${\mu}$Pa/Hz @1 yd','Fontsize',24)
title(['Signal Excess'],'FontSize',24)

    
%%
figure
semilogx(f/1000, DI,...
         f/1000, DI_hat);
xlabel('f, kHz','FontSize',18)
ylabel('DI, dB','FontSize',18)
title('Directivity Index for Planar Array','FontSize',18)
legendinfo = cell(1,length(R));
for idx = 1:length(R)
    legendinfo{idx} = ['DI Estimate R = ' num2str(R(idx)/1e3) ' kyds'];
end
legendinfo = [{['DI, f${_d}$ = ' num2str(fd/1e3) ' kHz, N = ' num2str(N) ', M = ' num2str(M)]},legendinfo];
legend(legendinfo,'Location','Northwest','FontSize',18);

%% Computing necessary aperture area for a set of ranges
R = linspace(10000,40000,100);     % Range vector
fd = [4000 6000];                   % Design Frequencies
TL_fd = TransmissionLoss(fd./1000, R, TL_flag); % Transmission Loss at fd
NL_fd = NL0 + ...
        (OceanNoise(fd, SS, ShipLevel, RainLevel, OT_flag, plot_flag))';
SL_fd = zeros(1,length(fd));        % Source level at fd
Ly_fixed = [Ly Ly*2];               % Aperture length in y-dir.

% Declaring variables
DI = zeros(length(fd),length(R));   % DI(fd,R)
Aperture = zeros(length(Ly_fixed), length(R), length(fd));  % Aperture Area
N = zeros(length(Ly_fixed), length(R), length(fd));     % Number of elements
Lx = zeros(length(Ly_fixed), length(R), length(fd));    % aperture length
for idf = 1:length(fd)
    SL_fd(idf) = SL(f == fd(idf));  % Loop over design frequency
    % Sonar Equation for Directivity Index assuming SE = 0
    % DI(fd,R) = TL(fd,R) + DT + NL(fd,R) - SL(fd)
    DI(idf,:) = TL_fd(idf,:) + DT + NL_fd(idf) - SL_fd(idf);
    DI(DI<3) = 3;
    DIlinear = 10.^(DI/10);
    for idy = 1:length(Ly_fixed)    % Loop over aperture length in y-dir.
        lambda = c/fd(idf);         % Design frequency wavelength
        dx = lambda/2; dy = lambda/2;   % assume Nyquist spacing
        M = floor(Ly_fixed(idy)/dy);    % Number of elements in y-dir.
        if M <= 1
            M = 2;  % Lower bound on M
        end
        % DI = pi/3 * 2*pi * cos(phi_s) * dx*dy/lambda^2 * sqrt(N^2 -1) *
        % sqrt(M^2 -1) Taken from Directivity TM by A. Nuttall pg. 20
        
        % Need to solve for required N, number of elements in x-direction
        N(idy,:,idf) = ceil(sqrt( (DIlinear(idf,:).*3.* lambda.^2./ ...
                            (sqrt(M^2-1)*2*pi^2*cos(phi_s)*dx*dy)).^2 + 1));
        Lx(idy,:,idf) = (N(idy,:,idf) -1) * dx; % Required length in x-dir.
        Aperture(idy,:,idf) = (N(idy,:,idf)-1) * dx * (M-1) * dy;
        
    end
end

plotLx1 = squeeze(Lx(1,:,:));
plotLx2 = squeeze(Lx(2,:,:));
figure
[ax,h1,h2] = plotyy(R./1e3,DI,R./1e3,plotLx1);
set(ax,'NextPlot','add','Fontsize',24)
hold on
plot(ax(2),R./1e3,plotLx2,'--');
text(1.25, 60, 'Sea State','FontSize',18)

xlabel('Range, kyds','FontSize',24);
ylabel('DI, [dB]','FontSize',24)
% legend(legendinfo,'Location','Northwest','FontSize',18);
% title(['Required DI and Aperture area vs. Range for various f_d and L_y'],'FontSize',24)
legend({['DI, $f_d$ = ' num2str(fd(1)/1e3) 'kHz'],...
        ['DI, $f_d$ = ' num2str(fd(2)/1e3) 'kHz'],...
        ['$L_y$ = ' num2str(Ly_fixed(1)) ' ft, $f_d$ = ' num2str(fd(1)/1e3) 'kHz'],...
        ['$L_y$ = ' num2str(Ly_fixed(1)) ' ft, $f_d$ = ' num2str(fd(2)/1e3) 'kHz'],...
        ['$L_y$ = ' num2str(Ly_fixed(2)) ' ft, $f_d$ = ' num2str(fd(1)/1e3) 'kHz'],...
        ['$L_y$ = ' num2str(Ly_fixed(2)) ' ft, $f_d$ = ' num2str(fd(2)/1e3) 'kHz']},...
        'Fontsize',24,'Location','Northwest')
str = 'L_x, ft';dim = [.925 .35 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Linestyle','none','Fontsize',24,'Fontname','Cambria');
   
   