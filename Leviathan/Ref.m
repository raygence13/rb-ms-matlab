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

close all
clear all

%-------------------SETUP PARAMETERS---------------------------------------
% Source Level model
fc = 1000;     % corner frequency for SL change
SL1k = 140;    % SL below fc Hz (fixed from f(1) : 1000)

% Frequency band
f0 = 10000;     % maximum frequency (receiver design), Hz
f = 20 : 1 : f0;    % frequency range, Hz

% Array spacing and DI
fd = 6000;      % desing frequency of receiver, Hz
c = 4900;       % sound speed, feet per sec

N = 36;         % Number of elements is one plane of panel receiver
M = 6;         % Number of elements is one plane of panel receiver
phi_s = 0;      % steering direction for DI calculation; 0 deg == broadside

% Range interval for transmission Loss computation
R = [10000 40000]; % Range, yards

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

% Compute the Detection Threshold
m = sqrt(det_indx) * s;     % mean level from detection index equation
tau = norminv(1-Pfa, mu0, s);    % detection threshold under H0 fro Pfa
Pd = 1 - normcdf(tau, m, s);     % Pd for varying mean levels 
mu1 = interp1(Pd, m, Pd1);	% find noise mean under H1 for desired Pd1
det_indx2 = (mu1/s)^2;	% compute detection index for desired Pfa and Pd1
d_dB = 5*log10(det_indx2);
%DT = d_dB - 5*log10(T*W) - 5*log10(NumLines);   % Detection threshold
DT = d_dB - 5*log10(TW);   % Detection threshold

% Compute Transmission Loss as function of Range and Frequency
TL = TransmissionLoss(f/1000, R, TL_flag);
TLc = TransmissionLoss(f/1000, R, 'c');

% Compute Ambient Noise Level
NL = NL0 + (OceanNoise(f, SS, ShipLevel, RainLevel, OT_flag, plot_flag))';
NL = repmat(NL, 1, length(R));

% Compute receiver parameters
lambda = c/fd;      % wavelength, feet
d = lambda/2;       % Element separation in receiver; feet
dx = d; dy = d;     % equal separation in x and y
Lx = (N - 1) * dx;  % array length in x
Ly = (M - 1) * dy;  % array length in y

% Compute array Directivity Index
W = ones(N,M);
DI = (DI_PlanarArray(N, M, f, fd, c, phi_s,'S'));
DI = repmat(DI, 1, length(R));

SL = repmat(SL, 1, length(R));

% SL at receiver input
SL_r0 = SL - TL - NL;	% SL at receiver hydrophone
% SIgnal Excess at broadband detector output
SE_r0 = SL - TL - NL + DI - DT; % Passive Signal Excess Equation
SEc_r0 = SL - TLc - NL + DI - DT;

%SE_r0 = SL - TL - NL - DT;  % SE - DI
DI_hat  = DT - SL + TL + NL;
DIc_hat  = DT - SL + TLc + NL;

for k = 1 : length(f)
    for k1 = 1 : length(R)
        if DI_hat(k, k1) < 0
            DI_hat(k, k1) = 0;
        end
    end
end

DIc_hat  = DT - SL + TLc + NL;
for k = 1 : length(f)
    for k1 = 1 : length(R)
        if DIc_hat(k, k1) < 0
            DIc_hat(k, k1) = 0;
        end
    end
end

FOM = SL - NL + DI - DT;    % Figure of Merit

% Plot results
figure
set(gcf, 'Position', [15 590 725 375])
semilogx(f/1000, TL(:, 1), f/1000,TL(:, end), f/1000, NL(:,1), ...
         f/1000, SL(:,1),'Linewidth',1.5);
v = axis;
axis([v(1:3) SL1k+10])
     legend(['TL, R = ' num2str(R(1))],['TL, R = ' num2str(R(end))], ...
       ['NL = SS ' num2str(SS) ' Ship Level ' num2str(ShipLevel)], 'SL', 'Location', 'NorthWest')
xlabel('Frequency, KHz','Fontsize',12)
ylabel('Level, dB // {\mu}Pa @1 yd','Fontsize',12)

figure
set(gcf, 'Position', [760 590 725 375])
semilogx(f/1000, SL_r0(:, 1),f/1000, SL_r0(:, end), 'Linewidth',1.5);
%axis([0 100 -160 20])
xlabel('f, KHz','Fontsize',12)
ylabel('SL{_r} , dB // {\mu}Pa/Hz @1 yd','Fontsize',12')
title(['Signal to Noise Ratio at Hydrophone: SL-TL-NL'])
legend(['R = ' num2str(R(1)) ' yds SS' num2str(SS)], ...
       ['R = ' num2str(R(end)) ' yds SS' num2str(SS) ], ...
       'Location', 'SouthWest')

figure
set(gcf, 'Position', [350 125 725 375])
semilogx(f/1000, SE_r0(:, 1), 'b', f/1000, SEc_r0(:, 1), 'b--', ...
         f/1000, SE_r0(:, end), 'g', f/1000, SEc_r0(:, end), 'g--', ...
         'Linewidth',1.5)
v = axis;
axis([v(1:3) v(4)+10])
xlabel('f, KHz','Fontsize',12)
ylabel('SE, dB // {\mu}Pa/Hz @1 yd','Fontsize',12')
%title(['Signal Excess'])
legend(['R = ' num2str(R(1)) ' yds Spherical Spreading'], ...
       ['R = ' num2str(R(1)) ' yds Cylindrical Spreading'], ...
       ['R = ' num2str(R(end)) ' yds Spherical Spreading'], ...
       ['R = ' num2str(R(end)) ' yds Cylindrical Spreading'], ...
        'Location', 'NorthWest')
% legend(['R = ' num2str(R(1)) ' yds SS' num2str(SS)], ...
%        ['R = ' num2str(R(end)) ' yds SS' num2str(SS) ' Ship Level ' num2str(ShipLevel)], ...
%         'Location', 'NorthWest')

    
figure
semilogx(f/1000, DI_hat, f/1000, DI(:,1), 'Linewidth',1.5); 
set(gcf, 'Position', [550 125 725 375])
xlabel('f, KHz','Fontsize',12)
ylabel('DI, dB','Fontsize',12)
%title('Directivity Index for Planar Array')
legend(['DI Estimate R = ' num2str(R(1)) ' yds'], ...
       ['DI Estimate R = ' num2str(R(end)) ' yds'], ...
       ['DI, f${_d}$ = ' num2str(fd) ' Hz, N = ' num2str(N) ', M = ' num2str(M)], ...
        'Location', 'NorthWest')

