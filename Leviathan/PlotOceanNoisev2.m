% This function will generate sea noise based on values taken from "Ambient
% Noise Standards for Acoustic Modeling and Analysis", Walt Sadowski, NUSC
% Techincal Document 7265, 1984.  This data is saved in a file called
% "SeaNoise.txt".
% PlotOcenNoisev2 was modified to plot all values down to 10 Hz.

clear all
%close all
X = load('SeaNoise.txt');           % load data
f = X(:, 1);                        % Frequencies Hz
SS0 = X(:, 2);  idx2=find(SS0~=0);  % Sea State 0
SS1 = X(:, 3);  idx3=find(SS1~=0);  % Sea State 1
SS2 = X(:, 4);  idx4=find(SS2~=0);  % Sea State 2
SS3 = X(:, 5);  idx5=find(SS3~=0);  % Sea State 3
SS4 = X(:, 6);  idx6=find(SS4~=0);  % Sea State 4
SS5 = X(:, 7);  idx7=find(SS5~=0);  % Sea State 5
SS6 = X(:, 8);  idx8=find(SS6~=0);  % Sea State 6
SL1 = X(:, 9);  idx9=find(SL1~=0);  % Shipping Level 1
SL2 = X(:, 10);  idx10=find(SL2~=0);  % Shipping Level 2
SL3 = X(:, 11);  idx11=find(SL3~=0);  % Shipping Level 3
SL4 = X(:, 12);  idx12=find(SL4~=0);  % Shipping Level 4
SL5 = X(:, 13);  idx13=find(SL5~=0);  % Shipping Level 5
SL6 = X(:, 14);  idx14=find(SL6~=0);  % Shipping Level 6
SL7 = X(:, 15);  idx15=find(SL7~=0);  % Shipping Level 7
OT = X(:, 16);   idx16=find(OT~=0);   % Ocean Turbulence
RL1 = X(:, 17);  idx17=find(RL1~=0);  % Rain Level, Intermediate
RL2 = X(:, 18);  idx18=find(RL2~=0);  % Rain Level, Medium
RL3 = X(:, 19);  idx19=find(RL3~=0);  % Rain Level, Heavy
TH = X(:, 20);  idx20=find(TH~=0);    % Thermal Noise

% Sea state values below where they peak is not same as Urick.  In Urick it
% drops off at about 4 dB per octave going down in frequency from where the
% peak is located.  Here they are flat from the peak down to 10 Hz

% SS5_1 = interp1(f(idx7), SS5(idx7), 5000);
% SS5_1 = interp1(f(idx7), SS5(idx7), [400 : 100 : 10000]);

% The following is the Beaufort Number - Sea-state relationship
% From Urick Table 1-1, Ambient Noise in The Sea, 1984
%                           Wind Speed
% Beaufort #    Sea State   Knots   m/s
%   0               0       <1      0-0.2
%   1               0.5     1-3     0.3-1.5
%   2               1       4-6     1.6-3.3
%   3               2       7-10    3.4-5.4
%   4               3       11-16   5.5-7.9
%   5               4       17-21   8.0-10.7
%   6               5       22-27   10.8-13.8
%   7               6       28-33   13.9-17.1
%   8               6       34-40   17.2-20.7

figure
set(1,'Position', [15 490 825 500])
semilogx(f(idx2), SS0(idx2), 'b', f(idx3), SS1(idx3), 'g', ...
         f(idx4), SS2(idx4), 'r', f(idx5), SS3(idx5), 'c', ...
         f(idx6), SS4(idx6), 'm', f(idx7), SS5(idx7), 'y', ...
         f(idx8), SS6(idx8), 'k', ...
         f(idx9), SL1(idx9), 'b--', f(idx10), SL2(idx10), 'g--', ...
         f(idx11), SL3(idx11), 'r--', f(idx12), SL4(idx12), 'c--', ...
         f(idx13), SL5(idx13), 'm--', f(idx14), SL6(idx14), 'y--', ...
         f(idx15), SL7(idx15), 'k--',...
         f(idx16), OT(idx16), 'b', f(idx17), RL1(idx17), 'b', ...
         f(idx18), RL2(idx18), 'g', f(idx19), RL3(idx19), 'r', ...
         f(idx20), TH(idx20), 'b', ...
         'LineWidth', 1.5);  grid

text(1.25, 60, 'Sea State','FontSize',14)
text(f(idx2(1))-4, SS0(idx2(1)), 'SS0','Color', 'b')
text(f(idx3(1))-4, SS1(idx3(1)), 'SS1','Color', 'g')
text(f(idx4(1))-4, SS2(idx4(1)), 'SS2','Color', 'r')
text(f(idx5(1))-4, SS3(idx5(1)), 'SS3','Color', 'c')
text(f(idx6(1))-4, SS4(idx6(1)), 'SS4','Color', 'm')
text(f(idx7(1))-4, SS5(idx7(1)), 'SS5','Color', 'y')
text(f(idx8(1))-4, SS6(idx8(1)), 'SS6','Color', 'k')

SL1mx = min(find(SL1 == max(SL1)));
SL2mx = min(find(SL2 == max(SL2)));
SL3mx = min(find(SL3 == max(SL3)));
SL4mx = min(find(SL4 == max(SL4)));
SL5mx = min(find(SL5 == max(SL5)));
SL6mx = min(find(SL6 == max(SL6)));
SL7mx = min(find(SL7 == max(SL7)));
text(20, 95,'Shipping Level','FontSize',14)
text(f(SL1mx), SL1(SL1mx), 'SL1','Color', 'b')
text(f(SL2mx), SL2(SL2mx), 'SL2','Color', 'g')
text(f(SL3mx), SL3(SL3mx), 'SL3','Color', 'r')
text(f(SL4mx), SL4(SL4mx), 'SL4','Color', 'c')
text(f(SL5mx), SL5(SL5mx), 'SL5','Color', 'm')
text(f(SL6mx), SL6(SL6mx), 'SL6','Color', 'y')
text(f(SL7mx), SL7(SL7mx), 'SL7','Color', 'k')

text(1.5, OT(idx16(6)), 'OT','FontSize',14)

text(700, 85, 'Rain Level','FontSize',14)
text(f(idx17(end)), RL1(idx17(end)), 'RL-INTER.','Color', 'b')
text(f(idx18(end)), RL2(idx18(end)), 'RL-MEDIUM','Color', 'g')
text(f(idx19(end)), RL3(idx19(end)), 'RL-HIGH','Color', 'r')

text(4000, TH(idx20(7)), 'THERMAL','FontSize',14)

axis([0 1e5 10 100])
xlabel('Frequency Hz','FontSize',14)
ylabel('Spectrum Level dB re 1 {\mu}Pa','FontSize',14)
title('Ocean Noise Spectra','FontSize',14)
% legend('SS0', 'SS1', 'SS2', 'SS3', 'SS4', 'SS5', 'SS6', ...
%        'SL1', 'SL2', 'SL3', 'SL4', 'SL5', 'SL6', 'SL7', ...
%        'OT','RI', 'RM', 'RH', 'TH', 'Location', 'EastOutside')
