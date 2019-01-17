% script to call functions for Leviathan IPT support
% RB 9MAR2018. RB 15MAR2018.
clear
present(0)
close all
clc
profile on

% Generate array element locations, parameters, and element time series.
bbflag = 1;
depth = -600*0.3048;    % depth [ft*m/ft]
% Ranges to evaluate: [5] (4572), [10] (9144), [20] (18288) [kYds] (m)
ranges = [4572, 9144, 18288];
% ranges = 4572;
thetaS = -asind(depth./ranges);
% Azimuth locations: 0, 15, 30, 45
phiS = [0, 15, 30, 45];
% phiS = 0;

%% Generate Element Time Series Data
tic
for rInd = 1:length(ranges)
    for phiInd = 1:length(phiS)
        [ETS, p] = GenETS(bbflag,depth,ranges(rInd),phiS(phiInd));
        
        %% Generate array geometry beampattern.
        phi_vec_BP      = -180:180; % azimuth vector for beampattern plots
        theta_vec_BP    = 0:90;     % elevation vector for beampattern plots
        f0_BP           = 680;      % frequency of beampattern
        [BP] = H_VS_Beampattern(p,phi_vec_BP,theta_vec_BP,f0_BP,phiS(phiInd),thetaS(rInd));
        
        %% Process the element time series.
        phi_search = -180:30:180;                       % search grid for azimuth bearings
        range_search = abs(depth):abs(depth)*10:(abs(depth)+30)*100;    % search grid for ranges
        theta_search = -asind(p.depth./range_search);   % consequential search grid for elevation
        MCTs = 100000;
        SNRin = -15;
        
        [ VSoutput, Houtput ] = ProcessETS(bbflag, SNRin, phi_search, theta_search, MCTs, ETS, p);
        if bbflag
            filename = ['BBresults_range' num2str(ranges(rInd)) '_phi' num2str(phiS(phiInd))];
        else
            filename = ['NBresults_range' num2str(ranges(rInd)) '_phi' num2str(phiS(phiInd))];
        end
        save(filename);
        toc
    end
end
profile report