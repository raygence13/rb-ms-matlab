

testName = 'Leviathan4xVS';     % Data Directory name
TowScope = 0;
Fs = 20e3;      % Sampling rate of array
%
% Setup some baseline parameters then build paths
%
basedir = '\Data\20Mar2015\';
cd(basedir)
addVstapaths
addpath('C:\Users\radienxe.bautista\Documents\MATLAB\Leviathan\AAXS_3.5')
addpath('C:\Users\radienxe.bautista\Documents\MATLAB\Leviathan\AAXS_3.5\ASIM')
addpath('C:\Users\radienxe.bautista\Documents\MATLAB\Leviathan\AAXS_3.5\Common')
addpath('C:\Users\radienxe.bautista\Documents\MATLAB\Leviathan\AAXS_3.5\UPSIM')
addpath('C:\Users\radienxe.bautista\Documents\MATLAB\Leviathan\AAXS_3.5\VSPP')
TestDirName = fullfile(basedir,filesep,testName);
mydir = pwd;

%% AEL Data Creation
RArrayLayout    = make_empty_rarraylayout;
RSArray         = make_empty_rsarray;
RNArray         = make_empty_rnarray;

mydir = pwd;
cd(fullfile(basedir,filesep,'AAXS_3.7.2',filesep,'ASIM',filesep,'asim_layout'))
Leviathan_array_layout;    % need to generate .m file to define parameters
cd(mydir)

% Fill RArrayLayout structure
RArrayLayout.array_name = 'Leviathan4xVSArray';
RArrayLayout.oflag = 0;
RArrayLayout.cable_file = fullfile(basedir,filesep,'SimParams',filesep,'VSTA13Xcable.txt');

% % Fill in RNArray parameters with 4x phone spacing
RNArray(1,1).isa = 1;
RNArray(1,1).type = 'h';
RNArray(1,1).Space_lE = or_sensors_locations;
RNArray(1,1).ls = 1;
RNArray(1,1).dl = 0;

RNArray(2,1).isa = 1;
RNArray(2,1).type = 'p'; % Pitch
RNArray(2,1).Space_lE = or_sensors_locations;
RNArray(2,1).ls = 1;
RNArray(2,1).dl = 0;

RNArray(3,1).isa = 1;
RNArray(3,1).type = 'r'; % Roll
RNArray(3,1).Space_lE = or_sensors_locations;
RNArray(3,1).ls = 1;
RNArray(3,1).dl = 0;

RNArray(4,1).isa = 1;
RNArray(4,1).type = 'h'; % Heading sensors on each subarray
RNArray(4,1).Space_lE = heading_sensors_locations;
RNArray(4,1).ls = 1;
RNArray(4,1).dl = 0;

RNArray(5,1).isa = 1;
RNArray(5,1).type = 'z'; % Depth sensors on each subarray
RNArray(5,1).Space_lE = depth_sensors_locations;
RNArray(5,1).ls = 1;
RNArray(5,1).dl = 0;

% Put in acoustic sensor locations
RSArray(1,1).Space_lE = hydrophone_locations;
RSArray(1,1).Tflag = [1 0 0 0]; % Pressure only
RSArray(2,1).Space_lE = or_sensors_locations;
RSArray(2,1).Tflag = [0 1 1 1]; % in-line, x-line hor, x-line ver

RSArray(1,1).dxc = 0; RSArray(1,1).dyc = 0; RSArray(1,1).dzc = 0;
RSArray(1,1).nhl = 1; RSArray(1,1).nvl = 1; RSArray(1,1).dxp = 0;
RSArray(1,1).dxil = 0; RSArray(1,1).oazil = 0; RSArray(1,1).odeil = 0; RSArray(1,1).dxxlh = 0;
RSArray(1,1).oazxlh = 90; RSArray(1,1).odexlh = 0;
RSArray(1,1).dxxlv = 0; RSArray(1,1).oazxlv = 0;
RSArray(1,1).odexlv = 90;

RSArray(2,1).dxc = 0; RSArray(2,1).dyc = 0; RSArray(2,1).dzc = 0;
RSArray(2,1).nhl = 1; RSArray(2,1).nvl = 1; RSArray(2,1).dxp = 0;
RSArray(2,1).dxil = 0; RSArray(2,1).oazil = 0; RSArray(2,1).odeil = 0; RSArray(2,1).dxxlh = 0;
RSArray(2,1).oazxlh = 90; RSArray(2,1).odexlh = 0;
RSArray(2,1).dxxlv = 0; RSArray(2,1).oazxlv = 0;
RSArray(2,1).odexlv = 90;

rev = cur_aaxs_rev;

if ~isdir(TestDirName)
    mkdir(TestDirName);
end
save(fullfile(TestDirName,filesep,'Base1xLayout'),'RArrayLayout', 'RNArray', 'RSArray','rev');

% Now fill out simulation control structure
mydir = pwd;
cd(fullfile(basedir, filesep,'SimParams'));
Leviathan4xVS_params

cd(mydir);

% Array Motion Parameters:
% Simulation Parameters
RArraySim.array_sim_name = testName;
RArraySim.tc = simulation_midpoint;
RArraySim.dsim = simulation_end;
RArraySim.motion = 'cable'; % uses cable model for motion
% cable and tow point motion parameters:
RACable.rref = 0;   % Range from array midpoint to tow ship
RACable.bref = 1;   % Initial array heading
RACable.T_Ti = time;    % time vector created in params script
RACable.S_Ti = speed;   % speed vector created in params script
RACable.C_Ti = course;  % course vector created in params script
RACable.Z_Ti = depth;   % depth vector created in params script

sim_control_filename = fullfile(TestDirName,filesep,'Sim_Control.mat');
save(sim_control_filename,'RArraySim','RACable','rev');

% Save all array layout stuff in a file for call to make_array_sim:
array_layout_filename = fullfile(TestDirName,filesep,['Base1xLayout.mat']);
save(array_layout_filename,'RArrayLayout','RNArray','RSarray','rev');

% make array shape simulation
mydir = pwd;
cd(TestDirName);
mkdir('ArraySim');
cd ArraySim; % must be in different directory than input files
make_asim_array(array_layout_filename, sim_control_filename); % Writes _array.mat file to TestDirName
array_filename = fullfile(pwd,filesep,testName,filesep,[testName '_array.mat']);

cd(mydir);
%%%__________________AEL Data Creation_____________%%%
%%%__________________Non-Acoustic_________________________%%%%
nagen_setup % will create the nagen_ctrl.mat
setup_nagen_run    % laod nagen_ctrl.mat and _nas.mat from ASIM

% Never do this again unless you want a new tow geometry

%% _______________ NOISE Data Creation
% This cell is for making a specific interferer noise projected onto the
% AEL we just defined
% flcose all; close all; clear all;
pack;

% Reload base names and set path
%testName = '';
% Setup some baseline parameters then build paths

basedir = '/Data/20Mar2015';
TestDirName = fullfile(basedir,filesep,testName);
array_filename = fullfile(TestDirName,filesep,'ArraySim',filesep,testName,filesep,[testName '_array.mat']);
% sim_control_filename = fullfile(TestDirName,filesep,'Sim_Control.mat');
mydir = pwd;
cd(fullfile(basedir,filesep,'AAXS_3.7.2'));
%setup_aaxs_paths;
cd(mydir);
addpath(fullfile(basedir,filesep,'VSTAXL_2.1'));

% Now fill out simulation control structure
mydir = pwd;
cd(fullfile(basedir,filesep,'SimParams'));
%vsta13x_uturn_params; %Creates time, course, speed, and depth variables
%for a U-turn 10 minute simulation
Leviathan4xVS_params
cd(mydir);

% Prepare to make some noise...
RDgen = make_empty_rdgen;
RDgen.sim_name = testName;
RDgen.tc = simulation_midpoint;
RDgen.dsim = simulation_end;
RDgen.tr = Fs;
RDgen.df = 0.5; % frame interval (sec). frequency resolution with go as 1/RDgen.df
RDgen.dw = 30; % simulation window length
RDgen.vssflag = 2; % vector sensors are accelerometers (m/s^2)
RDgen.array_file = array_filename;

RDgen.nrcf = 8; % Only recompute noise variance every 4 seconds (8 frames)
rev = cur_aaxs_rev;

RDgen.nDist.nfname = 'OmniNoise10dBperdecade';
RDgen.nDist.nftype = 'iso';
RDgen.nDist.ng = -10.0;
RDgen.nDist.Fn_Fn = [];
RDgen.nDist.noflag = 1;
% Save all noise simulation parameters:
noise_sim_filename1 = fullfile(TestDirName,filesep,'Omni_Noise_Sim_Ctrl.mat');
save(noise_sim_filename1,'RDgen','rev');

mydir = pwd;
cd(TestDirName);
mkdir('UPSim');
cd UPSim; % must be in different directory than input files
make_upsim_ds(noise_sim_filename1); % writes _array.mat file to TestDirName
cd(mydir);
%% ________________ NOISE Data Creation

%% IR Signal Data Creation
% Now make a constant bearing interferer with NB and BB noise components

%%%% Use this section to create targets
clear RDgen;
% prepare to make some noise...
RDgen = make_empty_rdgen;
RDgen.sim_name = [testName '_IR_Int'];
RDgen.tc = simulation_midpoint;
RDgen.dsim = simulation_end;
RDgen.tr = Fs;
RDgen.df = 0.5;
RDgen.dw = 60; % Simulation window length
RDgen.vssflag = 2; % Vector sensors are accelerometers (m/s^2)
RDgen.array_file = array_filename;

RDgen.nrcf = 8; % Only recompute noise variance very 4 seconds (8 frames)
rev= cur_aaxs_rev;

% set up an infinite range source moving through 1 deg in 8 seconds
RDgen.iSrc.name = 'NBSource';
RDgen.iSrc.tref = simulation_midpoint;
RDgen.iSrc.d = 0;
RDgen.iSrc.b = 145;
RDgen.iSrc.dd = 0;
RDgen.iSrc.db = 0;
% Single NB tone at 200 Hz
RDgen.iSrc.NB.name = '60HzTone';
RDgen.iSrc.NB.fc = 60.0;
RDgen.iSrc.NB.bw = 0.0;

% Setup BB from 200 Hz;
RDgen.iSrc.BB.name = 'BBsourceFc200HzSig15Hz';
RDgen.iSrc.BB.Fi_Fi = [184:2:216];
RDgen.iSrc.BB.SLdB_Fi = 75 -10/15*(abs(200-RDgen.iSrc.BB.Fi_Fi)); % 10dB down 15Hz from center

% Save all noise simulation parameters:
noise_sim_filename2 = fullfile(TestDirName,filesep,'IR_Interferer_60HzNB_200HzBB_145deg.mat');
save(noise_sim_filename2,'RDgen','rev');

mydir = pwd;
cd(TestDirName);
mkdir('UPSim');
cd UPSim; %must bein different directory than input files
make_upsim_ds(noise_sim_filename2); % writes _array.mat file to TestDirName
cd(mydir);

%% IR Signal Data Creation


%% FR Signal Data Creation
% Now make a constant bearing interferer with NB and BB noise components
%%%% Use this section to crate targets
clear RDgen;

% Prepare to make some noise...
RDgen = make_empty_rdgen;
RDgen.sim_name = [testName '_FR_Int'];
RDgen.tc = simulation_midpoint;
RDgen.dsim = simulation_end;
RDgen.tr = Fs;
RDgen.df = 0.5;
RDgen.dw = 60;
RDgen.vssflag = 2;
RDgen.array_file = array_filename;

RDgen.nrcf = 8;
rev = cur_aaxs_rev;

% Set up an infinite range source moving through 1 Deg in 8 seconds
RDgen.tSrc.name = '65DegNBplusBB';
RDgen.tSrc.tref = simulation_midpoint;
RDgen.tSrc.r = 3000; % Range
RDgen.tSrc.b = 65; % Constant bearing 65 deg starboard of ownship heading
RDgen.tSrc.z = 50; % depth
RDgen.tSrc.s = 3*0.5144;    % speed
RDgen.tSrc.c = 180; % course
RDgen.tSrc.w = 0;

% three NB Tone
RDgen.tSrc.NB(1).name = '75HzTone';
RDgen.tSrc.NB(1).fc = 75.0;
RDgen.tSrc.NB(1).bw = 0.8;

RDgen.tSrc.NB(2).name = '110HzTone';
RDgen.tSrc.NB(2).fc = 110.0;
RDgen.tSrc.NB(2).bw = 0.1;


RDgen.tSrc.NB(3).name = '300HzTone';
RDgen.tSrc.NB(3).fc = 300.0;
RDgen.tSrc.NB(3).bw = 3.0;

% Set up BB from 250 Hz;

RDgen.tSrc.BB.name = 'BBsourceFc250HzSig40Hz';
RDgen.tSrc.BB.Fi_Fi = [230:2:270];
RDgen.tSrc.BB.SLdB_Fi = 75 - 10/15*(abs(200-RDgen.tSrc.BB.Fi_Fi));

% Save all noise simulation parameters
noise_sim_filename2 = fullfile(TestDirName,filesep,'FR_Interferer_3NBtones_250HzBB_65Deg.mat');
save(noise_sim_filename2,'RDgen','rev');

mydir = pwd;
cd(TestDirName);
mkdir('UPSim')
cd UPSim; % Must be in different directory than input files
make_upsim_ds(noise_sim_filename2); % writes _array.mat file to TestDirName
cd(mydir);
%% FR Signal Data Creation
setup_SDS_run

%% convert AAXS data to VDS files
aaxs2vds13x

% Scale noise 70, IR BB 58 & NB72, FR BB 53 NB2 65 NB2 56 NB3 45
