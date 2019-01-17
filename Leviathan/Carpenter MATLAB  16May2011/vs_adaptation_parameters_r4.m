function ap = vs_adaptation_parameters_r4(EST, integration_time_secs)
%
% 22 March 2011 - rev 4 includes:
% 1. Accelerometer scale factors derived from DP tests 16-18 March 2011
% 2. Expanded high frequency from 4000 to 8000 Hz
% 3. New hydrophone sensitivity -165 for VS 209
% 4. New Accelerometer sensitivity 1.5 for VS 209
% 5.  Removes subsampling by 2

%% inputs:
% EST - eastern standard time (1) or daylight standard time (0)
% integration_time_secs

%% Integration & Analysis time
ap.integration_time_secs = integration_time_secs;
ap.analysis_duration_secs = ap.integration_time_secs;

%% Data collection system characteristics
ap.fs = 25000;
ap.data_collection_gain_db = 20;

%% Decimation factor (decimation performed in function
%% preprocess_vs_data_r2)
ap.ndec = 1;

%% Some factors
ap.ft2m = 0.304800609601219;
% Rho and c
ap.c = 1460;
ap.rho = 1.026;
ap.hyd_sens = -165; % dB//1V/1microPa
ap.accel_sens = 1.5; % V/g
ap.rho_c_factor = 20*log10(ap.c/1000*ap.rho*1e12); % See Ben Cray conversions sheet
ap.grav = 9.8;
ap.grav_factor = 20*log10(ap.grav); % See Ben Cray conversions sheet

%% Accelerometer scale factor for x and z channels (from DP 16-18 March 2011)
ap.accel_scale_factor = 1.23;

%% BPF characteristics
ap.fll = 2300;
ap.fhh = 8000;
ap.bw = ap.fhh - ap.fll;
ap.fc = (ap.fhh + ap.fll)/2;
ap.fs_ndec = ap.fs/ap.ndec; % Note: subsampling by 2 accounted for here
fl = ap.fll/(ap.fs_ndec/2);
fh = ap.fhh/(ap.fs_ndec/2);
ap.fir_length = 512;
ap.b_bpf = fir1(ap.fir_length,[fl fh]);

%% UTC ref: 1 January 1970 0:00 (see
% http://www.ldas-sw.ligo.caltech.edu/ligotools/faq/convert_time.html)
% Fractional number of days between UTC reference and Matlab reference
% UUV guys use 1 January 1970 0:00 as reference point
% Note: datenum('Jan-1-0000 00:00:00') returns the number 1
ap.utc_ref = datenum('1-Jan-1970', 'dd-mmm-yyyy'); % UTC ref: 1 January 1970 0:00 (see
%                          http://www.ldas-sw.ligo.caltech.edu/ligotools/faq/convert_time.html)
%                          Fractional number of days between UTC reference
%                          and Matlab reference
switch EST
    case 0
        ap.GMT_hr_diff = 4;      % Time difference in hours between local time and GMT
        %                       (this is 4 hours during daylight savings time and 5
        %                        hours during standard time)
    case 1
        ap.GMT_hr_diff = 5;      % Time difference in hours between local time and GMT
        %                       (this is 4 hours during daylight savings time and 5
        %                        hours during standard time)
end

%% CBF (1D setup)
%% Perform conventional BF with sampled spectral matrix
%% (ssm)
%% See m-file PMS394_vector_sensor_angle_estimates_ATF_Oct_2008
%% Beamforming Angles
%% See NAS report Appendix A 
ap.nbeams_theta_1D = 1; % Vertical angles
ap.nbeams_phi_1D = 360; % Horizontal angles
ap.phi_steer_deg_1D = linspace(-180,179,ap.nbeams_phi_1D);
ap.phi_steer_1D = ap.phi_steer_deg_1D*pi/180;
ap.theta_steer_deg_1D = 0;                  %linspace(0,180,nbeams_theta);
ap.theta_steer_1D = ap.theta_steer_deg_1D*pi/180;
ap.V_1D = zeros(4,ap.nbeams_phi_1D);

%% Thresholds
% Added divide by 2 on 3/3/2011 because of change to computation of dB
% level of bf output in main m-file
ap.front_back_threshold_db = 4/2; % Front to back threshold
ap.LOS_rate_threshold_db = 3; % Additional threshold for LOS rate calculation

%% Beamformer dimension
ap.cbf_1D2D_flag = 1;

%% Beam interpolation flag (0 = no interpolation, 1 = interpolation)
ap.beam_interp_flag = 1;

%%
ap.stern_angle_limits_degs = [-175 175];

for j1 = 1:ap.nbeams_theta_1D
    for j2 = 1:ap.nbeams_phi_1D
        ap.a_1D = [cos(ap.theta_steer_1D(j1)).*cos(ap.phi_steer_1D(j2));...
            cos(ap.theta_steer_1D(j1)).*sin(ap.phi_steer_1D(j2));...
            -sin(ap.theta_steer_1D(j1))];
        ap.v_1D = [1;ap.a_1D(1);ap.a_1D(2);ap.a_1D(3)];
        ap.V_1D(:,j2) = ap.v_1D;
    end
end

%% CBF (2D setup)
ap.nbeams_theta_2D = 180; % Vertical
ap.nbeams_phi_2D = 360; % Horizontal
ap.phi_steer_deg_2D = linspace(-180,179,ap.nbeams_phi_2D);
ap.phi_steer_2D = ap.phi_steer_deg_2D*pi/180;
ap.theta_steer_deg_2D = linspace(90,-90,ap.nbeams_theta_2D);
ap.theta_steer_2D = ap.theta_steer_deg_2D*pi/180;
for j1 = 1:ap.nbeams_theta_2D
    for j2 = 1:ap.nbeams_phi_2D
        ap.a_2D = [cos(ap.theta_steer_2D(j1)).*cos(ap.phi_steer_2D(j2));...
            cos(ap.theta_steer_2D(j1)).*sin(ap.phi_steer_2D(j2));...
            -sin(ap.theta_steer_2D(j1))];
    end
end

%% Vector sensor angle mount deltas.  These are determined by comparing vs
%% attitude data derived from the vs serial data streams and the attidue of
%% MARV.
% E.g., head = (heading_degs_mean + heading_delta_degs)*pi/180;
heading_delta_degs = 0;
pitch_delta_degs = 0;%2.5;%2.5;
roll_delta_degs = 0;%13;%13;
ap.heading_delta_degs = heading_delta_degs;
ap.pitch_delta_degs = pitch_delta_degs;
ap.roll_delta_degs = roll_delta_degs;




