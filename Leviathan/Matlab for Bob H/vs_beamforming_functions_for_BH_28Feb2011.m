% vs_beamforming_functions_for_BH_28Feb2011
%
% 19 October 2010

clear all
fname1 = mfilename;

commandwindow

%% Run number
run_number = 4;

%% Get sound boat GPS and 841 boat gps
    load('L:\_MARV_VS_Run_2_FEB23_2011\CAPS_SOUNDBOAT_GPS_20110223.mat')
    load('L:\_MARV_VS_Run_2_FEB23_2011\Hypack Nav Log 02232011\combined.mat')
    hypack_utc_time = hypack_841_zulu_hrs_to_UTC(2011,2,4,hypack.zulu_time_hrs);

%% Save Data flag (save data for hand-off to bearings only tracking
%% algorithm development)
save_data_flag = 0;

%% Compute power in whole band and band of interest (useful for run
%% analysis)
compute_power_flag = 1;

%% Generate adaptation parameters, etc.
% EST - eastern standard time (1) or daylight savings time (0)
EST = 1;
ap = vs_adaptation_parameters_r3(EST);

%% Get MARV NAV dat
start_folder1 = 'L:\_MARV_VS_Run_2_FEB23_2011\all_data\';
nav_folder1 = start_folder1;
%[nav_filename,nav_folder1] = uigetfile('nav*.mat*','NAV Data',start_folder1);
switch run_number
    case 1
nav_filename = 'nav_data_23-Feb-2011_081540.mat';
    case 2
nav_filename = 'nav_data_23-Feb-2011_110644.mat';
    case 3
nav_filename = 'nav_data_23-Feb-2011_131121.mat';
    case 4
nav_filename = 'nav_data_23-Feb-2011_142919.mat';
    case 5
nav_filename = 'nav_data_23-Feb-2011_150947.mat';
end
load([nav_folder1 nav_filename])

%% Fix utc_time, etc. (MARV typically has some samples where the UTVC times
%% are flakey)
fix_utc_time_threshold_secs = 2;
[nav_data, accept_operation_data, nav_utc_ref_time] = ...
    fix_utc_time_nuke_NaNs_r2(nav_data,accept_operation_data,fix_utc_time_threshold_secs);

%% Determine UTM Zone and convert UUV Lat/Long to UTM from combined file.
%% Requires Matlab Mapping toolbox.
utm_zone = utmzone(nav_data.lat_degs(1),nav_data.lon_degs(1));
% UTM mapping parameters for (see Mapping toolbox)
mstruct = defaultm('utm');
mstruct.zone = utm_zone;
mstruct.geoid = almanac('earth','geoid','meters','wgs84');
mstruct = defaultm(utm(mstruct));
[MARV_east,MARV_north] = mfwdtran(mstruct,nav_data.lat_degs,nav_data.lon_degs);


%% Vector Sensor Folder and file name
vs_folder = start_folder1;
switch run_number
    case 1 % Zig zag with stationary source
vs_filename_base = '20110223081541-0001_caps.dat';
stop_start = [4 113];
caps_nav_offset_secs = 0;
source_depth = 4*0.3048;
    case 2 % Zig zag with mobile source
vs_filename_base = '20110223110645-0001_caps.mat';
stop_start = [4 109];
caps_nav_offset_secs = 0;
source_depth = 4*0.3048;
    case 3 % Parallel source 1
vs_filename_base = '20110223131122-0001_caps.dat';
stop_start = [4 57];
caps_nav_offset_secs = 0;
source_depth = 4*0.3048;
    case 4 % Parallel source 2
vs_filename_base = '20110223142920-0001_caps.dat';
stop_start = [4 28];
caps_nav_offset_secs = 0;
source_depth = 4*0.3048;
    case 5 % Parallel source 3
vs_filename_base = '20110223150948-0001_caps.dat';
stop_start = [4 83];
caps_nav_offset_secs = 0;
source_depth = 4*0.3048;
end

%% Set aside memory
rough_number_of_time_steps = ceil((diff(stop_start)+1)*...
    60/ap.integration_time_secs);
bf_conven_1D_db = NaN*zeros(rough_number_of_time_steps,ap.nbeams_phi_1D);
running_time_secs = NaN*zeros(rough_number_of_time_steps,1);
power_per_interval_save.whole_band = NaN*zeros(rough_number_of_time_steps,1);
power_per_interval_save.bpf = NaN*zeros(rough_number_of_time_steps,1);
cbf_1D_angle_est_degs = NaN*zeros(rough_number_of_time_steps,1);
cbf_level = NaN*zeros(rough_number_of_time_steps,1);
max_stern_bf_conven_db = NaN*zeros(rough_number_of_time_steps,1);
HLOS_degs = NaN*zeros(rough_number_of_time_steps,1);
simple_HLOS_degs = NaN*zeros(rough_number_of_time_steps,1);
HLOS_degs_2 = NaN*zeros(rough_number_of_time_steps,1);
VLOS_degs = NaN*zeros(rough_number_of_time_steps,1);
MARV_east_mean = NaN*zeros(rough_number_of_time_steps,1);
MARV_north_mean = NaN*zeros(rough_number_of_time_steps,1);
heading_degs_mean = NaN*zeros(rough_number_of_time_steps,1);
roll_degs_mean = NaN*zeros(rough_number_of_time_steps,1);
pitch_degs_mean = NaN*zeros(rough_number_of_time_steps,1);
s_plus_n_to_n_ratio_cbf_db = NaN*zeros(rough_number_of_time_steps,1);
phi_est_inten = NaN*zeros(rough_number_of_time_steps,1);
theta_est_inten = NaN*zeros(rough_number_of_time_steps,1);
interference_flag = NaN*zeros(rough_number_of_time_steps,1);

%% Loop over 1 minute (nominal) VS data files
for m = 1:stop_start(1)
    last_file_num = m;
    if last_file_num < 10
        file_num = ['000',int2str(last_file_num)];
    elseif last_file_num >= 10 && last_file_num < 100
        file_num = ['00',int2str(last_file_num)];
    else
        file_num = ['0',int2str(last_file_num)];
    end
    vs_filename = [vs_filename_base(1:15) file_num vs_filename_base(end-8:end-4)];
    load([vs_folder vs_filename])
    if m == 1
        caps_time_start_UTC = caps_utc_time_start(vs_filename, ap);
    else
        pressure1 = pressure.s1;
        caps_time_start_UTC = caps_time_start_UTC + length(pressure1)/ap.fs;
    end
end
relative_caps_utc_time_start = caps_time_start_UTC - nav_utc_ref_time;

kk = 0;
file_cntr = 0;
for m = stop_start(1):stop_start(2)
    disp([m stop_start(2)])
    file_cntr = file_cntr + 1;
    
    last_file_num = m;
    if last_file_num < 10
        file_num = ['000',int2str(last_file_num)];
    elseif last_file_num >= 10 && last_file_num < 100
        file_num = ['00',int2str(last_file_num)];
    else
        file_num = ['0',int2str(last_file_num)];
    end
    vs_filename = [vs_filename_base(1:15) file_num vs_filename_base(end-8:end-4)];
        
    load([vs_folder vs_filename])
    pressure1 = pressure.s1;
    ax1 = ax.s1;
    ay1 = ay.s1;
    az1 = az.s1;
    
    nsamps_int = ap.integration_time_secs*ap.fs;
    file_size_secs = length(pressure1)/ap.fs;
    nloops = floor(length(pressure1)/nsamps_int);
    lost_time_secs = file_size_secs - nloops*ap.integration_time_secs;
    %% Loop over integration periods
    for n = 1:nloops
        kk = kk + 1;
        if kk == 1
            running_time_secs(kk) = relative_caps_utc_time_start;
        else
            running_time_secs(kk) = running_time_secs(kk-1) + ...
                ap.integration_time_secs;
        end
        if n == nloops
            running_time_secs(kk) = running_time_secs(kk) + lost_time_secs;
        end
        
        samps = 1:nsamps_int;
        samps = samps + (n-1)*nsamps_int;
        pressure = pressure1(samps);
        ax = ax1(samps);
        ay = ay1(samps);
        az = az1(samps);
        
        %% Preprocess data
        [dec_time,fs_cmplx,cmplx_env_sensor_signal,power_per_interval] = ...
            preprocess_vs_data_no_channel_mixup(ap,pressure,ax,ay,az,compute_power_flag);
        if compute_power_flag == 1
            power_per_interval_save.whole_band(kk) = power_per_interval.whole_band;
            power_per_interval_save.bpf(kk) = power_per_interval.bpf;
        end
                
        %% Conventional Beamformer (1D)
        bf_conven_1D = cbf_vs_1D(ap,dec_time,cmplx_env_sensor_signal);
        bf_conven_1D_db(kk,:) = 20*log10(abs(bf_conven_1D));
        
        %% Find angle (and level) of max 1D CBF output
        [cbf_1D_angle_est_degs(kk), cbf_level(kk), max_stern_bf_conven_db(kk), ...
            interference_flag(kk), s_plus_n_to_n_ratio_cbf_db(kk)] = ...
            oneD_bf_angle_estimate_with_stern_montoring(ap,bf_conven_1D_db(kk,:));
                
        %% compute_LOS_degs
        data = compute_LOS_degs(nav_data, ...
            MARV_east, MARV_north, running_time_secs(kk), ...
            nav_utc_ref_time, ap, cbf_1D_angle_est_degs(kk), ...
            0, ...
            ap.heading_delta_degs, ap.pitch_delta_degs, ap.roll_delta_degs);
        
        HLOS_degs(kk) = data.HLOS_degs;
        HLOS_degs_2(kk) = data.HLOS_degs_2;
        VLOS_degs(kk) = data.VLOS_degs;
        MARV_east_mean(kk) = data.MARV_east_mean;
        MARV_north_mean(kk) = data.MARV_north_mean;
        roll_degs_mean(kk) = data.roll_degs_mean;
        pitch_degs_mean(kk) = data.pitch_degs_mean;
        heading_degs_mean(kk) = data.heading_degs_mean;
        simple_HLOS_degs(kk) = heading_degs_mean(kk) + cbf_1D_angle_est_degs(kk);
                
    end
end

%% Set thresholds on Signal Plus Noise to Noise Ratio
%% Also determine where acoustic interference occurs
ii_low_cbf_level = find(s_plus_n_to_n_ratio_cbf_db >= ...
    ap.front_back_threshold_db & ...
interference_flag <= (numel(ap.phi_steer_deg_1D))-5);
ii_hi_cbf_level = find(s_plus_n_to_n_ratio_cbf_db >= ...
    ap.front_back_threshold_db + ap.LOS_rate_threshold_db & ...
    interference_flag <= (numel(ap.phi_steer_deg_1D))-5);

%% Sound boat & 841 location during files of interest
min_time = nav_utc_ref_time + running_time_secs(1);
max_time = nav_utc_ref_time + max(running_time_secs);
gps_utc_time = gps.utc_time;
ii_good_gps_time = find(gps_utc_time <= max_time & gps_utc_time >= min_time);
[source_east,source_north] = mfwdtran(mstruct,gps.lat(ii_good_gps_time),gps.lon(ii_good_gps_time));
source_east = interp1(gps_utc_time(ii_good_gps_time),source_east,...
    linspace(min_time,max_time,numel(running_time_secs)));
source_north = interp1(gps_utc_time(ii_good_gps_time),source_north,...
    linspace(min_time,max_time,numel(running_time_secs)));

ii_good_hypack_time = find(hypack_utc_time <= max_time & hypack_utc_time >= min_time);

%% Compute range to source from GPS
range_to_source_m = sqrt((MARV_east_mean-source_east').^2 + (MARV_north_mean-source_north').^2);

%% Compute geometric HLOS and bearing using GPS
HLOS_geometric_degs = 180/pi*atan2(-(MARV_east_mean-source_east'),-(MARV_north_mean-source_north'));
HLOS_geometric_degs_2 = HLOS_geometric_degs;
HLOS_geometric_rate = diff(HLOS_geometric_degs)./diff(running_time_secs);
ii = find(HLOS_geometric_degs < 0);
HLOS_geometric_degs(ii) = HLOS_geometric_degs(ii) + 360;

%% Image CBF Results
if 0
    figure(1); clf
    imagesc(running_time_secs-running_time_secs(1),(ap.phi_steer_deg_1D),...
        (bf_conven_1D_db)'-max(max(bf_conven_1D_db)))
    set(gca,'ydir','rev')
    colorbar
    xlabel('Time (secs)')
    ylabel('Horizontal Angle \eta (degs)')
    title(['CBF Output (Horizontal Plane): ',vs_filename],'interp','none','fontsize',8)
    pltfname1(fname1)
end

figure(3); clf
plot(running_time_secs-running_time_secs(1),cbf_level-max(max(bf_conven_1D_db)),...
    running_time_secs-running_time_secs(1),max_stern_bf_conven_db-max(max(bf_conven_1D_db)),'r')
legend('Max(CBF)','Stern (CBF)')
grid
xlabel('Running Time (secs)')
ylabel('dB')
title(['Front-to-back threshold (dB) = ',num2str(ap.front_back_threshold_db)])
pltfname1(fname1)

figure(4); clf
subplot(211)
plot(running_time_secs-running_time_secs(1),cbf_1D_angle_est_degs,'.',...
    running_time_secs(ii_low_cbf_level)-running_time_secs(1),...
    cbf_1D_angle_est_degs(ii_low_cbf_level),'r.',...
    running_time_secs-running_time_secs(1),heading_degs_mean,'c.')
legend('Ang. Est.','with detection & no interference','Mean(heading)')
set(gca,'ydir','rev')
ylabel('Bearing Estimate (degs)')
grid
title(['Integration Time (secs) = ',num2str(ap.integration_time_secs)])
subplot(212)
plot(running_time_secs-running_time_secs(1),HLOS_degs,'.',...
    running_time_secs(ii_low_cbf_level)-running_time_secs(1),HLOS_degs(ii_low_cbf_level),'r.')
hold on
plot(running_time_secs-running_time_secs(1),HLOS_geometric_degs,'g')
grid
legend('Estimate','Estimate & Det. Occurs','Ground Truth from GPS')
set(gca,'ydir','rev')
xlabel('Running Time (secs)')
ylabel('HLOS Estimate (degs)')
pltfname1(fname1)

figure(5); clf
plot(MARV_east_mean, MARV_north_mean,'k',...
    MARV_east_mean(ii_low_cbf_level), MARV_north_mean(ii_low_cbf_level),'b.',...
    MARV_east_mean(ii_hi_cbf_level), MARV_north_mean(ii_hi_cbf_level),'g.',...
    source_east,source_north,'cp',...
    hypack.boat841_east(ii_good_hypack_time),hypack.boat841_north(ii_good_hypack_time),'r.')
for kkk = 1:100:length(running_time_secs)
    text(MARV_east_mean(kkk), MARV_north_mean(kkk),...
        [' ',num2str(running_time_secs(kkk)-running_time_secs(1),'% 6.0f')])
    text(source_east(kkk),source_north(kkk),...
        [' ',num2str(running_time_secs(kkk)-running_time_secs(1),'% 6.0f')])
end
axis equal
legend('MARV','Detection Occurred','LOS Rate Calc.','Source','841')
xlabel('Easting (m)')
ylabel('Northing(m)')
grid
pltfname1(fname1)

figure(6); clf
subplot(211)
plot(running_time_secs-running_time_secs(1),10*log10(power_per_interval_save.whole_band))
grid
ylabel('Power (dB)')
title(['Whole Band'])
subplot(212)
plot(running_time_secs-running_time_secs(1),10*log10(power_per_interval_save.bpf))
grid
xlabel('Running Time (secs)')
ylabel('Power (dB)')
title([num2str(ap.fll),' - ',num2str(ap.fhh),' Hz'])
pltfname1(fname1)

figure(7); clf
subplot(311)
plot(nav_data.utc_time-caps_time_start_UTC,nav_data.roll_degs,...
running_time_secs-running_time_secs(1),roll_degs_mean,'r.')
legend('All data','Mean during T')
grid
ylabel(' Roll (degs)')
subplot(312)
plot(nav_data.utc_time-caps_time_start_UTC,nav_data.pitch_degs,...
running_time_secs-running_time_secs(1),pitch_degs_mean,'r.')
grid
ylabel('Pitch (degs)')
subplot(313)
plot(nav_data.utc_time-caps_time_start_UTC,nav_data.heading_degs,...
running_time_secs-running_time_secs(1),heading_degs_mean,'r.')
grid
ylabel('Yaw (degs)')
xlabel('Running Time (secs)')
pltfname1(fname1)

figure(8); clf
plot(running_time_secs-running_time_secs(1),interference_flag)
grid
xlabel('Running Time (secs)')
ylabel('Interference Counter')
title(['Interference Counter'])

figure(9); clf
plot(running_time_secs-running_time_secs(1),range_to_source_m)
grid
xlabel('Running Time (secs)')
ylabel('Relative Range (m)')
pltfname1(fname1)


%% Save data
if save_data_flag == 1
save(strcat(start_folder1,'\',vs_filename_base,'_1D_CBF_results_',num2str(ap.integration_time_secs)),...
    'running_time_secs','s_plus_n_to_n_ratio_cbf_db','cbf_1D_angle_est_degs',...
    'ii_low_cbf_level','ii_hi_cbf_level',...
    'HLOS_degs','HLOS_degs_2','roll_degs_mean','pitch_degs_mean',...
    'heading_degs_mean','MARV_east_mean','MARV_north_mean',...
    'source_east','source_north','ap','nav_data','caps_time_start_UTC')
end
