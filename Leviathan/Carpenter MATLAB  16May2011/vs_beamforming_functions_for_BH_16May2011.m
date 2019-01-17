% vs_beamforming_functions_for_BH_16May2011
%
% 23 March 2011

clear all
fname1 = mfilename;

commandwindow

%% Run number
run_number = 3;

%% Save data
save_data_flag = 0;

%% Study HLOS rate flag
study_HLOS_rate_flag = 0;

%% Get sound boat GPS and 841 boat gps
load('L:\_MARV_VS_Run_3\GPS_data\CAPS_SOUNDBOAT_GPS_20110323.mat')
load('L:\_MARV_VS_Run_3\Hypack Nav Log 03232011\001_1214.mat')
hypack_utc_time = hypack_841_zulu_hrs_to_UTC(2011,3,23,hypack.zulu_time_hrs);

%% Integration time
integration_time_secs = 2;

%% Compute power in whole band and band of interest (useful for run
%% analysis)
compute_power_flag = 0;

%% Generate TF plot (useful for run analysis)
compute_TF_flag = 1;

%% Generate adaptation parameters, etc.
% EST - eastern standard time (1) or daylight savings time (0)
EST = 0;
ap = vs_adaptation_parameters_r4(EST, integration_time_secs);

%% Get MARV NAV dat
start_folder1 = 'L:\_MARV_VS_Run_3\data\';
nav_folder1 = start_folder1;
%[nav_filename,nav_folder1] = uigetfile('nav*.mat*','NAV Data',start_folder1);
switch run_number
    case 1
        nav_filename = 'nav_data_23-Mar-2011_085704.mat';
    case 2
        nav_filename = 'nav_data_23-Mar-2011_104043.mat';
    case 3
        nav_filename = 'nav_data_23-Mar-2011_114405.mat';
    case 4
        nav_filename = 'nav_data_23-Mar-2011_132153.mat';
    case 5
        nav_filename = 'nav_data_23-Mar-2011_150252.mat';
end
load([nav_folder1 nav_filename])

%% Fix utc_time, etc. (MARV typically has some samples where the UTVC times
%% are flakey)
fix_utc_time_threshold_secs = .5;
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
        vs_filename_base = '20110323085706-0001_caps.dat';
        stop_start = [5 77]; % 77
        source_depth = 4*0.3048;
        gps_nav_offset_secs = 3600; % MARV clocks did not go to daylight savings time on this test
        caps_nav_offset_secs = 0;
        caps_nav_offset_secs = 0;% -19.5
        ap.heading_delta_degs = 0;
    case 2 % straight with moving source
        vs_filename_base = '20110323104045-0001_caps.dat';
        stop_start = [4 56]; % 56
        source_depth = 4*0.3048;
        gps_nav_offset_secs = 3600; % MARV clocks did not go to daylight savings time on this test
        caps_nav_offset_secs = 0;% -19.5
        ap.heading_delta_degs = 0;
    case 3 % straight with moving source
        vs_filename_base = '20110323114406-0001_caps.mat';
        stop_start = [4 83]; % 87
        source_depth = 4*0.3048;
        gps_nav_offset_secs = 3600; % MARV clocks did not go to daylight savings time on this test
        caps_nav_offset_secs = 0;% -19.5
        ap.heading_delta_degs = 0;
    case 4 % straight with moving source
        vs_filename_base = '20110323132155-0001_caps.mat';
        stop_start = [8 77]; %82
        source_depth = 4*0.3048;
        gps_nav_offset_secs = 3600; % MARV clocks did not go to daylight savings time on this test
        caps_nav_offset_secs = 0;% -19.5
        ap.heading_delta_degs = 0;
    case 5 % straight with moving source
        vs_filename_base = '20110323150254-0001_caps.mat';
        stop_start = [5 73]; % 75
        source_depth = 4*0.3048;
        gps_nav_offset_secs = 3600; % MARV clocks did not go to daylight savings time on this test
        caps_nav_offset_secs = 0;% -19.5
        ap.heading_delta_degs = 0;
end

nav_data.utc_time = nav_data.utc_time - caps_nav_offset_secs;
nav_utc_ref_time = nav_utc_ref_time - caps_nav_offset_secs;

%% Set aside memory
rough_number_of_time_steps = ceil((diff(stop_start)+1)*...
    60/ap.integration_time_secs);
bf_conven_db = NaN*zeros(rough_number_of_time_steps,ap.nbeams_phi_1D);
running_time_secs = NaN*zeros(rough_number_of_time_steps,1);
running_NAV_time_secs = NaN*zeros(rough_number_of_time_steps,1);
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
num_beams_close_2_max = NaN*zeros(rough_number_of_time_steps,1);
LOS_rate = NaN*zeros(rough_number_of_time_steps,1);
LOS_smooth = NaN*zeros(rough_number_of_time_steps,1);
smooth_time_out = NaN*zeros(rough_number_of_time_steps,1);
rate_time_out = NaN*zeros(rough_number_of_time_steps,1);
level_smooth = NaN*zeros(rough_number_of_time_steps,1);
level_rate = NaN*zeros(rough_number_of_time_steps,1);

%% Loop over 1 minute (nominal) VS data files
elapsed_time = 0;
for m = 1:stop_start(1)-1
    last_file_num = m;
    if last_file_num < 10
        file_num = ['000',int2str(last_file_num)];
    elseif last_file_num >= 10 && last_file_num < 100
        file_num = ['00',int2str(last_file_num)];
    else
        file_num = ['0',int2str(last_file_num)];
    end
    
    vs_filename = [vs_filename_base(1:15) file_num vs_filename_base(end-8:end-4)];
    disp(vs_filename)
    load([vs_folder vs_filename])
    elapsed_time = elapsed_time + length(pressure.s1)/ap.fs;
    if m == 1
        caps_time_start_UTC_first_file = caps_utc_time_start(vs_filename, ap);
    end
end
caps_time_start_UTC = caps_time_start_UTC_first_file + elapsed_time;
% Ellapsed CAPS time from start of NAV data
relative_caps_utc_time_start = caps_time_start_UTC - nav_utc_ref_time;

kk = 0; % Counter for the total number of integration periods
%%
ntaps_smooth = 27;
LOS_register_in = NaN*ones(ntaps_smooth,1);
time_register_in = NaN*ones(ntaps_smooth,1);
LOS_thresh = (ntaps_smooth+1)/2;
%
level_register_in = NaN*ones(ntaps_smooth,1);
%
LOS_rate_register_in = NaN*ones(5,1);
level_rate_register_in = NaN*ones(5,1);
time_rate_register_in = NaN*ones(5,1);
%%
interference_thresh = 5;
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
    disp(vs_filename)
    
    % Load VS file
    load([vs_folder vs_filename])
    pressure1 = pressure.s1;
    ax1 = ax.s1;
    ay1 = ay.s1;
    az1 = az.s1;
    
    nsamps_int = floor(ap.integration_time_secs*ap.fs);
    file_size_secs(file_cntr) = length(pressure1)/ap.fs;
    nloops = floor(length(pressure1)/nsamps_int);
    lost_time_secs = file_size_secs(file_cntr) - (nloops*nsamps_int/ap.fs);
    %% Loop over integration periods
    for n = 1:nloops
        kk = kk + 1;
        if kk == 1
            running_time_secs(kk) = caps_time_start_UTC;
            %running_NAV_time_secs(kk) =
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
        [dec_time,fs_cmplx,cmplx_env_sensor_signal,power_per_interval, TF] = ...
            preprocess_vs209_data_no_channel_mixup(ap,pressure,ax,ay,az,...
            compute_power_flag, compute_TF_flag);
        if compute_power_flag == 1
            power_per_interval_save.whole_band(kk) = power_per_interval.whole_band;
            power_per_interval_save.bpf(kk) = power_per_interval.bpf;
        end
        if compute_TF_flag == 1
            if kk == 1
                TF_save.Pxx = NaN*zeros(rough_number_of_time_steps,TF.nfft/2+1);
            end
            TF_save.Pxx(kk,:) = 10*log10((TF.Pxx).');
            TF_save.ff = TF.ff;
        end
        
        
        %% Conventional Beamformer
        % Positive theta points up (see NDIA presentation)
        bf_conven = cbf_vs(ap,dec_time,cmplx_env_sensor_signal);
        [imax_row,imax_col] = find(bf_conven == max(max(bf_conven)));
        imax_row = imax_row(1);
        imax_col = imax_col(1);
        max_phi(kk) = ap.phi_steer_deg_2D(imax_col);
        max_theta(kk) = ap.theta_steer_deg_2D(imax_row);
        
        
        if 0
            figure(1); clf
            imagesc(ap.phi_steer_deg_2D,ap.theta_steer_deg_2D,10*log10(bf_conven))
            set(gca,'ydir','rev')
            cb_axes = colorbar;
            set(get(cb_axes,'xlabel'),'string','dB','fontw','bold')
            set(gca,'ydir','normal')
            ylabel('')
            xlabel('')
            title(['CBF Output: ',vs_filename],'interp','none','fontsize',8)
            pltfname1(fname1)
            
            pause(.5)
        end
        
        switch ap.cbf_1D2D_flag
            case 1 % 1D
                bf_conven_db(kk,:) = 10*log10(abs(bf_conven));
            case 2 % 2D
                bf_conven_db(kk,:) = 10*log10(abs(bf_conven(imax_row,:)));
        end
        
        %% Find angle (and level) of max 1D CBF output
        [cbf_1D_angle_est_degs(kk), cbf_level(kk), max_stern_bf_conven_db(kk), ...
            num_beams_close_2_max(kk), s_plus_n_to_n_ratio_cbf_db(kk)] = ...
            oneD_bf_angle_estimate_with_stern_montoring(ap,bf_conven_db(kk,:));
        
        %% compute_LOS_degs
        data = compute_LOS_degs_r2(nav_data, ...
            MARV_east, MARV_north, running_time_secs(kk), ...
            nav_utc_ref_time, ap, cbf_1D_angle_est_degs(kk), 0);
        
        HLOS_degs(kk) = data.HLOS_degs;
        HLOS_degs_2(kk) = data.HLOS_degs_2;
        VLOS_degs(kk) = data.VLOS_degs;
        MARV_east_mean(kk) = data.MARV_east_mean;
        MARV_north_mean(kk) = data.MARV_north_mean;
        roll_degs_mean(kk) = data.roll_degs_mean;
        pitch_degs_mean(kk) = data.pitch_degs_mean;
        heading_degs_mean(kk) = data.heading_degs_mean;
        simple_HLOS_degs(kk) = heading_degs_mean(kk) + cbf_1D_angle_est_degs(kk);
        
        %%  Smooth LOS, Level Calculator (started 7 April 2011)
        if (s_plus_n_to_n_ratio_cbf_db(kk) >= ...
                ap.front_back_threshold_db + ap.LOS_rate_threshold_db) && ...
                (num_beams_close_2_max(kk) <= (numel(ap.phi_steer_deg_1D))-interference_thresh)
            LOS = HLOS_degs(kk);
            s_plus_n_level = cbf_level(kk);
        else
            LOS = NaN;
            s_plus_n_level = NaN;
        end
        
        [LOS_smooth(kk),LOS_register_out,time_register_out,smooth_time_out(kk),...
            level_smooth(kk), level_register_out] = ...
            smooth_LOS_calculator(ntaps_smooth, LOS, s_plus_n_level, running_time_secs(kk), ...
            ap.integration_time_secs, LOS_thresh, LOS_register_in, time_register_in, level_register_in);
        LOS_register_in = LOS_register_out;
        level_register_in = level_register_out;
        time_register_in = time_register_out;
        
        %% LOS, Level Rate Calculator
        [LOS_rate(kk),LOS_rate_register_out,time_rate_register_out,rate_time_out(kk),...
            level_rate(kk), level_rate_register_out] = ...
            LOS_rate_calculator_r2(LOS_smooth(kk), smooth_time_out(kk), ...
            ap.integration_time_secs, LOS_rate_register_in, time_rate_register_in,...
            level_smooth(kk), level_rate_register_in);
        LOS_rate_register_in = LOS_rate_register_out;
        level_rate_register_in = level_rate_register_out;
        time_rate_register_in = time_rate_register_out;
        
    end
end

%% Until timing issues are settled, we have to exclude data where heading
%% rate is excessive
ii_no_excessive_head_rate = find(abs(diff(heading_degs_mean)./diff(running_time_secs)) < 0.5);

%% Set thresholds on Signal Plus Noise to Noise Ratio
%% Also determine where acoustic interference occurs
ii_low_cbf_level = find(s_plus_n_to_n_ratio_cbf_db >= ...
    ap.front_back_threshold_db);
ii_hi_cbf_level = find(s_plus_n_to_n_ratio_cbf_db >= ...
    ap.front_back_threshold_db + ap.LOS_rate_threshold_db);
interference_thresh = 5;
ii_interference = find(num_beams_close_2_max >= (numel(ap.phi_steer_deg_1D))-interference_thresh);
ii_no_interference = find(num_beams_close_2_max <= (numel(ap.phi_steer_deg_1D))-interference_thresh);

%% Sound boat & 841 location during files of interest
min_time = running_time_secs(1);
max_time = max(running_time_secs);
gps_utc_time = gps.utc_time + gps_nav_offset_secs;
ii_good_gps_time = find(gps_utc_time <= max_time & gps_utc_time >= min_time);
[source_east,source_north] = mfwdtran(mstruct,gps.lat(ii_good_gps_time),gps.lon(ii_good_gps_time));
source_east = interp1(gps_utc_time(ii_good_gps_time),source_east,...
    linspace(min_time,max_time,numel(running_time_secs)));
source_north = interp1(gps_utc_time(ii_good_gps_time),source_north,...
    linspace(min_time,max_time,numel(running_time_secs)));

ii_good_hypack_time = find(hypack_utc_time <= max_time & hypack_utc_time >= min_time);

%% Compute source speed
source_speed_mpers = sqrt(diff(source_east).^2 + diff(source_north).^2)./diff(running_time_secs');

%% Compute geometric HLOS and bearing using GPS
HLOS_geometric_degs = 180/pi*atan2(-(MARV_east_mean-source_east'),-(MARV_north_mean-source_north'));
ii2 = find(HLOS_geometric_degs < 0);
HLOS_geometric_degs(ii2) = HLOS_geometric_degs(ii2) + 360;
HLOS_geometric_degs_2 = HLOS_geometric_degs;
%HLOS_geometric_rate = diff(HLOS_geometric_degs)./diff(running_time_secs);
HLOS_geometric_rate = filter([-1 8 0 -8 1],1,HLOS_geometric_degs)/(12*ap.integration_time_secs);

%% Compute range to source from GPS
range_to_source_m = sqrt((MARV_east_mean-source_east').^2 + (MARV_north_mean-source_north').^2);

%% Image CBF Results
if 1
    temp1 = ((bf_conven_db)'-max(max(bf_conven_db)))';
    ii_nonan = find(~isnan(running_time_secs));
    temp1 = temp1(ii_nonan,:);
    ii_nuke = find(temp1(:,180) >= max(max(temp1)) - 15);
    temp1(ii_nuke,:) = NaN;
    temp1 = temp1 - nanmax(nanmax(temp1));
    figure(1); clf
    imagesc((ap.phi_steer_deg_1D),running_time_secs(ii_nonan)-running_time_secs(1),...
        temp1)
    set(gca,'ydir','rev')
    cb_axes = colorbar;
    set(get(cb_axes,'xlabel'),'string','dB','fontw','bold')
    set(gca,'ydir','normal')
    ylabel('Time (secs)')
    xlabel('Azimuth Angle \eta (degs)')
    title(['CBF Output (Horizontal Plane): ',vs_filename],'interp','none','fontsize',8)
    pltfname1(fname1)
end

%% TF plot
if compute_TF_flag == 1
    iinan = find(~isnan(running_time_secs));
    figure(2); clf
    imagesc(TF_save.ff,running_time_secs(iinan)-running_time_secs(1),TF_save.Pxx(iinan,:));%
    xlabel('Frequency (Hz)')
    ylabel('Running Time (secs)')
    colorbar
    set(gca,'ydir','rev')
    title([num2str(ap.fll),' - ',num2str(ap.fhh),' Hz'])
    pltfname1(fname1)
end


figure(3); clf
plot(running_time_secs-running_time_secs(1),cbf_level-max(max(bf_conven_db)),...
    running_time_secs-running_time_secs(1),max_stern_bf_conven_db-max(max(bf_conven_db)),'r')
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
    running_time_secs(ii_interference)-running_time_secs(1),...
    cbf_1D_angle_est_degs(ii_interference),'go',...
    running_time_secs-running_time_secs(1),heading_degs_mean,'c.')
hl1 = legend('CBF (1D)','Detection Occurs','Interference','Mean(heading)');
set(hl1,'fontsize',7)
set(gca,'ydir','rev')
ylabel('Bearing Estimate (degs)')
grid
title(['Integration Time (secs) = ',num2str(ap.integration_time_secs)])
subplot(212)
plot(running_time_secs-running_time_secs(1),HLOS_degs,'.',...
    running_time_secs(ii_low_cbf_level)-running_time_secs(1),...
    HLOS_degs(ii_low_cbf_level),'r.',...
    running_time_secs(ii_interference)-running_time_secs(1),...
    HLOS_degs(ii_interference),'go',...
    running_time_secs-running_time_secs(1),HLOS_geometric_degs,'c')
grid
hl3 = legend('CBF','Det. Occurs','Interference','Ground Truth from GPS');
set(hl3,'fontsize',7)
set(gca,'ydir','rev')
xlabel('Running Time (secs)')
ylabel('HLOS Estimate (degs)')
title(['front_back_threshold_db = ',num2str(ap.front_back_threshold_db),...
    ' LOS_rate_threshold_db = ',num2str(ap.LOS_rate_threshold_db)],...
    'interp','none')
pltfname1(fname1)

figure(5); clf
plot(running_time_secs(ii_low_cbf_level)-running_time_secs(1),...
    HLOS_degs(ii_low_cbf_level),'r.',...
    running_time_secs(ii_interference)-running_time_secs(1),...
    HLOS_degs(ii_interference),'go',...
    running_time_secs-running_time_secs(1),HLOS_geometric_degs,'c',...
    smooth_time_out-running_time_secs(1),LOS_smooth,'bx')
grid
hl3 = legend('Det. Occurs','Interference','Ground Truth from GPS');
set(hl3,'fontsize',7)
set(gca,'ydir','rev')
xlabel('Running Time (secs)')
ylabel('HLOS Estimate (degs)')
title(['front_back_threshold_db = ',num2str(ap.front_back_threshold_db),...
    ' LOS_rate_threshold_db = ',num2str(ap.LOS_rate_threshold_db)],...
    'interp','none')
pltfname1(fname1)

figure(6); clf
plot(MARV_east_mean, MARV_north_mean,'k',...
    MARV_east_mean(ii_low_cbf_level), MARV_north_mean(ii_low_cbf_level),'b.',...
    source_east,source_north,'r.')
for kkk = 1:100:length(running_time_secs)
    text(MARV_east_mean(kkk), MARV_north_mean(kkk),...
        [' ',num2str(running_time_secs(kkk)-running_time_secs(1),'% 6.0f')])
end
axis equal
hl2 = legend('MARV','Detection Occurred','Source');
set(hl2,'fontsize',7)
xlabel('Easting (m)')
ylabel('Northing(m)')
grid
pltfname1(fname1)

figure(7); clf
plot(running_time_secs-running_time_secs(1),range_to_source_m,...
    running_time_secs(ii_low_cbf_level)-running_time_secs(1),...
    range_to_source_m(ii_low_cbf_level),'r.',...
    'linew',2)
grid
legend('','Detection Occurs')
xlabel('Running Time (secs)')
ylabel('Range (m)')
pltfname1(fname1)

if compute_power_flag == 1
    figure(8); clf
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
end

figure(9); clf
subplot(311)
plot(nav_data.utc_time-caps_time_start_UTC,nav_data.roll_degs,...
    running_time_secs-running_time_secs(1),roll_degs_mean,'r.')
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

%% Plot smooth level
figure(10); clf
plot(running_time_secs-running_time_secs(1),cbf_level,'.',...
    smooth_time_out-running_time_secs(1),level_smooth,'r.')
grid
hl5 = legend('level','Smooth');
set(hl5,'fontsize',7,'interp','none')
set(gca,'ydir','rev')
xlabel('Running Time (secs)')
ylabel('dB')
pltfname1(fname1)

%% Plot level rate
figure(11); clf
plot(rate_time_out(1:end)-running_time_secs(1),level_rate,'.',...
    'linew',2)
grid
xlabel('Running Time (secs)')
ylabel('Level Rate')
pltfname1(fname1)

%% Plot smooth LOS
figure(12); clf
plot(running_time_secs-running_time_secs(1),HLOS_degs,'.',...
    running_time_secs-running_time_secs(1),HLOS_geometric_degs,...
    smooth_time_out-running_time_secs(1),LOS_smooth,'r.')
grid
hl5 = legend('HLOS','HLOS_geometric_degs','Smooth');
set(hl5,'fontsize',7,'interp','none')
set(gca,'ydir','rev')
xlabel('Running Time (secs)')
ylabel('HLOS (degs)')
pltfname1(fname1)

%% Plot LOS rate
LOS_rate_2 = filter([-1 8 0 -8 1],1,LOS_smooth)/(12*ap.integration_time_secs);
figure(13); clf
plot(running_time_secs(1:end)-running_time_secs(1),HLOS_geometric_rate,...
    running_time_secs(1:end)-running_time_secs(1),LOS_rate_2,'yx',...
    rate_time_out(1:end)-running_time_secs(1),LOS_rate,'gx',...
    'linew',2)
grid
legend('Ground Truth','Estimated')
xlabel('Running Time (secs)')
ylabel('HLOS Rate (degs/sec)')
pltfname1(fname1)

if study_HLOS_rate_flag == 1
    %%  Study HLOS rate
    %
    figure(3); hold on
    plot(running_time_secs(ii_interference)-running_time_secs(1),...
        cbf_level(ii_interference)-max(max(bf_conven_db)),'g.')
    legend('Max(CBF)','Stern (CBF)','Interference')
    
    ii_limit_angle = find(abs(cbf_1D_angle_est_degs) <= 150);
    ii_hi_level_no_interference = intersect(ii_hi_cbf_level, ii_no_interference);
    ii_hi_level_no_interference = intersect(ii_hi_level_no_interference,ii_limit_angle);
    ii_hi_level_no_interference = intersect(ii_hi_level_no_interference,ii_no_excessive_head_rate);
    HLOS_degs_candidates = HLOS_degs(ii_hi_level_no_interference);
    running_time_secs_candidates = running_time_secs(ii_hi_level_no_interference);
    
    candidates = zeros(1,numel(HLOS_degs));
    candidates(ii_hi_level_no_interference) = 1;
    good_run = filter(ones(1,7),1,candidates);
    good_run_2 = zeros(size(good_run));
    ii = find(good_run >= 2);
    good_run_2(ii) = 1;
    good_run_3 = medfilt1(good_run_2,3);
    
    figure(12); clf
    subplot(211)
    plot(running_time_secs-running_time_secs(1),HLOS_degs)
    grid
    subplot(212)
    plot(running_time_secs-running_time_secs(1),candidates,'o',...
        running_time_secs-running_time_secs(1),good_run_3,'r.')
    grid
    
    jj = find(diff(good_run_3) ~= 0);
    n_legs = floor(numel(jj)/2);
    jj = jj(1:2*n_legs);
    jj = reshape(jj,2,n_legs);
    
    ntaps_mdfilt_1 = 15; % 15
    ntaps_medfilt_2 = round(ntaps_mdfilt_1/3);
    ntaps_lpf = 21; % 21
    for kkk = 1:n_legs
        mm = find(ii_hi_level_no_interference >= jj(1,kkk) & ...
            ii_hi_level_no_interference <= jj(2,kkk));
        mm1 = ii_hi_level_no_interference(mm);
        temp_data = HLOS_degs(mm1);
        numel(temp_data)
        if numel(temp_data) >= 4
            temp_time = running_time_secs(mm1);
            t1 = find(running_time_secs >= temp_time(1));
            t1 = t1(1);
            t2 = find(running_time_secs <= temp_time(end));
            t2 = t2(end);
            time = running_time_secs(t1:t2);
            %time = running_time_secs(jj(1,kkk):jj(2,kkk));
            data = interp1(temp_time,temp_data,time);
            %data = HLOS_degs(jj(1,kkk):jj(2,kkk));
            numel(data)
            if numel(data) >= 3*ntaps_lpf
                HLOS_medfilt = medfilt1(data,ntaps_mdfilt_1);
                HLOS_medfilt = medfilt1(HLOS_medfilt,ntaps_medfilt_2);
                HLOS_lowpass_interp = filtfilt(ones(ntaps_lpf,1)/ntaps_lpf,1,HLOS_medfilt);
                
                %% Compute LOS rate via filtering
                HLOS_rate_2 = filter([-1 8 0 -8 1],1,HLOS_lowpass_interp)/(12*ap.integration_time_secs);
                %HLOS_rate_2 = filter([8 -1 0 1 -8],1,HLOS_lowpass_interp)/(12*ap.integration_time_secs);
                %HLOS_rate_2 = filter([1 0 -1],1,HLOS_lowpass_interp)/(ap.integration_time_secs);
                HLOS_rate_2(1:4) = [];
                los_rate_time = time;
                los_rate_time(1:4) = [];
                %
                LOS_rate_2 = filter([-1 8 0 -8 1],1,LOS_smooth)/(12*ap.integration_time_secs);
                
                figure(13); clf
                plot(running_time_secs-running_time_secs(1),HLOS_degs,'.',...
                    running_time_secs_candidates-running_time_secs(1),HLOS_degs_candidates,'ro',...
                    time-running_time_secs(1),data,'cx',...
                    time-running_time_secs(1),HLOS_medfilt,'gx',...
                    time-running_time_secs(1),HLOS_lowpass_interp,'y.',...
                    running_time_secs-running_time_secs(1),HLOS_geometric_degs,...
                    smooth_time_out-running_time_secs(1),LOS_smooth,'k.')
                grid
                hl5 = legend('HLOS','HLOS candidates','Interp',...
                    'HLOS_medfilt','HLOS_lpf','HLOS_geometric_degs',...
                    'LOS smooth');
                set(hl5,'fontsize',7,'interp','none')
                set(gca,'ydir','rev')
                xlabel('Running Time (secs)')
                ylabel('HLOS (degs)')
                pltfname1(fname1)
                
                figure(14); clf
                plot(running_time_secs(1:end)-running_time_secs(1),HLOS_geometric_rate,...
                    los_rate_time-running_time_secs(1),HLOS_rate_2,'r',...
                    running_time_secs(1:end)-running_time_secs(1),LOS_rate_2,'yx',...
                    rate_time_out(1:end)-running_time_secs(1),LOS_rate,'gx',...
                    'linew',2)
                grid
                legend('Ground Truth','Estimated')
                xlabel('Running Time (secs)')
                ylabel('HLOS Rate (degs/sec)')
                pltfname1(fname1)
                
                figure(15); clf
                plt_index = t1+4:t2;
                plot(los_rate_time-running_time_secs(1),range_to_source_m(plt_index),...
                    'linew',2)
                grid
                xlabel('Running Time (secs)')
                ylabel('Range (m)')
                pltfname1(fname1)
                
                %% Plot for NDIA
                figure(16); clf
                plot(time-running_time_secs(1),data,'.',...
                    time-running_time_secs(1),HLOS_lowpass_interp,'g.',...
                    running_time_secs(plt_index)-running_time_secs(1),...
                    HLOS_geometric_degs(plt_index),'r')
                grid
                hl5 = legend('HLOS estimates',...
                    'Smoothed Estimates','Ground Truth');
                set(hl5,'fontsize',7,'interp','none')
                set(gca,'ydir','rev')
                xlabel('Running Time (secs)')
                ylabel('HLOS (degs)')
                pltfname1(fname1)
                
                figure(17); clf
                plot(running_time_secs(plt_index)-running_time_secs(1),...
                    HLOS_geometric_rate(plt_index),...
                    los_rate_time-running_time_secs(1),HLOS_rate_2,'r',...
                    'linew',2)
                grid
                legend('Ground Truth','Estimated')
                xlabel('Running Time (secs)')
                ylabel('HLOS Rate (degs/sec)')
                pltfname1(fname1)
                
                figure(18); clf
                plt_index = t1+4:t2;
                plot(MARV_east_mean(plt_index), MARV_north_mean(plt_index),'.',...
                    source_east(plt_index),source_north(plt_index),'r.')
                for kkk = 1:40:length(plt_index)
                    text(MARV_east_mean(plt_index(kkk)), MARV_north_mean(plt_index(kkk)),...
                        [' ',num2str(running_time_secs(plt_index(kkk))-running_time_secs(1),'% 6.0f')])
                    text(source_east(plt_index(kkk)),source_north(plt_index(kkk)),...
                        [' ',num2str(running_time_secs(plt_index(kkk))-running_time_secs(1),'% 6.0f')])
                end
                hold on
                plot([MARV_east_mean(plt_index(1)) source_east(plt_index(1))],...
                    [MARV_north_mean(plt_index(1)) source_north(plt_index(1))],'c',...
                    'linew',2)
                plot([MARV_east_mean(plt_index(end)) source_east(plt_index(end))],...
                    [MARV_north_mean(plt_index(end)) source_north(plt_index(end))],'c',...
                    'linew',2)
                axis equal
                hl2 = legend('MARV','Source');
                set(hl2,'fontsize',7)
                xlabel('Easting (m)')
                ylabel('Northing(m)')
                grid
                pltfname1(fname1)
                
                figure(19); clf
                plot(los_rate_time-running_time_secs(1),HLOS_rate_2,...
                    running_time_secs(plt_index)-running_time_secs(1),HLOS_geometric_rate(plt_index),'r',...
                    'linew',2)
                grid
                legend('Estimated','Ground Truth')
                xlabel('Running Time (secs)')
                ylabel('HLOS Rate (degs/sec)')
                pltfname1(fname1)
                
                
                
                pause
            end
        end
    end
end

%% Save data
if save_data_flag == 1
    save(strcat(start_folder1,'\',vs_filename_base(1:end-4),'_1D_CBF_results_',num2str(ap.integration_time_secs)),...
        'running_time_secs','s_plus_n_to_n_ratio_cbf_db','cbf_1D_angle_est_degs',...
        'ii_low_cbf_level','ii_hi_cbf_level','ii_interference',...
        'HLOS_degs','HLOS_degs_2','roll_degs_mean','pitch_degs_mean',...
        'heading_degs_mean','MARV_east_mean','MARV_north_mean',...
        'source_east','source_north','ap','nav_data','caps_time_start_UTC',...
        'range_to_source_m','fname1')
end
