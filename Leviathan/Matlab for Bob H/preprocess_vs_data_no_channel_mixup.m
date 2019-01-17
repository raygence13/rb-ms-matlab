function [dec_time,fs_cmplx,cmplx_env_sensor_signal, power_per_interval] = ...
    preprocess_vs_data_no_channel_mixup(ap,pressure,ax,ay,az, compute_power_flag)
%
% 19 October 2010

npts = length(pressure);
time = ((0:npts-1)/ap.fs);

if 0
    % 11/20/2009 - The great channel mix-up (kept in code as a reminder of
    % what can happen when care is not taken when connecting sensor to A2D
    % board)
    newx = ax;
    newy = az;
    newz = pressure;
    newp = ay;
    ax = newx;
    ay = newy;
    az = newz;
    pressure = newp;
end

%% NED Coords (see m-file PMS394_vector_DodgePond_31Aug2009_r2
%% & PowerPoint file VS_scale_factor_studies.pptx)
%% X+ out nose, Y+ out Stbd., Z+ Down (see Figure A.5 in NAS
%% Report)
temp = ax;
ax = az;
az = -temp;

%% Apply accelerometer sensitivity (see Wilcoxon spec sheet for VS-205)
ax = ax/ap.accel_sens;
ay = ay/ap.accel_sens;
az = az/ap.accel_sens;

%% Apply data collection gain
pressure = pressure*10^((-ap.data_collection_gain_db)/20);
ax = ax*10^((-ap.data_collection_gain_db)/20);
ay = ay*10^((-ap.data_collection_gain_db)/20);
az = az*10^((-ap.data_collection_gain_db)/20);

%% Apply HYDROPHONE RECEIVE SENSITIVITY
pressure = pressure*10^((-ap.hyd_sens)/20);

%% Subsample data to reduce data rate (just to ease memory limitations of
% Matlab)
time = time(1:ap.ndec:end);
temp = decimate(pressure,ap.ndec,'fir');
pressure = temp;
temp = decimate(ax,ap.ndec,'fir');
ax = temp;
temp = decimate(ay,ap.ndec,'fir');
ay = temp;
temp = decimate(az,ap.ndec,'fir');
az = temp;

%% Compute power in pressure signal in integration period
power_per_interval.whole_band = NaN;
if compute_power_flag == 1
    nsamps_int_period = ap.fs_ndec*ap.integration_time_secs;
    number_of_intervals = floor(length(pressure)/nsamps_int_period);
    power_per_interval.whole_band = pressure.^2;
    power_per_interval.whole_band = ...
        power_per_interval.whole_band(1:nsamps_int_period*number_of_intervals);
    power_per_interval.whole_band = ...
        reshape(power_per_interval.whole_band,nsamps_int_period,number_of_intervals);
    power_per_interval.whole_band = ...
        sum(power_per_interval.whole_band,1)/nsamps_int_period;
end

%% BPF the data
temp = fftfilt(ap.b_bpf,pressure');
pressure = temp;
temp = fftfilt(ap.b_bpf,ax');
ax = temp;
temp = fftfilt(ap.b_bpf,ay');
ay = temp;
temp = fftfilt(ap.b_bpf,az');
az = temp;

%% Compute power in BPF pressure signal in integration period
power_per_interval.bpf = NaN;
if compute_power_flag == 1
    power_per_interval.bpf = pressure.^2;
    power_per_interval.bpf = ...
        power_per_interval.bpf(1:nsamps_int_period*number_of_intervals);
    power_per_interval.bpf = ...
        reshape(power_per_interval.bpf,nsamps_int_period,number_of_intervals);
    power_per_interval.bpf = ...
        sum(power_per_interval.bpf,1)/nsamps_int_period;
end

%% Apply relative gain to pressure to allign with accel channels (see Ben
% Cray conversion sheet)
pressure = pressure*10^((-ap.rho_c_factor-ap.grav_factor)/20);

%% take the derivative of the pressure
p_deriv = differentiate_via_dft(pressure,ap.fs_ndec);

%% put data in form for quad demod function
sensor_signal = [];
sensor_signal(1,:,1) = -p_deriv; % Note negative sign
sensor_signal(1,:,2) = ax;
sensor_signal(1,:,3) = ay;
sensor_signal(1,:,4) = az;

%% Perform quadrature demodulation of BPF vector sensor signals
ndec = floor((ap.fs_ndec/(ap.fhh-ap.fll))/2);
fs_cmplx = ap.fs_ndec/ndec;
[dec_time,cmplx_env_sensor_signal] = quad_demod_VS(sensor_signal,time,ap.fc,1,ndec,60);

