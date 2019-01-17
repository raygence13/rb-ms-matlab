function data = compute_LOS_degs(nav_data, ...
    MARV_east, MARV_north, running_time_secs, nav_utc_ref_time, ap, ...
    az_angle_est_degs, el_angle_est_degs, ...
    heading_delta_degs, pitch_delta_degs, roll_delta_degs)

%% Start and stop times for integration period
start_time_utc = running_time_secs + nav_utc_ref_time;
stop_time_utc = start_time_utc + ap.integration_time_secs;
t_middle = (start_time_utc + stop_time_utc)/2;

%% Process attitude data to coincide with integrated time
lat_degs_mean = interp1(nav_data.utc_time,nav_data.lat_degs,t_middle);
lon_degs_mean = interp1(nav_data.utc_time,nav_data.lon_degs,t_middle);
MARV_east_mean = interp1(nav_data.utc_time,MARV_east,t_middle);
MARV_north_mean = interp1(nav_data.utc_time,MARV_north,t_middle);
roll_degs_mean = interp1(nav_data.utc_time,nav_data.roll_degs,t_middle);
pitch_degs_mean = interp1(nav_data.utc_time,nav_data.pitch_degs,t_middle);
heading_degs_mean = interp1(nav_data.utc_time,nav_data.heading_degs,t_middle);

%% Calculate horizontal line of sight (HLOS) from acoustic measurements and
%% UUV navigation data
% [d_cos_N,d_cos_E,d_cos_D] = to_NED(head,pitch,roll,az,el)
% Convert from antenna azimuth and elevation (taking into account platform
% heading, pitch and roll) to direction cosines in a North-East-Down coordinate
% system.  See S. Blackman chap 3., pg. 52
%
% INPUT
%	head - heading of platform in radians ( + is to the east of north)
%	pitch - pitch of platform in radians (+ is up)
%	roll - roll of platform in radians (+ is clockwise looking out through nose of vehicle)
%   az - sonar azimuth angle in radians (+ IS TO STARBOARD)
%   el - sonar elevation angle in radians (-90 DEGS IS STRAIGHT DOWN)
%
% OUTPUT
%	d_cos_N - direction cosine with respect to North axis
%	d_cos_E - direction cosine with respect to East axis
%	d_cos_D - direction cosine with respect to Down axis
% 22 March 2010 - see NAS report pg. 16
% Correct for VS offsets
% 2/2/2011: negative sign change to + (assumes delta angles = serial - nav)
head = (heading_degs_mean + heading_delta_degs)*pi/180;
pitch = (pitch_degs_mean + pitch_delta_degs)*pi/180;
roll = (roll_degs_mean + roll_delta_degs)*pi/180;

eul_vect = [roll pitch head];
nav2body = eulr2dcm(eul_vect);
body2nav = nav2body';
% (use nomenclature fom NAS report, pg. 16)
eta = pi/180*az_angle_est_degs;
epsilon = pi/180*el_angle_est_degs;
u_R = [cos(epsilon)*cos(eta);cos(epsilon)*sin(eta);-sin(epsilon)];
temp = body2nav*u_R;
d_cos_N_CBF = temp(1);
d_cos_E_CBF = temp(2);
d_cos_D_CBF = temp(3);
HLOS_degs = 180/pi*atan2(d_cos_E_CBF,d_cos_N_CBF);
HLOS_degs_2 = HLOS_degs;
ii = find(HLOS_degs < 0);
HLOS_degs(ii) = 360 + HLOS_degs(ii);
VLOS_degs = 180/pi*asin(d_cos_D_CBF);

%% Output data structure
data.lat_degs_mean = lat_degs_mean;
data.lon_degs_mean = lon_degs_mean;
data.MARV_east_mean = MARV_east_mean;
data.MARV_north_mean = MARV_north_mean;
data.roll_degs_mean = roll_degs_mean;
data.pitch_degs_mean = pitch_degs_mean;
data.heading_degs_mean = heading_degs_mean;
data.az_angle_est_degs = az_angle_est_degs;
data.el_angle_est_degs = el_angle_est_degs;
data.HLOS_degs = HLOS_degs;
data.HLOS_degs_2 = HLOS_degs_2;
data.VLOS_degs = VLOS_degs;