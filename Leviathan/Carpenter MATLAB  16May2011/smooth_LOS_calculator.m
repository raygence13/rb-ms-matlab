function [LOS_smooth, LOS_register_out, time_register_out, time_out,...
    level_smooth, level_register_out] = ...
    smooth_LOS_calculator(ntaps_smooth, LOS, level, time_in, integration_time_secs, ...
    LOS_thresh, LOS_register_in, time_register_in, level_register_in)

%% Make sure input resgister is a column vector
LOS_register_in = LOS_register_in(:);
level_register_in = level_register_in(:);
time_register_in = time_register_in(:);

%% Shift input register by 1
LOS_register_in = circshift(LOS_register_in,1);
level_register_in = circshift(level_register_in,1);
time_register_in = circshift(time_register_in,1);

%% Replace element 1 in input register by current LOS input
LOS_register_in(1) = LOS;
level_register_in(1) = level;
time_register_in(1) = time_in;

%% Find out how many non-NaNs
ii = ~isnan(LOS_register_in);
non_NaN_cntr = sum(ii);
ii1 = find(ii == 1);

%% Interpolate LOS values if LOS_thresh is met, median filter, LPF,
%% else set rate to NaN
if non_NaN_cntr >= LOS_thresh
    LOS_interp = interp1(time_register_in(ii1),LOS_register_in(ii1),...
        time_register_in,'linear','extrap');
    level_interp = interp1(time_register_in(ii1),level_register_in(ii1),...
        time_register_in,'linear','extrap');
    % Median filter twice
    temp1 = medfilt1(LOS_interp,7);
    temp1 = medfilt1(temp1,3);
    % LPF
    nn = floor(ntaps_smooth/3);
    LOS_smooth = filtfilt(ones(nn,1)/nn,1,temp1);
    LOS_smooth = LOS_smooth((ntaps_smooth+1)/2);
    level_smooth = filtfilt(ones(nn,1)/nn,1,level_interp);
    level_smooth = level_smooth((ntaps_smooth+1)/2);
else
    LOS_smooth = NaN;
    level_smooth = NaN;
end

LOS_register_out = LOS_register_in;
level_register_out = level_register_in;
time_register_out = time_register_in;
time_out = time_in - integration_time_secs*(ntaps_smooth+1)/2;

