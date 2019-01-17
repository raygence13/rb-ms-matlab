function [LOS_smooth, LOS_rate,LOS_register_out,time_register_out,time_out] = ...
    LOS_rate_calculator(LOS, time_in, integration_time_secs, ...
    LOS_thresh, LOS_register_in, time_register_in)

%% Make sure input resgister is a column vector
LOS_register_in = LOS_register_in(:);
time_register_in = time_register_in(:);

%% Shift input register by 1
LOS_register_in = circshift(LOS_register_in,1);
time_register_in = circshift(time_register_in,1);

%% Replace element 1 in input register by current LOS input
LOS_register_in(1) = LOS;
time_register_in(1) = time_in;

%% Find out how many non-NaNs
ii = ~isnan(LOS_register_in);
non_NaN_cntr = sum(ii);
ii1 = find(ii == 1);

%% Interpolate LOS values if LOS_thresh is met, median filter, LPF,
%% and compute derivative, else set rate to NaN
if non_NaN_cntr >= LOS_thresh
    LOS_interp = interp1(time_register_in(ii1),LOS_register_in(ii1),...
        time_register_in,'linear','extrap');
    % Median filter
    temp1 = medfilt1(LOS_interp,3);
    % LPF
    LOS_smooth = filtfilt(ones(3,1)/3,1,temp1);
    % Compute derivative
    h = integration_time_secs;
    LOS_rate = (-LOS_smooth(7) + 8*LOS_smooth(6) - 8*LOS_smooth(4) + LOS_smooth(3))/(12*h);
    LOS_smooth = LOS_smooth((numel(LOS_smooth)+1)/2);
else
    LOS_smooth = NaN;
    LOS_rate = NaN;
end

LOS_register_out = LOS_register_in;
time_register_out = time_register_in;
time_out = time_in - integration_time_secs*(numel(time_register_in)+1)/2;

