function [LOS_rate,LOS_register_out,time_register_out,time_out,...
    level_rate, level_register_out] = ...
    LOS_rate_calculator_r2(LOS, time_in, integration_time_secs, ...
    LOS_register_in, time_register_in, level, level_register_in)

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

%% Compute derivative
% HLOS_rate_2 = filter([-1 8 0 -8 1],1,HLOS_lowpass_interp)/(12*ap.integration_time_secs)
LOS_rate = (-LOS_register_in(1) + 8*LOS_register_in(2) - ...
    8*LOS_register_in(4) + LOS_register_in(5))/(12*integration_time_secs);
LOS_register_out = LOS_register_in;
%
%level_rate = (-level_register_in(1) + 8*level_register_in(2) - ...
    %8*level_register_in(4) + level_register_in(5))/(12*integration_time_secs);
p = polyfit(time_register_in,level_register_in,1);
level_rate = p(1);
level_register_out = level_register_in;

time_register_out = time_register_in;
time_out = time_in - integration_time_secs*(numel(time_register_out)+1)/2;

