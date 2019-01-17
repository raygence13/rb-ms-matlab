function [angle_est_degs, level, max_stern_bf_db, ...
    num_beams_close_2_max, s_plus_n_to_n_ratio_db] = ...
    oneD_bf_angle_estimate_with_stern_montoring(ap,bf_1D_db)

%% Find Conventional BF response towards stern
nc = numel(ap.phi_steer_deg_1D);
ang1 = ap.stern_angle_limits_degs(1);
ang2 = ap.stern_angle_limits_degs(2);
ii_stern = find(ap.phi_steer_deg_1D <= ang1 | ap.phi_steer_deg_1D >= ang2);
stern_bf = 10.^(bf_1D_db(:,ii_stern)/20);
max_stern_bf_db = (20*log10(max(abs(stern_bf)')))';

%% Interference Flag (detect acomms or tracking pinger interference)
max_beam_response_db = max(bf_1D_db);
ii = find(bf_1D_db >= max_beam_response_db - 3);
num_beams_close_2_max = numel(ii);

%% Without beam interpolation
if ap.beam_interp_flag == 0
    %% Find angle of max BF output (and level)
    angle_est_degs = zeros(size(bf_1D_db,1),1);
    level = zeros(size(bf_1D_db,1),1);
    for k = 1: size(bf_1D_db,1)
        ii = find(bf_1D_db(k,:) == max(bf_1D_db(k,:)));
        angle_est_degs(k) = ap.phi_steer_deg_1D(ii(1));
        level(k) = bf_1D_db(k,ii(1));
    end
end

%% Find angle of max BF output via beam interpolation (and level)
if ap.beam_interp_flag == 1
    angle_est_degs = zeros(size(bf_1D_db,1),1);
    level = zeros(size(bf_1D_db,1),1);
    steer_angles_sine_space = sin(pi/180*ap.phi_steer_deg_1D);
    bf_1D = 10.^(bf_1D_db/20);
    for k = 1: size(bf_1D_db,1)
        % Find location of peak value and 2 adjacent bins
        ii = find(bf_1D_db(k,:) == max(bf_1D_db(k,:)));
        ii = ii(1);
        aa = ap.phi_steer_deg_1D(ii);
        if ii == 1
            data = bf_1D(k,1:3);
            angles = steer_angles_sine_space(1:3);
        elseif ii == nc
            data = bf_1D(k,nc-2:nc);
            angles = steer_angles_sine_space(nc-2:nc);
        else
            data = bf_1D(k,ii-1:ii+1);
            angles = steer_angles_sine_space(ii-1:ii+1);
        end
        
        % Least squares fit
        H = zeros(3,3);
        for kk = 1:3
            H(kk,:) = [1 angles(kk) angles(kk)*angles(kk)];
        end
        % Invert H'*H
        temp = invert_3by3(H'*H); % Hard wired 
        %temp = inv(H'*H);  % Matlab function
        angle_est = temp*H'*data';
        angle_est_degs(k) = real(180/pi*asin(complex(-angle_est(2)/(2*angle_est(3)))));
        if aa > 90
            angle_est_degs(k) = 180 - angle_est_degs(k);
        end
        if aa < -90
            angle_est_degs(k) = -(180 + angle_est_degs(k));
        end
        level(k) = bf_1D_db(k,ii);
    end
end

s_plus_n_to_n_ratio_db = level - max_stern_bf_db;
