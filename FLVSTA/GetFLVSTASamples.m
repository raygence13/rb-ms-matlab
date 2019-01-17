function [ VSData, RollData, sampsLeft ] = GetFLVSTASamples( fid, NumSamples )
%GetFLVSTASamples
% Should be 96 Hydrophone elements and 96 VS Elements
% Inputs
% fid: obtained from fopen
% NumSamples: Number of time samples to read

% Outputs
% VSData: [H;X;Y;Z]
% H: Hydrophone data, [Pa](t)
% X,Y,Z: Accelerometer data, [m/s^2](t)
% RollData: [XYZ]1, [XYZ]2, ..., [XYZ]96
% RollData in units of Deg(t)
% Gains: [Hfwd,Haft,Afwd,Aaft]
% Gains in units of [V/Pa](t) for hydrophone, [V/g](t) for accelerometer


%%
%Determine how many samples are in the file
%
    curPos= ftell(fid);                 % Get current Position
    fseek(fid, 0, 'eof');               % Seek to End of File
    bytesRemain= ftell(fid) - curPos;   % Calculate bytes remaining
    sampsLeft= bytesRemain/(4 * 484);   % Bytes per Samples
    fseek(fid, curPos, 'bof');
    if sampsLeft > NumSamples
        SamplesToRead= NumSamples;
    else
        SamplesToRead= sampsLeft;
    end

    %% 
    % If there are samples read 
    %
    if SamplesToRead > 0
        %% Read a chunk of Raw Data
        rawData= fread(fid, [484 SamplesToRead], 'float');
    
        %% Pull Out The Gains Hydro Fwd, Hydro Aft, Accel Fwd, Accel Aft
        Gains= 10.^(rawData(481:484,:)./10);

        %% Roll Data VS1,VS2.....VS96
        RollData= rawData(385:480,:);
        VS_Cal = 0.64;          % Vector Sensor Calibration Gain V/g
%         H_Cal = 10^(-190/10)*1e6;   % Hydrophone Calibration Gain V/Pa
        %%VSData H,X,Y,Z
        Data  = rawData(1:384,:);
%         H_V     = Data(1:4:end,:); % using nominal calibration values
%         X_V     = Data(2:4:end,:)/VS_Cal;
%         Y_V     = Data(3:4:end,:)/VS_Cal;
%         Z_V     = Data(4:4:end,:)/VS_Cal;
        H_V     = Data(1:4:end,:);  % will use measured values outside
        X_V     = Data(2:4:end,:);
        Y_V     = Data(3:4:end,:);
        Z_V     = Data(4:4:end,:);
        
        Hfwd = bsxfun(@times,H_V(1:48,:),1./Gains(1,:));
        Haft = bsxfun(@times,H_V(49:end,:),1./Gains(2,:));
        H = [Hfwd;Haft];
        
        Xfwd = bsxfun(@times,X_V(1:48,:),1./Gains(3,:));
        Xaft = bsxfun(@times,X_V(49:end,:),1./Gains(4,:));
        X = [Xfwd;Xaft]*9.81;
        
        Yfwd = bsxfun(@times,Y_V(1:48,:),1./Gains(3,:));
        Yaft = bsxfun(@times,Y_V(49:end,:),1./Gains(4,:));
        Y = [Yfwd;Yaft]*9.81;
        
        Zfwd = bsxfun(@times,Z_V(1:48,:),1./Gains(3,:));
        Zaft = bsxfun(@times,Z_V(49:end,:),1./Gains(4,:));
        Z = [Zfwd;Zaft]*9.81;
        
        VSData = [H;X;Y;Z];
        
    else
        fprintf('ERROR: No Samples Left To Read\n');
    end
end