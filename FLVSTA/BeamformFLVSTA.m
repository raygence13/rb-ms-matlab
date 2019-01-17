close all
clear
%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headersize  = 104;
Fs          = 12500;        	% Sampling Frequency
Nchans      = 384;             	% H + X + Y + Z
c           = 1500;            	% medium sound speed [m/s]
rho         = 1026;          	% medium density [kg/m^3]

d           = 0.1406525;        % channel spacing [m]
idx_z       = 1:4:Nchans/2;     % 48 fwd VS elements
idx_y       = idx_z + 1;
idx_x       = idx_y + 1;
idx_h       = idx_x + 1;

Aperture    = d*48*4;           % Aperture length [m]
prop_time   = round(Aperture/c*Fs*10);       % propagation time [bins]
nPts = prop_time;              % block size

p = [zeros(Nchans/8,1)...
    (0:(Nchans/8)-1)'*d...
    zeros(Nchans/8,1)].';       % sensor position vector [x y z] 3 x Nchans/8

% f0  = 1000:100:2025;                     % frequency of interest [Hz]
f0 = 1000;
debugflag = 1;

tfft  = 2^nextpow2(2*Fs);       % time series fft size
bins  = 1:tfft/2;               % bin number for first half of frequency vector
f0idx = floor(f0.*tfft./Fs);      % frequency of interest index [bin#]
% f0_100nPts = round(Fs/f0*100);  % number of samples for 100 periods of f0

time  = ((1:nPts)*(1/Fs));      % time vector for block data
f     = Fs*(0:tfft/2-1)/tfft;   % frequency vector
overlap = .5;                   % block overlap
OverLapPts = round(nPts*(1-overlap));   % overlap size
h_window = hanning(nPts);       % time series window


% 3-dimensional beamforming
phi = -180:.5:180;      % degrees
l_phi = length(phi);    
theta = -10:10;
l_theta = length(theta);

ux = sind(phi).'*cosd(theta);   
ux = ux(:);     % x-component of unit vector
uy = cosd(phi).'*cosd(theta);
uy = uy(:);     % y-component of unit vector
uz = repmat(sind(theta),[1 length(phi)]);
uz = uz(:);     % z-component of unit vector

u = ones(3,l_phi*l_theta);      % unit vector
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);
kp = 2*pi/c*f0*p.'*u;           % steered wave-number

% 1-dimensional beamforming
ucbf = cosd(phi);
W = exp(-1i*2*pi/c*f0*d*(0:(Nchans/8-1))'*ucbf);

[fname, pname] = uigetfile('*.rst; ', 'select input RST data file');
if(~ischar(fname))
    return;
end
[nullpath filename theext]=fileparts(fname);
cd(pname);
input_fullname = [pname fname];

fidin    = fopen(input_fullname,'r','l');

fseek(fidin,0,'eof');
filebytes       = ftell(fidin);
timeSampsAvail  = floor((filebytes - headersize) / (4 * Nchans) );
fseek(fidin,0,'bof');
fseek(fidin,headersize,'cof');

numBlocks = floor(timeSampsAvail/nPts);

data = fread(fidin,[Nchans,nPts],'float');
H_BFB1 = zeros(1,l_phi,4);
H_BFB3 = zeros(l_phi,l_theta,4);
VS_BFB = zeros(l_phi,l_theta,4);

if debugflag
    numBlocks = 1;
end

for block = 1:numBlocks
    
    data = data - repmat(mean(data,2),1,nPts); % Remove Mean
    blockdata = bsxfun(@times,data,h_window.');
    
    for blockavg = 1:4
        fftdata = fft(blockdata,tfft,2)./tfft;
        fftdata = fftdata(:,bins);
        PSDData = fftdata.*conj(fftdata).*2;
        
        %  Temporal Frequency
        h  = blockdata(idx_h,:);        % Hydrophone time series [Pa (s)]
        H_freq = fft(h,tfft,2)./tfft;   % Hydrophone spectrum at each sensor
        H_freq = H_freq(:,bins);        % Hydrophone positive spectrum at each sensor
        H_PSD = H_freq.*conj(H_freq);   % Hydrophone PSD at each sensor
        H_avgPSD = mean(H_PSD,1);
        
        % beamforming data
        
        % 3-d Hydrophone beamforming
        H_Snaps     = (H_freq(:,f0idx));    % Hydrophone snapshots at f0
        H_Steer     = exp(-1i*kp);          % Steering matrix
        H_Summed    = (H_Steer'*H_Snaps);   % Summed (over sensors) response vs. beams
        H_BF3 = reshape(H_Summed,l_phi,l_theta);
        H_SpatialPSD3 = H_BF3.*conj(H_BF3); % Spatial PSD
        H_BFB3(:,:,blockavg) = H_SpatialPSD3./max(max(H_SpatialPSD3));
        
        % 1-d Hydrophone beamforming
        H_BF1 = ((H_Snaps)'*W);
        H_SpatialPSD1 = H_BF1.*conj(H_BF1);
        H_BFB1(:,:,blockavg) = H_SpatialPSD1./max(max(H_SpatialPSD1));
        
        % XYZ beamforming
        x_a  = blockdata(idx_x,:)*9.81; % x channels acceleration over time [m/s^2 (s)]
        X_a = fft(x_a,tfft,2)./tfft;    % X Acceleration spectrum
        X_p = X_a.*rho.*c./(1i*2*pi*f0);% X Pressure spectrum
        X_p = X_p(:,bins);
        X_snaps = X_p(:,f0idx);
        Xweighted = X_snaps*ux';
        
        y_a  = blockdata(idx_y,:)*9.81; % y channels acceleration over time [m/s^2 (s)]
        Y_a = fft(y_a,tfft,2)./tfft;    % Y Acceleration spectrum
        Y_p = Y_a.*rho.*c./(1i*2*pi*f0);% Y Pressure spectrum
        Y_p = Y_p(:,bins);
        Y_snaps = Y_p(:,f0idx);
        Yweighted = Y_snaps*uy';
        
        z_a  = blockdata(idx_z,:)*9.81; % z channels acceleration over time [m/s^2 (s)]
        Z_a = fft(z_a,tfft,2)./tfft;    % Z Acceleration spectrum
        Z_p = Z_a.*rho.*c./(1i*2*pi*f0);% Z Pressure spectrum
        Z_p = Z_p(:,bins);
        Z_snaps = Z_p(:,f0idx);
        Zweighted = Z_snaps*uz';
        
        XYZweighted = Xweighted + Yweighted + Zweighted;
        alpha = 0.5;

        VS_Delayed = exp(-1i*(alpha*H_Steer+(1-alpha)*XYZweighted));
        VS_Sum = sum(VS_Delayed,1);
        VS_Beamformed = reshape(VS_Sum,l_phi,l_theta);
        VS_SpatialPSD = VS_Beamformed.*conj(VS_Beamformed);
        VS_BFB(:,:,blockavg) = VS_SpatialPSD./max(max(VS_SpatialPSD));
        
        if blockavg == 4
            H_BF1avg = mean(H_BFB1,3);
            H_BF3avg = mean(H_BFB3,3);
            VS_BFavg = mean(VS_BFB,3);
            
            figure(2)
            subplot(221)
            imagesc(phi,theta,10*log10(H_BF3avg).');
            colorbar
            xlabel('\phi'); ylabel('\theta')
            title(['Hydrophone Only 3-d Beam Response, f_0 = ' num2str(f0) ' Hz'])
%             pause(0.5)
            
            subplot(223)
            imagesc(phi,theta,10*log10(VS_BFavg).');
            xlabel('\phi'); ylabel('\theta')
            title('VS Response')
            colorbar
            
            subplot(2,2,[2 4])
            thetaslice = find(theta==10);
            plot(phi,10*log10(abs(H_BF3avg(:,thetaslice))))
            grid on
            xlabel('\phi'); ylabel('dB')
            title(['Hydrophone Only Beam Response DE slice, \theta = ' num2str(theta(thetaslice)) ])
            
            figure(3)
            plot(phi,10*log10(abs(H_BF1avg)))
            grid on
            xlabel('\phi, deg'); ylabel('dB')
            title(['Hydrophone Only 1-d Beam Response, f_0 = ' num2str(f0) ' Hz'])
%             pause
        end
        %read next chunk of data.
        data  = [data(:,OverLapPts +1: nPts) fread(fidin,...
            [Nchans OverLapPts],'float')];
    end
    %     keep_for_next_block = data(:,shiftedblockstart_idx:end);
end