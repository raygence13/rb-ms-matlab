close all
clear
present(0)
%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
debugflag = 0;                  % iterates only one block of data if 1


headersize  = 104;
Fs          = 12500;        	% Sampling Frequency
Nchans      = 384;             	% H + X + Y + Z
Nsens       = Nchans/8;         % number of H only
c           = 1500;            	% medium sound speed [m/s]
rho         = 1026;          	% medium density [kg/m^3]

d           = 0.1406525;        % channel spacing [m]
idx_z       = 1:4:Nchans/2;     % 48 fwd VS elements
idx_y       = idx_z + 1;
idx_x       = idx_y + 1;
idx_h       = idx_x + 1;

Aperture    = d*48*4;           % Aperture length [m]
prop_time   = round(Aperture/c*Fs*10);       % propagation time [bins]
nPts = prop_time;               % block size [bins]
% 
% p = [zeros(Nsens,1)...
%     (0:Nsens-1)'*d...
%     zeros(Nsens,1)].';       % sensor position vector [x y z] 3 x Nsens

p = [zeros(Nsens,1)...
    zeros(Nsens,1)...
    (0:Nsens-1)'*d].';       % sensor position vector [x y z] 3 x Nsens
f0  = 1000:25:6000;          % frequency of interest [Hz]
% f0 = 1050;
l_f0 = length(f0);

tfft  = 2^nextpow2(2*Fs);       % time series fft size
bins  = 1:tfft/2;               % bin number for first half of frequency vector
f0idx = floor(f0.*tfft./Fs);    % frequency of interest index [bin#]
% f0_100nPts = round(Fs/f0*100);  % number of samples for 100 periods of f0

time  = ((1:nPts)*(1/Fs));      % time vector for block data
f     = Fs*(0:tfft/2-1)/tfft;   % frequency vector
overlap = .5;                   % block overlap
OverLapPts = round(nPts*(1-overlap));   % overlap size
h_window = hanning(nPts);       % time series window


% 3-dimensional beamforming parameters
phi = -180:180;      % degrees
l_phi = length(phi);    
theta = -90:5:90;
l_theta = length(theta);

ux = sind(phi).'*cosd(theta);   
ux = ux(:);     % x-component of unit vector
uy = repmat(sind(theta),[1 length(phi)]);
uy = uy(:);     % z-component of unit vector
uz = cosd(phi).'*cosd(theta);
uz = uz(:);     % y-component of unit vector

% ux = repmat(cosd(theta),[1 l_phi]);   
% ux = ux(:);     % x-component of unit vector
% uy = sind(phi).'*sind(theta);
% uy = uy(:);     % z-component of unit vector
% uz = cosd(phi).'*sind(theta);
% uz = uz(:);     % y-component of unit vector

% 3-dimensional unit vector
u = zeros(3,l_phi*l_theta);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);

% 1-dimensional beamforming parameter
ucbf = cosd(phi);

% [fname, pname] = uigetfile('*.rst; ', 'select input RST data file');
[fname, pname] = uigetfile('*.rst', 'select input RST data file');
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


if debugflag
    numBlocks = 1;
else
    numBlocks = floor(timeSampsAvail/nPts);
end

data = fread(fidin,[Nchans,nPts],'float');

H_BFB1 = zeros(l_phi,l_f0);
H_BFB3_Data = zeros(l_phi,l_theta,l_f0);
H_BFB3_Sim = zeros(l_phi,l_theta,l_f0);
VS_BFB = zeros(l_phi,l_theta);

% Fake target
phi_tgt = 0;
theta_tgt = 0;
u_tgt = [sind(phi_tgt)*cosd(theta_tgt);
    sind(theta_tgt);
    cosd(phi_tgt)*cosd(theta_tgt)];
tau = (p.'*u_tgt)./c;

for block = 1:numBlocks
    
    data = data - repmat(mean(data,2),1,nPts); % Remove Mean
    blockdata = bsxfun(@times,data,h_window.');
    
%     for block_ind = 1:blocks_to_avg
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
        for indf = 1:l_f0
            % Data snaps
            H_Data     = H_freq(:,f0idx(indf));            % Hydrophone snapshots at f0

            % Simulated Snaps
            H_Sim = exp(-1i*2*pi*f0(indf)*tau);
            
            % 3-d Hydrophone beamforming
            H_SM        = exp(-1i*2*pi/c*f0(indf)*p.'*u);   % Steering matrix
            H_DS_Data    = H_SM'*H_Data;                    % Delayed & Summed (over sensors) response vs. beams
            H_DS_Sim    = H_SM'*H_Sim;                      % Delayed & Summed (over sensors) response vs. beams
            
            H_BF3_Data       = reshape(H_DS_Data,l_phi,l_theta);
            H_BF3_Sim       = reshape(H_DS_Sim,l_phi,l_theta);
            
            H_SpatialPSD3_Data = H_BF3_Data.*conj(H_BF3_Data); % Spatial PSD
            H_SpatialPSD3_Sim = H_BF3_Sim.*conj(H_BF3_Sim); % Spatial PSD
            
            Normalizer_Data = max(max(H_SpatialPSD3_Data));
            Normalizer_Sim = max(max(H_SpatialPSD3_Sim));
            
            H_BFB3_Data(:,:,indf) = H_SpatialPSD3_Data./Normalizer_Data;
            H_BFB3_Sim(:,:,indf) = H_SpatialPSD3_Sim./Normalizer_Sim;
            
            % 1-d Hydrophone beamforming
            W = exp(-1i*2*pi/c*f0(indf)*d*(0:(Nsens -1))'*ucbf);
            H_BF1 = W'*H_Data;
            H_SpatialPSD1 = (H_BF1.*conj(H_BF1)).';
            Normalizer_Data = max(H_SpatialPSD1);
            H_BFB1(:,indf) = H_SpatialPSD1./Normalizer_Data;
            
        end
        
        % plots for 1D beamforming
            fplot = find(f0==2025);
            figure(1)
            subplot(211)
            imagesc(phi,f0,10*log10(abs(H_BFB1)).')
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Frequency, Hz','FontSize',18);
            title('1-D Beamforming','FontSize',18)
            colorbar
            caxis([-40 0])
            
            subplot(212)
            plot(phi,10*log10(abs(H_BFB1(:,fplot))))
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('dB','FontSize',18)
            title(['FRAZ slice at ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
            xlim([phi(1) phi(end)])
            
            % plots for 3D beamforming
            thetaslice = round(l_theta/2);
            pH_BFB3_fraz = squeeze(H_BFB3_Data(:,thetaslice,:));
            figure(2)
            subplot(211)
            imagesc(phi,f0,10*log10(abs(pH_BFB3_fraz)).')
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Frequency, Hz','FontSize',18);
            title(['Data FRAZ DE slice at $\theta$ = ' num2str(theta(thetaslice)) ' deg'],'FontSize',18)
            colorbar
            caxis([-40 0])

            pH_BFB3_DEAZ = squeeze(H_BFB3_Data(:,:,fplot));
            subplot(212)
            imagesc(phi,theta,10*log10(abs(pH_BFB3_DEAZ)).')
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('DE $\theta$, deg','FontSize',18);
            title(['DE vs AZ at $f_o$ = ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
            colorbar
            caxis([-40 0])
            
            if debugflag
            thetaslice = find(theta==theta_tgt);
            pH_BFB3_fraz = squeeze(H_BFB3_Sim(:,thetaslice,:));
            figure(3)
            subplot(211)
            imagesc(phi,f0,10*log10(abs(pH_BFB3_fraz)).')
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('Frequency, Hz','FontSize',18);
            title(['Simulated FRAZ DE slice at $\theta$ = ' num2str(theta(thetaslice)) ', $\phi$ = ' num2str(phi_tgt) ' deg'],'FontSize',18)
            colorbar
            caxis([-40 0])

            pH_BFB3_DEAZ = squeeze(H_BFB3_Sim(:,:,fplot));
            subplot(212)
            imagesc(phi,theta,10*log10(abs(pH_BFB3_DEAZ)).')
            xlabel('AZ $\phi$, deg','FontSize',18); ylabel('DE $\theta$, deg','FontSize',18);
            title(['DE vs AZ at $f_o$ = ' num2str(f0(fplot)) ' Hz'],'FontSize',18)
            colorbar
            caxis([-40 0])
            end
        % read next chunk of data.
        data  = [data(:,OverLapPts +1: nPts) fread(fidin,...
            [Nchans OverLapPts],'float')];
%     end
    %     keep_for_next_block = data(:,shiftedblockstart_idx:end);
end