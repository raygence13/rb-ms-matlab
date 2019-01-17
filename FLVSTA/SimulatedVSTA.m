% VSTA Beampattern
clearvars

Nchans = 48;
Fs = 12500;         % [Hz] Sampling frequency
d   = 0.1406525;    % [m] spacing
p   = [zeros(Nchans,1) zeros(Nchans,1) (0:(Nchans)-1)'*d].'; % position vector[m]
c   = 1500; % [m/s]
f0  = 2025; % [Hz] frequency of interest

% Steering angles
phi         = (-180:180).*pi/180;    % [rad]
l_phi       = length(phi);
theta       = (-180:180).*pi/180;      % [rad]
l_theta     = length(theta);

ux = repmat(cos(theta),[1 length(phi)]);
uy = sin(theta).'*sin(phi);
uz = sin(theta).'*cos(phi);

% ux = sin(phi).'*cos(theta);
% uy = cos(phi).'*cos(theta);
% uz = ones(1,l_phi)'*sin(theta);

ux = ux(:);
uy = uy(:);
uz = uz(:);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);

up = u.'*p;
kp = 2*pi/c*f0*up.';

amplitude = 1;

% Target info
phi_tgt = 30*pi/180;
theta_tgt = 90*pi/180;
[ h_data ] = MakeArrayData(p,Nchans,f0,c,Fs,theta_tgt,phi_tgt,amplitude);

t = (0:length(h_data)-1)/Fs;
Nfft = 2^nextpow2(length(h_data));
f_axis     = Fs*(0:Nfft/2-1)/Nfft;

H_freq = fft(h_data,Nfft,2); H_freq = H_freq(:,1:Nfft/2);

H_snaps = phase(H_freq(:,round(f0*Nfft/Fs)));
X_snaps = (H_snaps)*cos(theta_tgt);
Y_snaps = (H_snaps)*sin(theta_tgt)*sin(phi_tgt);
Z_snaps = (H_snaps)*sin(theta_tgt)*cos(phi_tgt);

alpha = .5;

% XYZsnaps = [X_snaps,Y_snaps,Z_snaps];
Xweighted = X_snaps*ux.';Yweighted = Y_snaps*uy.';Zweighted = Z_snaps*uz.';
XYZweighted = Xweighted + Yweighted + Zweighted;

% VS = alpha*H_snaps + (1-alpha)*XYZweighted;
% VSdelayed = bsxfun(@minus,VS,kp);
H_delayed = bsxfun(@minus,H_snaps,kp);
% VSresponse = alpha*H_snaps + (1-alpha).*XYZcombined;
figure(1)
clf
% plot(f_axis,10*log10(H_freq.*conj(H_freq)))
plot(t,h_data)
% Response = exp(-1i*Hdelayed);
Response = exp(-1i*(alpha*H_delayed+(1-alpha)*XYZweighted));
SumResponse = sum(Response,1);
BeamformedH = reshape(SumResponse,l_theta,l_phi);
figure(2)
clf
imagesc(theta*180/pi,phi*180/pi,10*log10(BeamformedH.*conj(BeamformedH)));
zlim([-40 0])
xlabel('$\phi$'); ylabel('$\theta$')

theta_ind = find(theta==theta_tgt);
phi_ind = find(phi==phi_tgt);
figure(3)
subplot(211)
plot(phi*180/pi,20*log10(abs(BeamformedH(theta_ind,:))))
xlabel(['$\phi$,$\theta=$' num2str(theta_tgt*180/pi)])
subplot(212)
plot(theta*180/pi,20*log10(abs(BeamformedH(:,phi_ind))))
xlabel(['$\theta$, $\phi=$' num2str(phi_tgt*180/pi)])