%% WAA Beampattern Script
% Ray Bautista
% 27July2016
%%%%-----------------------------------------------------------------------
%%%%    Standard MATLAB Coordinate System
%%%%    U = [sin(AZ)cos(DE);
%%%%        cos(AZ)cos(DE);
%%%%        sin(DE)];
%%%%-----------------------------------------------------------------------

clearvars
close all
% clc
present(0)

profile on
%% Plane Wave Parameters
f = 2000;       % Frequency of plane wave [Hz]
c = 1500;       % Speed of sound [m/s]
lambda = c/f;   % Wave length [m]

%% Array Parameters for Cylinder Sector of M Rings with N elements per Ring

%%%%-----------------------------------------------------------------------
%%%%    Look Direction
%%%%-----------------------------------------------------------------------
theta0 =  0*pi/180;
% D/E MRA [rad]
phi0 = 90*pi/180;
% Azimuth MRA [rad]

%%%%-----------------------------------------------------------------------
%%%%    Initializing array parameters
%%%%-----------------------------------------------------------------------
M = 40;                 % Number of Rings in Cylinder
p1 = [zeros(1,M);(-(M-1)/2:(M-1)/2)*lambda/2;zeros(1,M)];
p2 = [zeros(1,M);zeros(1,M);(-(M-1)/2:(M-1)/2)*lambda/2];
% p3 = [(-(M-1)/2:(M-1)/2)*lambda/2;zeros(1,M);zeros(1,M)];
p = [p1,p2];
figure
plot3(p(1,:),p(2,:),p(3,:),'.');
hold on
plot3(p2(1,:),p2(2,:),p2(3,:),'.');
%% Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Initializing sweeping source direction variables
%%%%-----------------------------------------------------------------------
phi = (-180:.25:180).*pi/180;    % phi vector to plot BP against [rad]
l_phi = length(phi);
theta = (-90:.25:90).*pi/180;    % theta vector to plot BP against [rad]
l_theta = length(theta);

k0 = 2*pi/lambda.*...    % Wavenumber LOOK direction [rad/m] 
    [sin(phi0).*cos(theta0);    % x-component of look-direction
    cos(phi0).*cos(theta0);     % y-component of look-direction
    sin(theta0);];...            % z-component of look-direction
             % [3 x length(phi)]
k0p = k0.'*p;


ux = sin(phi).'*cos(theta);
uy = cos(phi).'*cos(theta);
uz = repmat(sin(theta),[l_phi 1]);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);
tic
kp = 2*pi/c*f.*(u.'*p);
k_norm = sqrt(sum((2*pi/c*f.*u).^2,1));
p_norm = sqrt(sum(p.^2,1));
kp_norm = k_norm.'*p_norm;

cos_er = kp./kp_norm;
cos_er(acos(cos_er) > pi/2) = 0;
cos_er(cos_er~=0) = 2;
toc
%%

delay = bsxfun(@minus,kp,k0p);
Response = cos_er.*exp(-1i.*delay);
%     Response = bsxfun(@times,Response,W); % [length(phi) x N*M]
    % weighting/baffling/delay portion of beamforming
%     BP = sum(Response,1)/(N*M);
% 
Beampattern = mean(Response,2);
Beampattern = Beampattern./max(Beampattern);
Beampattern = reshape(Beampattern,l_phi,l_theta);    % sum portion of beamforming
% end
toc

%% Plotting Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Finding look-direction index to plot DE/AZ slices
%%%%-----------------------------------------------------------------------
if phi0<0
    phi_LDI = find(abs(phi)<=abs(phi0),1);
else
    phi_LDI = find(phi>=phi0,1);
end

if theta0<0
    theta_LDI = find(abs(theta)<=abs(theta0),1); 
else
    theta_LDI = find(theta>=theta0,1);
end

%%%%-----------------------------------------------------------------------
%%%%    BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(4)
imagesc(phi*180/pi,theta*180/pi,20*log10(abs(Beampattern)).')

set(gca,'YDir','normal')
colorbar
lowerlim = -12; upperlim = 20*log10(max(max(abs(Beampattern))));
caxis([lowerlim upperlim]);
xlabel('\phi'); ylabel('\theta'); zlabel('[dB]')
title('Beampattern Plot')

%%%%-----------------------------------------------------------------------
%%%%    1D-Slice BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(5)
subplot(211)
plot(theta.*180/pi,20*log10(abs(Beampattern(phi_LDI,:))).')
xlabel(['\theta in [deg], for \phi_0 = ' num2str(phi(phi_LDI)*180/pi) '\circ'])
ylabel('[dB]')
title('D/E and AZ Slice Plots for Look Direction')
axis([theta(1)*180/pi theta(end)*180/pi lowerlim upperlim])

subplot(212)
plot(phi.*180/pi,20*log10(abs(Beampattern(:,theta_LDI))).')
axis([phi(1)*180/pi phi(end)*180/pi lowerlim upperlim])
xlabel(['\phi in [deg], for \theta_0 = ' num2str(theta(theta_LDI)*180/pi) '\circ'])
ylabel('[dB]')

profile report
