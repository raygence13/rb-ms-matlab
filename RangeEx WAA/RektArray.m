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
theta0 =  60*pi/180;
% D/E MRA [rad]
phi0 = 90*pi/180;
% Azimuth MRA [rad]

%%%%-----------------------------------------------------------------------
%%%%    Initializing array parameters
%%%%-----------------------------------------------------------------------
Lx = 20;
Ly = 20;
L = Lx*Ly;

%%%%-----------------------------------------------------------------------
%%%%    Defining element locations on first ring
%%%%-----------------------------------------------------------------------
d = lambda/2;  % spacing for x-direction of 1st planar array

px = (0:Lx-1)'*d;
p = repmat(px,[Ly 1]);

for ind = 1:Ly
    p((1:Lx) + (ind-1)*Lx,2) = (ind-1)*d;
end
p = p.';
figure(1)
hold on
plot(p(1,:),p(2,:),'bo')
axis equal
%% Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Initializing sweeping source direction variables
%%%%-----------------------------------------------------------------------
phi = (0:180).*pi/180;    % phi vector to plot BP against [rad]
l_phi = length(phi);
theta = (0:180).*pi/180;    % theta vector to plot BP against [rad]
l_theta = length(theta);

k0 = 2*pi/lambda.*...           % Wavenumber LOOK direction [rad/m] 
    [sin(theta0).*cos(phi0);    % x-component of look-direction
    sin(theta0).*sin(phi0)];     % y-component of look-direction
%     cos(theta0);];...           % z-component of look-direction
             % [3 x length(phi)]
k0p1 = k0'*p;

%%%%-----------------------------------------------------------------------
%%%%    VECTORIZED BEAMFORMER
%%%%-----------------------------------------------------------------------
% p1_norm = sqrt(sum(p1.^2,1));
% |p| for finding cos(ang) to implement baffling, [1 x N*M]


ux = sin(theta).'*cos(phi);
uy = sin(theta).'*sin(phi);
% uz = repmat(cos(theta),[l_phi 1]);
u(1,:) = ux(:);
u(2,:) = uy(:);
% u(3,:) = uz(:);
tic
kp1 = 2*pi/c*f.*(u'*p);
toc
%%
%     k_norm = sqrt(sum((2*pi/c*f.*u).^2,1));
    % |k| for finding cos(ang) to implement baffling, [1 x length(phi)]
    
%     kp1_norm = k_norm.'*p1_norm;  % |k||p|, [length(phi) x N*M]
%     cos_er = (kp1)./kp1_norm; % (k*p) / (|k||p|), [length(phi) x N*M]
    % calculates cos(ang) between element and source for element resp.
%     cos_er(acos(cos_er) > (pi/2)) = 0;
    % including baffling
delay1 = bsxfun(@minus,kp1,k0p1);
%  Response = cos_er.*exp(-1i.*delay);
 Response1 = exp(-1i.*delay1);
% 
Beampattern1 = mean(Response1,2);
Beampattern1 = Beampattern1./max(Beampattern1);
Beampattern1 = reshape(Beampattern1,l_theta,l_phi);    % sum portion of beamforming


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
figure(2)
% imagesc(phi*180/pi,theta*180/pi,20*log10(abs(Beampattern1)).')
imagesc(theta.*180/pi,phi.*180/pi,20*log10(abs(Beampattern1)).')

set(gca,'YDir','normal')
colorbar
lowerlim = -40; upperlim = 0;
caxis([lowerlim upperlim]);
xlabel('\theta'); ylabel('\phi'); zlabel('[dB]')
title('Beampattern Plot')

%%%%-----------------------------------------------------------------------
%%%%    1D-Slice BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(5)
subplot(211)
plot(theta.*180/pi,20*log10(abs(Beampattern1(phi_LDI,:))).')
xlabel(['\theta in [deg], for \phi_0 = ' num2str(phi(phi_LDI)*180/pi) '\circ'])
ylabel('[dB]')
title('D/E and AZ Slice Plots for Look Direction')
axis([theta(1)*180/pi theta(end)*180/pi lowerlim upperlim])

subplot(212)
plot(phi.*180/pi,20*log10(abs(Beampattern1(:,theta_LDI))).')
axis([phi(1)*180/pi phi(end)*180/pi lowerlim upperlim])
xlabel(['\phi in [deg], for \theta_0 = ' num2str(theta(theta_LDI)*180/pi) '\circ'])
ylabel('[dB]')

profile report
