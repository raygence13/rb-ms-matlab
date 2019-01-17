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
phi0 = 0*pi/180;
% Azimuth MRA [rad]

%%%%-----------------------------------------------------------------------
%%%%    Initializing array parameters
%%%%-----------------------------------------------------------------------
M1 = 9;
N1 = 7;
M2 = 4;
N2 = 5;

%%%%-----------------------------------------------------------------------
%%%%    Defining element locations on first ring
%%%%-----------------------------------------------------------------------
dx1 = N2*lambda/2;  % spacing for x-direction of 1st planar array
px1 = (0:N1-1)*dx1;
px1 = repmat(px1,[1 M1]);

dy1 = M2*lambda/2;  % spacing for y-direction of 1st planar array
py1 = (0:M1-1)*dy1;
py1 = repmat(py1,[1 N1]);

p1 = [px1; py1; zeros(1,M1*N1)];
% p1 = [px1;py1];

dx2 = N1*lambda/2;  % spacing for x-direction of 2nd planar array
px2 = (0:N2-1)*dx2;
px2 = repmat(px2,[1 M2]);

dy2 = M1*lambda/2;  % spacing for y-direction of 2nd planar array
py2 = (0:M2-1)*dy2;
py2 = repmat(py2,[1 N2]);
p2 = [px2;py2;zeros(1,M2*N2)];
% p2 = [px2;py2];

figure(1)
hold on
plot(p1(1,:),p1(2,:),'bo')
plot(p2(1,:),p2(2,:),'go')
%% Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Initializing sweeping source direction variables
%%%%-----------------------------------------------------------------------
phi = (-180:180).*pi/180;    % phi vector to plot BP against [rad]
l_phi = length(phi);
theta = (-180:180).*pi/180;    % theta vector to plot BP against [rad]
l_theta = length(theta);

k0 = 2*pi/lambda.*...           % Wavenumber LOOK direction [rad/m] 
    [sin(theta0).*cos(phi0);    % x-component of look-direction
    sin(theta0).*sin(phi0);     % y-component of look-direction
    cos(theta0);];...           % z-component of look-direction
             % [3 x length(phi)]
k0p1 = k0.'*p1;
k0p2 = k0.'*p2;

%%%%-----------------------------------------------------------------------
%%%%    VECTORIZED BEAMFORMER
%%%%-----------------------------------------------------------------------
% p1_norm = sqrt(sum(p1.^2,1));
% |p| for finding cos(ang) to implement baffling, [1 x N*M]


ux = sin(theta).'*cos(phi);
uy = sin(theta).'*sin(phi);
uz = repmat(cos(theta),[l_phi 1]);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);
tic
kp1 = 2*pi/c*f.*(u.'*p1);
kp2 = 2*pi/c*f.*(u.'*p2);
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
delay2 = bsxfun(@minus,kp2,k0p2);
%  Response = cos_er.*exp(-1i.*delay);
 Response1 = exp(-1i.*delay1);
 Response2 = exp(-1i.*delay2);
%  Response = bsxfun(@times,W.',Response);
%     Response = bsxfun(@times,Response,W); % [length(phi) x N*M]
    % weighting/baffling/delay portion of beamforming
%     BP = sum(Response,1)/(N*M);
% 
Beampattern1 = mean(Response1,2);
Beampattern1 = Beampattern1./max(Beampattern1);
Beampattern1 = reshape(Beampattern1,l_phi,l_theta);    % sum portion of beamforming

Beampattern2 = mean(Response2,2);
Beampattern2 = Beampattern2./max(Beampattern2);
Beampattern2 = reshape(Beampattern2,l_phi,l_theta);    % sum portion of beamforming

Beampattern = Beampattern1.*conj(Beampattern2);
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
figure(2)
imagesc(phi*180/pi,theta*180/pi,20*log10(abs(Beampattern1)).')

set(gca,'YDir','normal')
colorbar
lowerlim = -40; upperlim = 0;
caxis([lowerlim upperlim]);
xlabel('\phi'); ylabel('\theta'); zlabel('[dB]')
title('Beampattern Plot')

%%%%-----------------------------------------------------------------------
%%%%    BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(3)
imagesc(phi*180/pi,theta*180/pi,20*log10(abs(Beampattern2)).')

set(gca,'YDir','normal')
colorbar
lowerlim = -40; upperlim = 0;
caxis([lowerlim upperlim]);
xlabel('\phi'); ylabel('\theta'); zlabel('[dB]')
title('Beampattern Plot')
%%%%-----------------------------------------------------------------------
%%%%    BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(4)
imagesc(phi*180/pi,theta*180/pi,20*log10(abs(Beampattern)).')

set(gca,'YDir','normal')
colorbar
lowerlim = -40; upperlim = 0;
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
