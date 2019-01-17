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
theta0 =  -60*pi/180;
% D/E MRA [rad]
phi0 = 90*pi/180;
% Azimuth MRA [rad]

%%%%-----------------------------------------------------------------------
%%%%    Initializing array parameters
%%%%-----------------------------------------------------------------------
M = 20;                 % Number of Rings in Cylinder
N = 30;                 % Number of Elements per Ring
sector = 60*pi/180;     % Sector of circle array covers [rad]
dphi = sector/N;        % Element angular spacing [rad]
r = lambda/2/dphi;      % Radius to gaurantee lambda/2 arc spacing [m]

DEshift = 70*pi/180;
% D/E angular shift [rad]
shift = (sector+DEshift-dphi)/2;
% Gaurantees array is centered around DEshift/2 [rad]

%%%%-----------------------------------------------------------------------
%%%%    Defining element locations on first ring
%%%%-----------------------------------------------------------------------
p = [r.*cos((0:N-1).*dphi-shift);...    % x-component of element position
    -ones(1,N)*(M-1)*lambda/4;...       % y-component of element position
    r.*sin((0:N-1).*dphi-shift)];       % z-component of element position
% First ring element position vector [m],[3 x N]

FC = [r.*cos(0:.05:2*pi);...
    -ones(1,126)*(M-1)*lambda/4;...
    r.*(sin(0:.05:2*pi))];
% Equivalent radius Full Cylinder

% 1D-Plot to ensure elements sit along circle
figure(1)
plot(FC(1,:),FC(3,:),'--')
hold on
plot(p(1,:,1),p(3,:,1),'k.')
axis equal
xlabel('x-direction'); ylabel('z-direction');
title('1D-Plot of Element Positions')

%%%%-----------------------------------------------------------------------
%%%%    Creating cylinder of M rings
%%%%    Vectorized element positions
%%%%    Every p(2,ring:ring+N-1) contains incrementing ring space in y-dir
%%%%-----------------------------------------------------------------------
tic
if M > 1
    p = repmat(p,[1 M]);
    FC = repmat(FC,[1 1 M]);
    figure(2)
    
    for ring = 1:M-1
        FC(2,:,ring+1) = FC(2,:,ring) + lambda/2; % adding ring space
        plot3(FC(1,:,ring),FC(2,:,ring),FC(3,:,ring),':')
        axis equal
        hold on
    end
    
    for ring = N:N:N*M-1
        p(2,ring+1:ring+N) = p(2,ring+1-N) + lambda/2; % adding ring space
        % Full position vector [m], [3 x N*M]
        plot3(p(1,ring+1-N:ring),p(2,ring+1-N:ring),p(3,ring+1-N:ring),'k.')
        hold on
        
    end
    xlabel('x in [m]'); ylabel('y in [m]'); zlabel('z in [m]')
    title('Element Positions')
end
toc

%% Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Initializing sweeping source direction variables
%%%%-----------------------------------------------------------------------
phi = (0:.5:180).*pi/180;    % phi vector to plot BP against [rad]
l_phi = length(phi);
theta = (-90:.5:50).*pi/180;    % theta vector to plot BP against [rad]
l_theta = length(theta);

k0 = 2*pi/lambda.*...    % Wavenumber LOOK direction [rad/m]
    [sin(phi0).*cos(theta0);    % x-component of look-direction
    cos(phi0).*cos(theta0);     % y-component of look-direction
    sin(theta0);];...            % z-component of look-direction
    % [3 x length(phi)]
k0p = k0.'*p;
%%%%-----------------------------------------------------------------------
%%%%    Shading
%%%%-----------------------------------------------------------------------
wr = taylorwin(N);  % shading in row direction, [1 x N]
wc = taylorwin(M);  % shading in column direction, [1 x M]
w = wr*wc.';        % calculates weights for entire array [N x M]

% W = repmat(w(:).',l_phi*l_theta,1);
W = w(:);
% vectorized window for efficiency [length(phi) x N*M]

%%%%-----------------------------------------------------------------------
%%%%    VECTORIZED BEAMFORMER
%%%%-----------------------------------------------------------------------
p_norm = sqrt(sum(p.^2,1));
% |p| for finding cos(ang) to implement baffling, [1 x N*M]


ux = sin(phi).'*cos(theta);
uy = cos(phi).'*cos(theta);
uz = repmat(sin(theta),[l_phi 1]);
u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);
tic
kp = 2*pi/c*f.*(u.'*p);
toc
%%
k_norm = sqrt(sum((2*pi/c*f.*u).^2,1));
% |k| for finding cos(ang) to implement baffling, [1 x length(phi)]

kp_norm = k_norm.'*p_norm;  % |k||p|, [length(phi) x N*M]
cos_er = (kp)./kp_norm; % (k*p) / (|k||p|), [length(phi) x N*M]
% calculates cos(ang) between element and source for element resp.
cos_er(acos(cos_er) > (pi/2)) = 0;
% including baffling
delay = bsxfun(@minus,kp,k0p);
Response = cos_er.*exp(-1i.*delay);
Response = bsxfun(@times,W.',Response);
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
lowerlim = -70; upperlim = 20*log10(max(max(abs(Beampattern))));
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
