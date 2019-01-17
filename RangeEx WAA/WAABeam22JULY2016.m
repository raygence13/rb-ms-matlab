%% WAA Beampattern Script
% Ray Bautista
% 22July2016
%%%%-----------------------------------------------------------------------
%%%%    Standard MATLAB Coordinate System
%%%%    U = [sin(AZ)cos(DE);
%%%%        cos(AZ)cos(DE);
%%%%        sin(DE)];
%%%%-----------------------------------------------------------------------

clearvars
close all
clc
present(0)

profile on
%% Plane Wave Parameters
f0 = 2000;      % Frequency of plane wave [Hz]
c = 1500;       % Speed of sound [m/s]
lambda = c/f0;  % Wave length [m]

%% Array Parameters for Cylinder Sector of M Rings with N elements per Ring

%%%%-----------------------------------------------------------------------
%%%%    Look Direction
%%%%-----------------------------------------------------------------------
disp('Array configured such that phi=90[deg] is normal (broadside) to array')

phi0 = input('Azimuth MRA between +- 180[deg]= ')*pi/180;
while abs(phi0)>pi      % While statement to ensure user input is satisfactory
    phi0 =  input('Azimuth MRA between +- 180[deg] = ')*pi/180;
end
% Azimuth MRA [rad]

theta0 =  input('Depth/Elevation MRA between +- 90[deg] = ')*pi/180;
while abs(theta0)>pi/2  % While statement to ensure user input is satisfactory
    theta0 =  input('Depth/Elevation MRA between +- 90[deg] = ')*pi/180;
end
% D/E MRA [rad]


%%%%-----------------------------------------------------------------------
%%%%    Initializing array parameters
%%%%-----------------------------------------------------------------------
M = 32;                 % Number of Rings in Cylinder
N = 14;                 % Number of Elements per Ring
sector = 60*pi/180;     % Sector of circle array covers [rad]
ring_space = lambda/2;
dphi = sector/N;        % Element angular spacing [rad]
r = lambda/2/dphi;      % Radius to gaurantee lambda/2 arc spacing [m]


DEshift = -input('DE Angular shift of array (for best response, input DE of MRA): ')*pi/90;
% D/E angular shift [rad]
shift = (sector+DEshift-dphi)/2;
% Gaurantees array is centered around DEshift/2 [rad]
%% Creating Element Positions
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

%%%%-----------------------------------------------------------------------
%%%%    Creating cylinder of M rings with FAILED elements
%%%%    Vectorized element positions
%%%%    Every p(2,ring:ring+N-1) contains incrementing ring space in y-dir
%%%%-----------------------------------------------------------------------
failure_type = input('Select failure type: none(1), random(2), or section(3): ');
switch failure_type
    case 1  % NO Failures
        failed = ones(1,N*M);
    
    case 2  % Random failed elements
        q = input('Probability of a failed sensor: q = '); % Probability of a failed sensor
        while q > 1
            q = input('q must be less than 1: q = ');
        end
        failed = rand(1,N*M) > q;
        disp(['Total number of failed elements = ' num2str(N*M - sum(failed))])
    
    case 3  % Section of failed elements
        disp('Total number of failed elements = num_c*num_r')
        num_c = input('Number of columns in failed section: num_c = ');
        while num_c > M
            num_c = input(['Must be less than ' num2str(M) ]);
        end
        
        num_r = input('Number of rows in failed section: num_r = ');
        while num_r > N
            num_r = input(['Must be less than ' num2str(N) ]);
        end
        col_failed = randi([0,M-num_c])*N;
        % column index for beginning of failed column elements
        row_failed = randi([0,N-num_r]);
        % row index for beginning of failed row elements
        start_failed = col_failed+row_failed;
        
        failed = ones(1,N*M);
        for f = 1:num_c
            failed( (start_failed + 1:start_failed+num_r) + (f-1)*N ) = 0;
            % Inserting num_r 0's every Nth index, looped num_c times
        end
        
    otherwise  % NO Failures
        failed = ones(1,N*M);
end

% 1D-Plot to ensure elements sit along circle
figure(1)
plot(FC(1,:),FC(3,:),'--')
hold on
plot(p(1,:,1),p(3,:,1),'k.')
axis equal
xlabel('x-direction'); ylabel('z-direction');
title('1D-Plot of Element Positions')

if M > 1
    p = repmat(p,[1 M]);
    
    FC = repmat(FC,[1 1 M]);
    figure(2)
    
    for ring = 1:M
        FC(2,:,ring+1) = FC(2,:,ring) + lambda/2; % adding ring space
        plot3(FC(1,:,ring),FC(2,:,ring),FC(3,:,ring),':')
        axis equal
        hold on
    end
    
    for ring = N:N:N*M-1
        p(2,ring+1:ring+N) = p(2,ring+1-N) + lambda/2; % adding ring space
        % Full position vector [m], [3 x N*M]
    end
    
    p = bsxfun(@times,p,failed);    % Removing elements
    for ring = N:N:N*M
        plot3(p(1,ring+1-N:ring),p(2,ring+1-N:ring),p(3,ring+1-N:ring),'k.')
        hold on
    end
    xlabel('x in [m]'); ylabel('y in [m]'); zlabel('z in [m]')
    title('Element Positions')
end

%% Beampattern

%%%%-----------------------------------------------------------------------
%%%%    Initializing sweeping source direction variables
%%%%-----------------------------------------------------------------------
phi = (-180:.5:180).*pi/180;    % phi vector to plot BP against [rad]
theta = (-90:.5:90).*pi/180;    % theta vector to plot BP against [rad]

%%%%-----------------------------------------------------------------------
%%%%    Shading
%%%%-----------------------------------------------------------------------
window = input('Array Shading: Rectangular(1), Taylor(2), Chebyshev(3): ');
switch window
    case 1      % Rectangular window (no shading)
        wr = ones(N,1);     % shading in row direction, [N x 1]
        wc = ones(M,1);     % shading in column direction, [M x 1]
    case 2      % Taylor window
        wr = taylorwin(N);  % shading in row direction, [N x 1]
        wc = taylorwin(M);  % shading in column direction, [M x 1]
    case 3      % Chebyshev window
        R = abs(input('Side lobe level [dB]= '));
        wr = chebwin(N,R);  % shading in row direction, [N x 1]
        wc = chebwin(M,R);  % shading in column direction, [M x 1]
    otherwise   % Rectangular window (no shading)
        wr = ones(N,1);
        wc = ones(M,1);
end

w = wr*wc.';
W = repmat(w(:).',length(phi),1);
% vectorized window for efficiency [length(phi) x N*M]

%%%%-----------------------------------------------------------------------
%%%%    VECTORIZED BEAMFORMER
%%%%-----------------------------------------------------------------------
f_rel = .5:.5:2.5; % relative frequency
BP = zeros(length(theta),length(phi),length(f_rel)); 
% Initialize BeamPattern matrix
p_norm = sqrt(sum(p.^2,1));
% |p| for finding cos(ang) to implement baffling, [1 x N*M]

tic
for f = 1:length(f_rel)         % Sweeping through frequency
    for DE = 1:length(theta)    % Sweeping through DE, looking at all AZ
        k0 = repmat(2*pi/(lambda/f_rel(f)).*...% Wavenumber LOOK direction [rad/m] 
        [sin(phi0).*cos(theta0);    % x-component of look-direction
        cos(phi0).*cos(theta0);     % y-component of look-direction
        sin(theta0);]...            % z-component of look-direction
        ,[1 length(phi)]);          % [3 x length(phi)]
    
        k = 2*pi/(lambda/f_rel(f)).*...     % Wavenumber SOURCE direction
            [sin(phi)*cos(theta(DE));       % x-component of source-direction
             cos(phi)*cos(theta(DE));       % y-component of source-direction
             repmat(sin(theta(DE))...       % z-component of source-direction
            ,[1 length(phi)])];             % [3 x length(phi)]
        k_norm = sqrt(sum(k.^2,1));
        % |k| for finding cos(ang) to implement baffling, [1 x length(phi)]
        
        kp_norm = k_norm.'*p_norm;  % |k||p|, [length(phi) x N*M]
        cos_er = (k.'*p)./kp_norm;  % (k*p) / (|k||p|), [length(phi) x N*M]
        % calculates cos(ang) between element and source for element resp.
        cos_er(acos(cos_er) > (pi/2)) = 0;
        % including baffling
        cos_er(isnan(cos_er)) = 0;
        
        Response = W.*cos_er.*exp(-1i*(k-k0).'*p); % [length(phi) x N*M]
        % weighting/baffling/delay portion of beamforming
        BP(DE,:,f) = mean(Response,2);
        % sum portion of beamforming
    end
    BP(:,:,f) = BP(:,:,f)./max(max(BP(:,:,f)));
    % Normalizing BP
end
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


f_ind = input(['Plot frequency index between 1 and ' num2str(length(f_rel)) ': ']);
% User defined for plotting certain frequency

%%%%-----------------------------------------------------------------------
%%%%    Calculating Pointing Error
%%%%-----------------------------------------------------------------------
BPtemp = BP(:,:,f_ind);
% Used to calculate pointing error, sometimes not in LOOK direction
[~,BPidx] = max(BPtemp(:));
% Finding peak response in BP
[theta_max_idx,phi_max_idx] = ind2sub(size(BPtemp),BPidx);
% Finding theta and phi index for peak response to calculate pointing error

theta_pointing_error = abs(theta(theta_LDI)-theta(theta_max_idx))*180/pi;
phi_pointing_error = abs(phi(phi_LDI)-phi(phi_max_idx))*180/pi;

%%%%-----------------------------------------------------------------------
%%%%    Calculating Directivity Index
%%%%-----------------------------------------------------------------------
DI = computeDI(abs(BPtemp).^2,theta,phi);
disp(['DI(theta0 = ' num2str(theta0*180/pi)...
    ',phi0 = ' num2str(phi0*180/pi)...
    ',f = ' num2str(f0*f_rel(f_ind))...
    ') = ' num2str(DI) ' dB'])

%%%%-----------------------------------------------------------------------
%%%%    Finding Max Sidelobe height level
%%%%-----------------------------------------------------------------------
lowerlim = -70;
threshold = input(['Peak Side Lobe Height Threshold between 0 and ' num2str(lowerlim) '[dB]: ']);
while lowerlim > threshold || threshold > 0
    threshold = input(['Peak Side Lobe Height Threshold between 0 and ' num2str(lowerlim) '[dB]: ']);
end
NPeaks = input('Number of Side Lobes to find: ')+1;

[theta_pks,theta_ind] = findpeaks(20*log10(abs(BP(:,phi_LDI,f_ind))),...
    'MinPeakHeight',threshold,'SortStr','descend','NPeaks',NPeaks);
disp(['Theta Peak Side Lobe Levels [dB]= ' num2str(theta_pks(2:end).')])
disp([' at theta [deg]= ' num2str(theta(theta_ind(2:end).')*180/pi)]);

[phi_pks,phi_ind] = findpeaks(20*log10(abs(BP(theta_LDI,:,f_ind))),...
    'MinPeakHeight',threshold,'SortStr','descend','NPeaks',NPeaks);
disp(['Phi Peak Side Lobe Levels [dB]= ' num2str(phi_pks(2:end))])
disp([' at phi [deg]= ' num2str(phi(phi_ind(2:end))*180/pi)]);


%%%%-----------------------------------------------------------------------
%%%%    BP plot vs DE/AZ
%%%%-----------------------------------------------------------------------
figure(4)
imagesc(phi*180/pi,theta*180/pi,20*log10(abs(BPtemp)))
set(gca,'YDir','normal')
colorbar
upperlim = 20*log10(max(max(abs(BPtemp))));
caxis([lowerlim upperlim]);
xlabel('\phi'); ylabel('\theta'); zlabel('[dB]')
title(['Beampattern Plot for f = ' num2str(f0*f_rel(f_ind)) ' Hz'])

%%%%-----------------------------------------------------------------------
%%%%    1D-Slice BP(theta0,phi0,f_ind) plot vs DE and vs AZ
%%%%-----------------------------------------------------------------------
figure(5)
subplot(211)
plot(theta.*180/pi,20*log10(abs(BPtemp(:,phi_LDI))));
hold on
p5b = plot(theta(theta_ind(2:end))*180/pi,theta_pks(2:end),'kv','MarkerSize',5);
xlabel(['\theta in [deg], for \phi_0 = ' num2str(phi(phi_LDI)*180/pi) '\circ'])
ylabel('[dB]')
title(['D/E and AZ Slice Plots for Look Direction for f = ' num2str(f0*f_rel(f_ind)) ' Hz'])
axis([-90 90 lowerlim upperlim])
legend(p5b,'Peak Side Lobes')

subplot(212)
plot(phi.*180/pi,20*log10(abs(BPtemp(theta_LDI,:))));
hold on
p5d = plot(phi(phi_ind(2:end))*180/pi,phi_pks(2:end),'kv','MarkerSize',5);
axis([-180 180 lowerlim upperlim])
xlabel(['\phi in [deg], for \theta_0 = ' num2str(theta(theta_LDI)*180/pi) '\circ'])
ylabel('[dB]')
profile report

%%%%-----------------------------------------------------------------------
%%%%    1D-Slice BP(theta_max,phi_max,f_ind) plot vs DE and vs AZ
%%%%-----------------------------------------------------------------------
if theta_pointing_error > 0 || phi_pointing_error > 0
figure(6)
subplot(211)
plot(theta.*180/pi,20*log10(abs(BP(:,phi_max_idx))));
xlabel(['\theta in [deg], for \phi_{MAX} = ' num2str(phi(phi_max_idx)*180/pi) '\circ'])
ylabel('[dB]')
title(['D/E and AZ Slice Plots for Max Response for f = ' num2str(f0*f_rel(f_ind)) ' Hz'])
axis([-90 90 lowerlim upperlim])

subplot(212)
p5c = plot(phi.*180/pi,20*log10(abs(BPtemp(theta_max_idx,:))));
axis([-180 180 lowerlim upperlim])
xlabel(['\phi in [deg], for \theta_{MAX} = ' num2str(theta(theta_max_idx)*180/pi) '\circ'])
ylabel('[dB]')
end

profile report