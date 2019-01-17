% VSLA beamformer
clear;

%%%%% Flags

show_scan_plots=1;

%%%%% Environment and array paramters

% Speed of sound (m/s)
c=1500;
% Design frequency (Hz) and wavelength (m)
f_design=400; lambda_design=c/f_design;
% Number of elements
N_e=20;
% Element positions (uniformaly spaced at one half design wavelength
% along z-axis; assume hydrophones and geophones are colocated)
p_e=[zeros(1,N_e); zeros(1,N_e); 0.5*lambda_design*(0:N_e-1)];
% Center array about the origin (for aesthetics only)
p_e(3,:)=p_e(3,:)-(p_e(3,N_e)-p_e(3,1))/2;

%%%%% Targer parameters

% Target frequency (Hz) and amplitude
f_tgt=150; A_tgt=8;
% Target course (deg -> rad) and speed (m/s)
cse= -45*pi/180; spd=5;
% Initial range (m) and bearing (deg -> rad)
rng(1)=3000; brg(1)=95*pi/180;
% Sceneario length (min -> sec)
T=10*60;
% Scan interval (sec) and times (sec)
dt=4; t=0:dt:T; N_t=length(t); T=t(N_t);
% Target positions (assumes target in y-z plane)
p_tgt=repmat([0;rng(1)*[sin(brg(1));cos(brg(1))]],[1,N_t])+...
    [0;spd*[sin(cse);cos(cse)]]*(t-t(1));
% Range and bearing to target
rng=sqrt(p_tgt(2,:).^2+p_tgt(3,:).^2);
brg=atan2(p_tgt(2,:),p_tgt(3,:)); % TODO: map to [-pi, pi)

figure(1); clf;
hdl=plot(p_e(2,:),p_e(3,:),'b.',p_tgt(2,:),p_tgt(3,:),'r.');
axis equal; grid on;
xlabel('y (m)'); ylabel('z (m)');
title([int2str(N_e),' Element Stationary Line Array and CV Target']);
legend(hdl,{'Array','Target'});

figure(2); clf;
plot(brg*180/pi,t,'g-'); xlim([-1,1]*180); grid on;
xlabel('AZ (deg)'); ylabel('t (s)');
title('Azimuthal Bearing to Target');

%%%%% CBF parameters

% Sampling times for each scan
f_samp=2*f_design;
t_samp=linspace(0,dt-dt/(f_samp*dt-1),f_samp*dt);
N_samp=length(t_samp);

% FFT size and frequencies (Hz)
N_fft=2^nextpow2(length(t_samp));
f_fft=f_samp*(0:N_fft-1)'/N_fft;

% Number of beams
N_beams=51;
% Beam angles (assume equispaced in angle for now)
b_beams=linspace(-pi+2*pi/N_beams/2,pi-2*pi/N_beams/2,N_beams);
% Beam unit vectors
u_beams=[zeros(1,N_beams); sin(b_beams); cos(b_beams)];

% Window function for elment time series
w_fft=window(@hamming,N_samp);

% Taylor weights
w_e=taylorwin(N_e);

% Relative hydrophone/geophone weight (alpha = 1 for hydrophone only)
alpha=0.5; % 0.5;

%%%%% Run the CBF
tau=zeros(1,N_e);
h_e=zeros(N_samp,N_e); H_e=zeros(N_fft,N_e);
x_e=zeros(N_samp,N_e); X_e=zeros(N_fft,N_e);
y_e=zeros(N_samp,N_e); Y_e=zeros(N_fft,N_e);
z_e=zeros(N_samp,N_e); Z_e=zeros(N_fft,N_e);
bb=zeros(N_t,N_beams);
for n_t=1:N_t
    fprintf(1,'Scan %3i of %3i\n',n_t,N_t);
    
    %%%%% Generate element data (time series)
    
    % Unit vector to target
    u_tgt=p_tgt(:,n_t)/norm(p_tgt(:,n_t));
    % Generate element time series and transform to frequency domain
    for n_e=1:N_e
        % Time-of-arrival delay to element
        tau(n_e)=p_e(:,n_e)'*u_tgt/c;
        % Element time series (note: -tau gives mirror image)
        pressure_e=A_tgt*sin(2*pi*f_tgt*(t(n_t)+t_samp+tau(n_e)));
        h_e(:,n_e)=pressure_e; % "+ noise"
        x_e(:,n_e)=pressure_e*u_tgt(1); % "+ noise"
        y_e(:,n_e)=pressure_e*u_tgt(2); % "+ noise"
        z_e(:,n_e)=pressure_e*u_tgt(3); % "+ noise"
        % Transform to frequency domain
        H_e(:,n_e)=fft(w_fft.*h_e(:,n_e),N_fft);
        X_e(:,n_e)=fft(w_fft.*x_e(:,n_e),N_fft);
        Y_e(:,n_e)=fft(w_fft.*y_e(:,n_e),N_fft);
        Z_e(:,n_e)=fft(w_fft.*z_e(:,n_e),N_fft);
    end
    
    %%%%% Beamform element data
    
    fraz=zeros([N_fft,N_beams]);
    for n_b=1:N_beams
        for n_e=1:N_e
            tau_b=p_e(:,n_e)'*u_beams(:,n_b)/c;
            fraz(:,n_b)=fraz(:,n_b)+exp(-1i*2*pi*f_fft*tau_b).*(...
                alpha*w_e(n_e)*H_e(:,n_e)+(1-alpha)*(...
                    w_e(n_e)*X_e(:,n_e)*0+...
                    w_e(n_e)*Y_e(:,n_e)*sin(b_beams(n_b))+...
                    w_e(n_e)*Z_e(:,n_e)*cos(b_beams(n_b))));
        end
    end
    
    idx=1:N_fft/2;
    pwr=fraz(idx,:).*conj(fraz(idx,:));
    bb(n_t,:)=sum(pwr);
    
    if show_scan_plots
        % Plot 5 cycles of time series for 2 elements
        k=find(t_samp <= 5/f_tgt); e=[1,5];
        figure(91); clf;
        hdl=plot(t_samp(k),h_e(k,e(1)),'k-',t_samp(k),h_e(k,e(2)),'k--');
        xlabel('t (s)'); xlim([t_samp(1),t_samp(k(end))]);
        ylabel(['e_',int2str(e(1)),'(t), e_',int2str(e(2)),'(t)']);
        % Plot power spectrum on first element
        figure(92); clf;
        plot(f_fft(idx),10*log10(H_e(idx,1).*conj(H_e(idx,1))));
        xlim([f_fft(1),f_fft(N_fft/2)]);
        xlabel('f (Hz)'); ylabel('dB');
        % Plot FRAZ
        figure(93); clf;
        imagesc(b_beams*180/pi,f_fft(idx),10*log10(pwr));
        xlabel('AZ (deg)'); ylabel('f (Hz)');
        title(sprintf('Scan %3i of %3i\n',n_t,N_t));
        % Plot broadband output
        figure(94); clf;
        bar(b_beams*180/pi,10*log10(bb(n_t,:)),1); xlim([-180,180]);
        hold on; plot(brg(n_t)*[1,1]*180/pi,ylim,'r-o'); hold off;
        xlabel('AZ (deg)'); ylabel('dB');
        title(sprintf('Scan %3i of %3i\n',n_t,N_t));
        drawnow
    end
    
end

figure(3); clf;
imagesc(b_beams*180/pi,t,bb); set(gca,'YDir','normal');
hold on;
plot(brg*180/pi,t,'g-');
plot((-brg)*180/pi,t,'g--');
hold off;
colormap(1-gray);
xlabel('AZ (deg)'); ylabel('t (s)');
if alpha == 1
    title('Hydrophone Only (\alpha = 1)');
else
    title(['\alpha = ',num2str(alpha)]);
end

