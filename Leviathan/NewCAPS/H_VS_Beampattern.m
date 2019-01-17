function [BP] = H_VS_Beampattern(p,phi_vec,theta_vec,f0,phi0,theta0)
%[BP] = H_VS_Beampattern(p,phi_vec,theta_vec,phi0,theta0).
%Function to generate beampatterns using hydrophone only and vector
% sensors based on array geometry defined by p.AEL
% RB 9MAR2018. RB 15MAR2018.
% Inputs:
% p:                structure containing general parameters
% p.AEL:            Array Element Locations [3 x Nsensors] [m]
% p.Nsensors:       number of sensors
% p.MaxDistance:    max distance between sensor pairs in geometry [m]
% p.depth:          depth of array [m]
% p.c:              sound speed of water [m/s]
% p.rho:            density [kg/m^3]
% f0:               center frequency
%
% phi_vec:      user defined phi vector (deg)
% theta_vec:    user defined theta vector (deg)
% phi0:         Azimuth look direction (deg)
% theta0:       D/E look direction [sin(deg)]
%
% Outputs:
% BP:       structure containing beampatterns
% BP.H:     2D hydrophone beampattern
% BP.Haz:   1D hydrophone beampattern along theta0 slice
% BP.Hde:   1D hydrophone beampattern along phi0 slice
% BP.VS:    2D vector sensor beampattern
% BP.VSaz:  1D vector sensor beampattern along theta0 slice
% BP.VSde:  1D vector sensor beampattern along phi0 slice

if max(abs(theta_vec)) <= 1
    theta_vec = asind(theta_vec);
end
l_phi = length(phi_vec);
l_theta = length(theta_vec);

% steering vectors
ux = cosd(phi_vec)'*cosd(theta_vec);
ux = ux(:);     % x-component of unit vector
uy = sind(phi_vec)'*cosd(theta_vec);
uy = uy(:);     % z-component of unit vector
uz = ones(l_phi,1)*sind(theta_vec);
uz = uz(:);     % y-component of unit vector

% 3-dimensional unit vector
u = zeros(3,l_phi*l_theta);
u(1,:) = ux;
u(2,:) = uy;
u(3,:) = uz;

% finding indices for phi0 and theta0
if phi0<0
    phi_LDI = find(abs(phi_vec)<=abs(phi0),1);
else
    phi_LDI = find(phi_vec>=phi0,1);
end

if theta0<0
    theta_LDI = find(abs(theta_vec)<=abs(theta0),1);
else
    theta_LDI = find(theta_vec>=theta0,1);
end

% look direction unit vector
u0 = [cosd(phi0)*cosd(theta0);
      sind(phi0)*cosd(theta0);
      sind(theta0)];
  
% time delay for look direction
tau0 = p.AEL'*u0/p.c;

% hydrophone main response axis
H_MRA = exp(-1i*2*pi*f0*tau0)/p.Nsensors;

% moving target in space
H_TgtSteer = exp(-1i*2*pi*f0/p.c*p.AEL'*u);

% hydrophone beampattern
BP.H = H_MRA'*H_TgtSteer; BP.H = reshape(BP.H,l_phi,l_theta);
% hydrophone beampattern slices
BP.Haz = BP.H(:,theta_LDI); BP.Hde = BP.H(phi_LDI,:);

% x-component of pressure main response axis
Px_MRA = H_MRA*u0(1);
% y-component of pressure main response axis
Py_MRA = H_MRA*u0(2);
% z-component of pressure main response axis
Pz_MRA = H_MRA*u0(3);

% combined hydrophone and xyz MRAs
VS_MRA = [H_MRA; Px_MRA; Py_MRA; Pz_MRA]/sqrt(p.Nsensors);

% moving target in space for xyz channels
VSx_TgtSteer = bsxfun(@times,H_TgtSteer,ux.');
VSy_TgtSteer = bsxfun(@times,H_TgtSteer,uy.');
VSz_TgtSteer = bsxfun(@times,H_TgtSteer,uz.');

% vector sensor beampattern
BP.VS = VS_MRA'*[H_TgtSteer; VSx_TgtSteer; VSy_TgtSteer; VSz_TgtSteer];
BP.VS = reshape(BP.VS,l_phi,l_theta);
% vector sensor beampattern slices
BP.VSaz = BP.VS(:,theta_LDI); BP.VSde = BP.VS(phi_LDI,:);


% figure('position',[1950 10 1600 800])
% subplot(2,2,[2 4])
% imagesc(phi_vec,theta_vec,20*log10(abs(BP.H')))
% set(gca,'Ydir','normal')
% caxis([-40 0]); xlabel('Port$\leftarrow\phi\rightarrow$Stbd'); ylabel('Down$\leftarrow\theta\rightarrow$Up'); title(['2D H Beampattern, $f_0 = $' num2str(f0) ' Hz'])
% hcb = colorbar; xlabel(hcb, 'dB')
% 
% subplot(2,2,1)
% plot(theta_vec,10*log10(abs(BP.Hde)))
% xlim([theta_vec(1) theta_vec(end)]);
% xlabel('Down$\leftarrow\theta\rightarrow$Up'); ylabel('[dB]');
% title(['D/E H Beampattern along $\phi = $' num2str(phi0) '$^{\circ}$']);
% 
% subplot(2,2,3)
% plot(phi_vec,10*log10(abs(BP.Haz)))
% xlim([phi_vec(1) phi_vec(end)]); ylim([-40 0])
% xlabel('Port$\leftarrow\phi\rightarrow$Stbd'); ylabel('[dB]');
% title(['AZ H Beampattern along $\theta = $' num2str(round(theta0*10)/10) '$^{\circ}$']);


% figure('position',[1950 10 1600 800])
% subplot(2,2,[2 4])
% imagesc(phi_vec,theta_vec,20*log10(abs(BP.VS')))
% set(gca,'Ydir','normal')
% caxis([-40 0]); xlabel('Port$\leftarrow\phi\rightarrow$Stbd'); ylabel('Down$\leftarrow\theta\rightarrow$Up');
% title(['2D VS Beampattern, $f_0 = $' num2str(f0) ' Hz'])
% hcb = colorbar; xlabel(hcb, 'dB')
% 
% subplot(2,2,1)
% plot(theta_vec,10*log10(abs(BP.VSde)))
% xlim([theta_vec(1) theta_vec(end)]);
% xlabel('Down$\leftarrow\theta\rightarrow$Up'); ylabel('[dB]');
% title(['D/E VS Beampattern along $\phi = $' num2str(phi0) '$^{\circ}$']);
% 
% subplot(2,2,3)
% plot(phi_vec,10*log10(abs(BP.VSaz)))
% xlim([phi_vec(1) phi_vec(end)]); ylim([-40 0])
% xlabel('Port$\leftarrow\phi\rightarrow$Stbd'); ylabel('[dB]');
% title(['AZ VS Beampattern along $\theta = $' num2str(round(theta0*10)/10) '$^{\circ}$']);

end

