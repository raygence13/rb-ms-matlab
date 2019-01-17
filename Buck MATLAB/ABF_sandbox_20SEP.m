clearvars
close all
clc
present(0)
profile on

%% Need to create center-symmetric CSA
% requires M and N odd
M = 4;
N = 5;
beta = 3;

subA = zeros(floor((beta*M*N-1)/2),1);
subA(1:N:end) = 1;
M_array = [flipud(subA);subA(2:end)];

subB = zeros(floor((beta*M*N-1)/2),1);
subB(1:M:end) = 1;
N_array = [flipud(subB);subB(2:end)];

CSA = double(or(M_array,N_array));
sensorind = ((-(length(CSA)-1)/2)):(((length(CSA)-1)/2));

L = length(CSA); % Number of Sensors in ULA
figure(1)
subplot(211)
% stem(sensorind,CSA,'k')
stem(sensorind,M_array,'b')
hold on
stem(sensorind,N_array,'color',[0 .5 0])
xlabel('Sensor Index')
SampleSize = 2;


%% Inputs
% Creating Signal of Interest (SoI)
var1    = 1;    % variance of SoI
u1      = 0;    % SoI direction of arrival (DoA)
v1      = exp(1i*pi*(0:L-1)'*u1);  % SoI replica vector
b1      = sqrt(var1/2)*randn(1,SampleSize)...   % SoI Samples, [1 x SampleSize]
          + 1i*sqrt(var1/2)*randn(1,SampleSize);
SoI_Data    = v1*b1;     % SoI Data, [L x SampleSize]
ECM_SoI     = var1*(v1*v1');    % Ensemble Covariance Matrix (ECM) for SoI
SCM_SoI     = (SoI_Data*SoI_Data')/SampleSize;  % Sample Covariance Matrix (SCM) for SoI

% Creating Interferer (Int)
var2    = var1*100; % variance of Int
u2      = 3/L;      % Int DoA (in CBF peak side lobe)
v2      = exp(1i*pi*(0:L-1)'*u2);  % Int replica vector
b2      = sqrt(var2/2)*randn(1,SampleSize)...   % Int Samples, [1 x SampleSize]
          + 1i*sqrt(var2/2)*randn(1,SampleSize);
Int_Data    = v2*b2;     % Int Data, [L x SampleSize]
ECM_Int     = var2*(v2*v2');
SCM_Int     = (Int_Data*Int_Data')/SampleSize;  % Int SCM

% Creating uncorrelated Gaussian Noise
varW        = var1;
Noise_Data  = sqrt(varW/2)*randn(L,SampleSize)...
              + 1i*sqrt(varW/2)*randn(L,SampleSize);
ECM_Noise   = varW*eye(L);
SCM_Noise   = (Noise_Data*Noise_Data')/SampleSize;

x = SoI_Data + Int_Data + Noise_Data;   % Input vector
xM = bsxfun(@times,x,M_array);
xN = bsxfun(@times,x,N_array);
ECM = ECM_SoI + ECM_Int + ECM_Noise;
S_x = (x*x')/SampleSize;      % Data SCM
S_csa = (xM*xN')/SampleSize;

CSApp_coarray = conv(M_array,flipud(N_array));
contiguous_min_lag = floor( (beta-1)*M*N + (M+N)/2 -1);

lag = (-(L-1)):(L-1);
cont_lag_region = -contiguous_min_lag:contiguous_min_lag;
CSApp_coarray_cont = ...
    CSApp_coarray(find(lag == -contiguous_min_lag):find(lag == contiguous_min_lag));
figure(1)
subplot(212)
stem(lag,CSApp_coarray,'k')
xlabel('$\gamma$')
rCSApp = zeros(size(cont_lag_region));
for ind = 1:length(cont_lag_region)
    rCSApp(ind) = sum(diag(S_csa,cont_lag_region(ind)))/CSApp_coarray_cont(ind);
end
S_csa_aug = toeplitz(rCSApp((length(rCSApp)+1)/2:-1:1),rCSApp((length(rCSApp)+1)/2:end));
S_csa_smoothed = S_csa_aug.^2/(length(cont_lag_region));


%% Outputs
uNfft = 2048;
ur = linspace(0,1,uNfft/2);
ul = linspace(-1,0-1/uNfft,uNfft/2);
u = [ul ur];

u0      = 0;    % Look direction
v0      = exp(1i*pi*(0:L-1)'*u0);	% manifold vector

%% MVDR
DL          = 5*varW*eye(L);  % Diagonal loading factor
EMI         = (ECM)\v0;   % ECM inverted
wMVDR_EMI   = EMI/(v0.'*EMI);               % MVDR weight vector
EBP_MVDR = fftshift(fft(wMVDR_EMI,uNfft));  % Ensemble beampattern for MVDR

% SMI         = (SCM_Int+SCM_Noise)\v0;
SMI = S_x\v0;
wMVDR_SMI   = SMI/(v0.'*SMI);
SBP_MVDR = fftshift(fft(wMVDR_SMI,uNfft));

SMIcsa = S_csa_aug\v0(1:contiguous_min_lag+1);

% SMI_DL         = (SCM_Int+SCM_Noise+DL)\v0;
SMI_DL          = (S_x + DL)\v0;
wMVDR_SMI_DL   = SMI_DL/(v0.'*SMI_DL);
SBP_DLMVDR = fftshift(fft(wMVDR_SMI_DL,uNfft));

%% DMR
[EnsembleEigVec,EnsembleEigVal] = eig(ECM);
[gamma,I] = sort(diag(EnsembleEigVal),1,'descend');
Xi = EnsembleEigVec(:,I);

Xi_dominant = gamma(1)*Xi(:,1)*Xi(:,1)';
Xi_noise = zeros(L,L,L-1);
for noise_ind = 1:L-1
    Xi_noise(:,:,noise_ind) = varW*Xi(:,1+noise_ind)*Xi(:,1+noise_ind)';
end
Xi_noise = sum(Xi_noise,3);

Sigma_DMR = Xi_dominant + Xi_noise;
Ens_wDMR = (Sigma_DMR\v0)/(v0.'*(Sigma_DMR\v0));
Ens_BP_DMR = fftshift(fft(Ens_wDMR,uNfft));

% [SampleEigVec,SampleEigVal] = eig(S_csa_smoothed);
% [g,I] = sort(diag(SampleEigVal),1,'descend');
% E = SampleEigVec(:,I);
% 
% g_norm = 10*log10(g/max(g));
% D = sum(g_norm > -10);
% varW_hat = sum(g_norm(D+1:end))*(SampleSize/(SampleSize-1))*(1/(L-D));
% 
% E_dominant = zeros(contiguous_min_lag+1,contiguous_min_lag+1,D);
% for int_ind = 1:D
%     E_dominant(:,:,int_ind) = g(int_ind)*E(:,int_ind)*E(:,int_ind)';
% end
% E_dominant = sum(E_dominant,3);
% 
% E_noise = zeros(contiguous_min_lag+1,contiguous_min_lag+1,contiguous_min_lag-D);
% for noise_ind = 1:contiguous_min_lag-D
%     E_noise(:,:,noise_ind) = varW_hat*E(:,D+noise_ind)*E(:,D+noise_ind)';
% end
% E_noise = sum(E_noise,3);
% 
% S_DMR = E_dominant + E_noise;
% 
% wDMR = (S_DMR\v0(1:length(S_DMR)))/(v0(1:length(S_DMR)).'*(S_DMR\v0(1:length(S_DMR))));
% BP_DMR = fftshift(fft(wDMR,uNfft));

% ULA DMR Beamforming
[SampleEigVec,SampleEigVal] = eig(S_x);
[g,I] = sort(diag(SampleEigVal),1,'descend');
E = SampleEigVec(:,I);

g_norm = 10*log10(g/max(g));
D = sum(g_norm > -10);
varW_hat = sum(g_norm(D+1:end))*(SampleSize/(SampleSize-1))*(1/(L-D));

E_dominant = zeros(L,L,D);
for int_ind = 1:D
    E_dominant(:,:,int_ind) = g(int_ind)*E(:,int_ind)*E(:,int_ind)';
end
E_dominant = sum(E_dominant,3);

E_noise = zeros(L,L,L-D);
for noise_ind = 1:L-D
    E_noise(:,:,noise_ind) = varW_hat*E(:,D+noise_ind)*E(:,D+noise_ind)';
end
E_noise = sum(E_noise,3);

S_DMR = E_dominant + E_noise;

wDMR = (S_DMR\v0)/(v0.'*(S_DMR\v0));
BP_DMR = fftshift(fft(wDMR,uNfft));
%%
figure(2)
subplot(211)
stem(u2,10*log10(var2));
hold on
stem(u1,10*log10(var1))
xlim([-1 1])

subplot(212)
EMVDR_plot = plot(u,10*log10(abs(EBP_MVDR)));
hold on
SMVDR_plot = plot(u,10*log10(abs(SBP_MVDR)),'r');
DLMVDR_plot = plot(u,10*log10(abs(SBP_DLMVDR)),'k');
EDMR_plot = plot(u,10*log10(abs(Ens_BP_DMR)),'g--');
SDMR_plot = plot(u,10*log10(abs(BP_DMR)),'m');
ylim([-60 30])
legend([EMVDR_plot,SMVDR_plot, DLMVDR_plot,EDMR_plot, SDMR_plot],'Ensemble MVDR','SMI MVDR','DL MVDR','Ensemble DMR','Sample DMR')
