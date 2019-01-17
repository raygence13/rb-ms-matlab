%% MATLAB script to be concurrent with M.S. Thesis
% R. Bautista
% 1 May 2017
clearvars
close all
clc
present(0)
profile on

%%
% Color code for each undersampling factor
plotcolor.ULA       = [1 0 0];          % ULA: red
plotcolor.M         = [0 .25 .9];       % M: blue
plotcolor.N         = [0 0.6 0.2];      % N: green
plotcolor.CSAcbf    = [1 .6 0];         % CSAcbf: orange
plotcolor.CSApp     = [0 0 0];          % CSApp: purple
plotcolor.sim       = [.25 .25 .25];    % sim: almost black

DEBUG = 1;
doROC = 0;
islocal = 1;

% Correlation length values
% alpha = [0 .1 .25 .3 .5 .6 .75 .85 .9 .99];
% alpha = 0:.01:.99;
alpha = 0;
     
%% Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Initialize input/output noise powers as a function of alpha
%%%%    for memory purposes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dummy_to_collapse_variable_initializiation = 1:1 
% Input Signal Powers
simulated.ULA.InSigPow      = zeros(1,length(alpha));
simulated.M.InSigPow        = zeros(1,length(alpha));
simulated.N.InSigPow        = zeros(1,length(alpha));
simulated.CSAcbf.InSigPow   = zeros(1,length(alpha));
simulated.CSApp.InSigPow    = zeros(1,length(alpha));

% Input Noise Powers
derived.InNoisePow          = zeros(1,length(alpha));
simulated.ULA.InNoisePow    = zeros(1,length(alpha));
simulated.M.InNoisePow      = zeros(1,length(alpha));
simulated.N.InNoisePow      = zeros(1,length(alpha));
simulated.CSAcbf.InNoisePow = zeros(1,length(alpha));
simulated.CSApp.InNoisePow  = zeros(1,length(alpha));

% Input SNRs:
derived.InSNR           = zeros(1,length(alpha));
simulated.ULA.InSNR     = zeros(1,length(alpha));
simulated.M.InSNR       = zeros(1,length(alpha));
simulated.N.InSNR       = zeros(1,length(alpha));
simulated.CSAcbf.InSNR  = zeros(1,length(alpha));
simulated.CSApp.InSNR   = zeros(1,length(alpha));

% Output Signal Powers:
simulated.ULA.OutSigPow     = zeros(1,length(alpha));
simulated.M.OutSigPow       = zeros(1,length(alpha));
simulated.N.OutSigPow       = zeros(1,length(alpha));
simulated.CSAcbf.OutSigPow  = zeros(1,length(alpha));
simulated.CSApp.OutSigPow   = zeros(1,length(alpha));

% Output Noise Powers:
derived.ULA.OutNoisePow     = zeros(1,length(alpha));
simulated.ULA.OutNoisePow   = zeros(1,length(alpha));

derived.M.OutNoisePow      = zeros(1,length(alpha));
simulated.M.OutNoisePow    = zeros(1,length(alpha));

derived.N.OutNoisePow      = zeros(1,length(alpha));
simulated.N.OutNoisePow    = zeros(1,length(alpha));

derived.CSAcbf.OutNoisePow      = zeros(1,length(alpha));
simulated.CSAcbf.OutNoisePow    = zeros(1,length(alpha));

derived.CSApp.OutNoisePow = zeros(1,length(alpha));
simulated.CSApp.OutNoisePow = zeros(1,length(alpha));

% Output Signal + Noise Powers:
derived.ULA.OutDataSigNoisePow = zeros(1,length(alpha));
ULA_OutDataSigNoisePow = zeros(1,length(alpha));

derived.M.OutDataSigNoisePow = zeros(1,length(alpha));
simulated.M.OutDataSigNoisePow = zeros(1,length(alpha));

derived.N.OutDataSigNoisePow = zeros(1,length(alpha));
simulated.N.OutDataSigNoisePow = zeros(1,length(alpha));

derived.CSAcbf.OutDataSigNoisePow = zeros(1,length(alpha));
simulated.CSAcbf.OutDataSigNoisePow = zeros(1,length(alpha));

derived.CSApp.OutDataSigNoisePow = zeros(1,length(alpha));
simulated.CSApp.OutDataSigNoisePow = zeros(1,length(alpha));

% Output SNRs:
derived.ULA.OutSNR = zeros(1,length(alpha));
simulated.ULA.OutSNR = zeros(1,length(alpha));
simulated.ULA.PowOutNoisePow  = zeros(1,length(alpha));

derived.M.OutSNR = zeros(1,length(alpha));
simulated.M.OutSNR = zeros(1,length(alpha));
simulated.M.PowOutNoisePow  = zeros(1,length(alpha));

derived.N.OutSNR = zeros(1,length(alpha));
simulated.N.OutSNR = zeros(1,length(alpha));
simulated.N.PowOutNoisePow  = zeros(1,length(alpha));

derived.CSAcbf.OutSNR = zeros(1,length(alpha));
simulated.CSAcbf.OutSNR = zeros(1,length(alpha));

derived.CSApp.PTS             = cell(1,1);
derived.CSApp.PowOutNoisePow  = zeros(1,length(alpha));
simulated.CSApp.PowOutNoisePow  = zeros(1,length(alpha));
derived.CSApp.OutSNR          = zeros(1,length(alpha));
simulated.CSApp.OutSNR          = zeros(1,length(alpha));

% Noise Only PDFs and CDFs
derived.ULA.NoisePDF     = cell(1,length(alpha));
derived.ULA.xNoisePDF    = cell(1,length(alpha));
simulated.ULA.NoisePDF     = cell(1,length(alpha));
simulated.ULA.xNoisePDF    = cell(1,length(alpha));
derived.ULA.NoiseCDF  = cell(1,length(alpha));
derived.ULA.xNoiseCDF = cell(1,length(alpha));
simulated.ULA.NoiseCDF  = cell(1,length(alpha));
simulated.ULA.xNoiseCDF = cell(1,length(alpha));

derived.M.NoisePDF   = cell(1,length(alpha));
derived.M.xNoisePDF  = cell(1,length(alpha));
simulated.M.NoisePDF   = cell(1,length(alpha));
simulated.M.xNoisePDF  = cell(1,length(alpha));
derived.M.NoiseCDF  = cell(1,length(alpha));
derived.M.xNoiseCDF = cell(1,length(alpha));
simulated.M.NoiseCDF  = cell(1,length(alpha));
simulated.M.xNoiseCDF = cell(1,length(alpha));

derived.N.NoisePDF  = cell(1,length(alpha));
derived.N.xNoisePDF = cell(1,length(alpha));
simulated.N.NoisePDF  = cell(1,length(alpha));
simulated.N.xNoisePDF = cell(1,length(alpha));
derived.N.NoiseCDF  = cell(1,length(alpha));
derived.N.xNoiseCDF = cell(1,length(alpha));
simulated.N.NoiseCDF  = cell(1,length(alpha));
simulated.N.xNoiseCDF = cell(1,length(alpha));

derived.CSAcbf.NoisePDF  = cell(1,length(alpha));
derived.CSAcbf.xNoisePDF = cell(1,length(alpha));
simulated.CSAcbf.NoisePDF  = cell(1,length(alpha));
simulated.CSAcbf.xNoisePDF = cell(1,length(alpha));
derived.CSAcbf.NoiseCDF  = cell(1,length(alpha));
derived.CSAcbf.xNoiseCDF = cell(1,length(alpha));
simulated.CSAcbf.NoiseCDF  = cell(1,length(alpha));
simulated.CSAcbf.xNoiseCDF = cell(1,length(alpha));

derived.CSApp.rhoNoise = zeros(1,length(alpha));
simulated.CSApp.rhoNoise = zeros(1,length(alpha));
derived.CSApp.NoisePDF  = cell(1,length(alpha));
derived.CSApp.xNoisePDF = cell(1,length(alpha));
simulated.CSApp.NoisePDF  = cell(1,length(alpha));
simulated.CSApp.xNoisePDF = cell(1,length(alpha));
derived.CSApp.NoiseCDF  = cell(1,length(alpha));
derived.CSApp.xNoiseCDF = cell(1,length(alpha));
simulated.CSApp.NoiseCDF  = cell(1,length(alpha));
simulated.CSApp.xNoiseCDF = cell(1,length(alpha));

% Signal + Noise PDFs and CDFs
derived.ULA.SigNoisePDF     = cell(1,length(alpha));
derived.ULA.xSigNoisePDF    = cell(1,length(alpha));
simulated.ULA.SigNoisePDF     = cell(1,length(alpha));
simulated.ULA.xSigNoisePDF    = cell(1,length(alpha));
derived.ULA.SigNoiseCDF  = cell(1,length(alpha));
derived.ULA.xSigNoiseCDF = cell(1,length(alpha));
simulated.ULA.SigNoiseCDF  = cell(1,length(alpha));
simulated.ULA.xSigNoiseCDF = cell(1,length(alpha));

derived.ULA.Pfa = cell(1,length(alpha));
derived.ULA.Pd = cell(1,length(alpha));

derived.M.SigNoisePDF   = cell(1,length(alpha));
derived.M.xSigNoisePDF  = cell(1,length(alpha));
simulated.M.SigNoisePDF   = cell(1,length(alpha));
simulated.M.xSigNoisePDF  = cell(1,length(alpha));
derived.M.SigNoiseCDF  = cell(1,length(alpha));
derived.M.xSigNoiseCDF = cell(1,length(alpha));
simulated.M.SigNoiseCDF  = cell(1,length(alpha));
simulated.M.xSigNoiseCDF = cell(1,length(alpha));

derived.N.SigNoisePDF  = cell(1,length(alpha));
derived.N.xSigNoisePDF = cell(1,length(alpha));
simulated.N.SigNoisePDF  = cell(1,length(alpha));
simulated.N.xSigNoisePDF = cell(1,length(alpha));
derived.N.SigNoiseCDF  = cell(1,length(alpha));
derived.N.xSigNoiseCDF = cell(1,length(alpha));
simulated.N.SigNoiseCDF  = cell(1,length(alpha));
simulated.N.xSigNoiseCDF = cell(1,length(alpha));

derived.CSAcbf.SigNoisePDF  = cell(1,length(alpha));
derived.CSAcbf.xSigNoisePDF = cell(1,length(alpha));
simulated.CSAcbf.SigNoisePDF  = cell(1,length(alpha));
simulated.CSAcbf.xSigNoisePDF = cell(1,length(alpha));
derived.CSAcbf.SigNoiseCDF  = cell(1,length(alpha));
derived.CSAcbf.xSigNoiseCDF = cell(1,length(alpha));
simulated.CSAcbf.SigNoiseCDF  = cell(1,length(alpha));
simulated.CSAcbf.xSigNoiseCDF = cell(1,length(alpha));

derived.CSAcbf.Pfa = cell(1,length(alpha));
derived.CSAcbf.Pd = cell(1,length(alpha));

derived.CSApp.rhoSigNoise = zeros(1,length(alpha));
simulated.CSApp.rhoSigNoise = zeros(1,length(alpha));
derived.CSApp.SigNoisePDF  = cell(1,length(alpha));
derived.CSApp.xSigNoisePDF = cell(1,length(alpha));
simulated.CSApp.SigNoisePDF  = cell(1,length(alpha));
simulated.CSApp.xSigNoisePDF = cell(1,length(alpha));
derived.CSApp.SigNoiseCDF  = cell(1,length(alpha));
derived.CSApp.xSigNoiseCDF = cell(1,length(alpha));
simulated.CSApp.SigNoiseCDF  = cell(1,length(alpha));
simulated.CSApp.xSigNoiseCDF = cell(1,length(alpha));

derived.CSApp.Pfa = cell(1,length(alpha));
derived.CSApp.Pd = cell(1,length(alpha));

% Sig + Noise & Noise Only KSTest
derived.ULA.NoiseKSTest     = zeros(1,length(alpha));
derived.M.NoiseKSTest  = zeros(1,length(alpha));
derived.N.NoiseKSTest  = zeros(1,length(alpha));
derived.CSAcbf.NoiseKSTest  = zeros(1,length(alpha));
derived.CSApp.NoiseKSTest   = zeros(1,length(alpha));

derived.ULA.SigNoiseKSTest      = zeros(1,length(alpha));
derived.M.SigNoiseKSTest   = zeros(1,length(alpha));
derived.N.SigNoiseKSTest   = zeros(1,length(alpha));
derived.CSAcbf.SigNoiseKSTest   = zeros(1,length(alpha));
derived.CSApp.SigNoiseKSTest    = zeros(1,length(alpha));

% Biased Noise PSDs
truePSDNoise                = cell(1,length(alpha));
derived.ULA.NoisePSD        = cell(1,length(alpha));
simulated.ULA.NoisePSD      = cell(1,length(alpha));
derived.M.NoisePSD     = cell(1,length(alpha));
simulated.M.NoisePSD   = cell(1,length(alpha));
derived.N.NoisePSD     = cell(1,length(alpha));
simulated.N.NoisePSD   = cell(1,length(alpha));
derived.CSAcbf.NoisePSD     = cell(1,length(alpha));
simulated.CSAcbf.NoisePSD   = cell(1,length(alpha));
derived.CSApp.NoisePSD      = cell(1,length(alpha));
simulated.CSApp.NoisePSD    = cell(1,length(alpha));

% Biased Signal + Noise PSDs
truePSDSigNoise                 = cell(1,length(alpha));
derived.ULA.SigNoisePSD     	= cell(1,length(alpha));
simulated.ULA.SigNoisePSD       = cell(1,length(alpha));
derived.M.SigNoisePSD      = cell(1,length(alpha));
simulated.M.SigNoisePSD    = cell(1,length(alpha));
derived.N.SigNoisePSD      = cell(1,length(alpha));
simulated.N.SigNoisePSD    = cell(1,length(alpha));
derived.CSAcbf.SigNoisePSD      = cell(1,length(alpha));
simulated.CSAcbf.SigNoisePSD    = cell(1,length(alpha));
derived.CSApp.SigNoisePSD       = cell(1,length(alpha));
simulated.CSApp.SigNoisePSD     = cell(1,length(alpha));


derived.ULA.Beampattern     = cell(1,1);
derived.M.Beampattern  = cell(1,1);
derived.N.Beampattern  = cell(1,1);
derived.CSAcbf.Beampattern  = cell(1,1);
derived.CSApp.Beampattern   = cell(1,1);
end

%% Array Parameters

% undersampling factors
M = 2;              % M-undersampling
N = 3;              % N-undersampling
beta = 10;          % extension factor/number of shared sensors
L = M*N*beta;       % Number of sensors in full ULA aperture
l = 0:L-1;          % sensor index

% Spatial frequency u = sin(theta)
uNfft = 2048;
ur = linspace(0,1,uNfft/2);
ul = linspace(-1,0-1/uNfft,uNfft/2);
u = [ul ur];

% sensor position vector for ULA
derived.ULA.ArrayPos = ones(1,L);

% sensor position vectors for subarrays
derived.M.ArrayPos = zeros(1,L);
derived.M.ArrayPos(1:M:end) = 1;
I_M = sum(derived.M.ArrayPos);  % number of sensors in M-subarray

derived.N.ArrayPos = zeros(1,L);
derived.N.ArrayPos(1:N:end) = 1;
I_N = sum(derived.N.ArrayPos);  % number of sensors in N-subarray

% combined subarrays/coprime array
derived.CSAcbf.ArrayPos = double(or(derived.M.ArrayPos,derived.N.ArrayPos));
I_csa = sum(derived.CSAcbf.ArrayPos); % number of unique sensors in CSA
beta = sum(and(derived.M.ArrayPos,derived.N.ArrayPos));

% auto-correlation of ULA array sensor position Wieght Function
derived.ULA.coarray = conv(derived.ULA.ArrayPos,fliplr(derived.ULA.ArrayPos));
% auto-correlation of M sensor position Wieght Function
derived.M.coarray = conv(derived.M.ArrayPos,fliplr(derived.M.ArrayPos));
% auto-correlation of M sensor position Wieght Function
derived.N.coarray = conv(derived.N.ArrayPos,fliplr(derived.N.ArrayPos));
% auto-correlation of CSA array sensor position Wieght Function
derived.CSAcbf.coarray = conv(derived.CSAcbf.ArrayPos,fliplr(derived.CSAcbf.ArrayPos));
% cross-correlation of M and N sensor position Wieght Function
derived.CSApp.coarray = conv(derived.M.ArrayPos,fliplr(derived.N.ArrayPos));

% lag index for coarray
gamma = -(L-1):(L-1);

% signal direction
uS = 0;
vS = exp(-1i*pi*(0:L-1)'*uS);
uI = 3/L;
vI = exp(-1i*pi*(0:L-1)'*uI);
% replica vectors for each geometry: vULA,vM,vN,vCSA
% vULA: [L x 1]
% vM: [L x 1]
vM = vS.*derived.M.ArrayPos.';
% vN: [L x 1]
vN = vS.*derived.N.ArrayPos.';
% vCSA: zero-filled [L x 1]
vCSA = vS.*derived.CSAcbf.ArrayPos.';

% weight vectors for each geometry: wULA,wM,wN,wCSA
% wULA: [L x 1]
wULA = 1/L*vS;
% wM: [Me x 1]
wM = 1/I_M*vM;
% wN: [Ne x 1]
wN = 1/I_N*vN;
% wCSA: zero-filled [L x 1]
wCSA = 1/I_csa*vCSA;

varS = 1;   % signal power
% Analytical Signal Only Covariance Matrix for each processor
RULA_Sig = vS*vS'*varS;
rULA_Sig = zeros(size(gamma));

RM_Sig = vM*vM'*varS;
rM_Sig = zeros(size(gamma));

RN_Sig = vN*vN'*varS;
rN_Sig = zeros(size(gamma));

RCSAcbf_Sig = vCSA*vCSA'*varS;
rCSAcbf_Sig = zeros(size(gamma));

RCSApp_Sig = vM*vN'*varS;
rCSApp_Sig = zeros(size(gamma));

for ind = 1:length(gamma)
    rULA_Sig(ind) = sum(diag(RULA_Sig,gamma(ind)));
    rM_Sig(ind) = sum(diag(RM_Sig,gamma(ind)));
    rN_Sig(ind) = sum(diag(RN_Sig,gamma(ind)));
    rCSAcbf_Sig(ind) = sum(diag(RCSAcbf_Sig,gamma(ind)));
    rCSApp_Sig(ind) = sum(diag(RCSApp_Sig,gamma(ind)));
end

% Analytical Signal Only Correlation Function for each processor
rULA_Sig = rULA_Sig/L;
rM_Sig = rM_Sig/I_M;
rN_Sig = rN_Sig/I_N;
rCSAcbf_Sig = rCSAcbf_Sig/I_csa;
rCSApp_Sig = rCSApp_Sig/beta;

SigPSD = fftshift(fft(rULA_Sig,uNfft));

derived.ULA.Beampattern = fftshift(fft(rULA_Sig/L/varS,uNfft));
derived.M.Beampattern = fftshift(fft(rM_Sig/L*M/varS,uNfft));
derived.N.Beampattern = fftshift(fft(rN_Sig/L*N/varS,uNfft));
derived.CSAcbf.Beampattern = fftshift(fft(rCSAcbf_Sig/I_csa/varS,uNfft));
derived.CSApp.Beampattern = fftshift(fft(rCSApp_Sig/L/varS,uNfft));

% creating Power Terms Sequence (PTS) in alpha^{|aM-bN|+|dM-eN|} and alpha[M|a-d|+N|b-e|]
% will ultimately be used as alpha^derived.CSApp.PTS
% creating Steer Direction Terms (SDT) in
% exp(-j\pi(aM-bN)u0)exp(j\pi(dM-eN)u0)
% exp(-j\pi(aM-dM)u0)exp(j\pi(bN-eN)u0)
[derived.CSApp.PTS,SDT] = CSAProdPowTermsSeq(L,I_M,I_N,M,N,beta,uS);

for a_ind = 1:length(alpha)
    %% Alpha portion
    threshold = -15;     % Sets value for extrasamples/sensors for steady state
    [b,a,extrasamples] = ARProcessCoefficients( alpha(a_ind),1,threshold );
    
    SampleSize = 2*L+3;  % number of snapshots
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      creating noise
    %%%%    simulated.*.InDataNoise contains noisy data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % white noise filter input
    varW = 1;           % white noise variance
    % Spatially uncorrelated noise n[l]
    n = sqrt(varW/2)*randn(L+extrasamples,SampleSize)...
        + 1i*sqrt(varW/2)*randn(L+extrasamples,SampleSize);
    
    % Spatially correlated noise c[l] = sqrt{1-\alpha^2}*n[l] + \alpha*c[l-1]
    c_os = filter(b*sqrt(1-alpha(a_ind)^2),a,n);    % Oversampled process [L+extrasamples x SampleSize]
    
    ULA_InDataNoise = c_os(end-L+1:end,:);                  % Last L samples [L x SampleSize]
    M_InDataNoise = bsxfun(@times,ULA_InDataNoise,derived.M.ArrayPos.');
    N_InDataNoise = bsxfun(@times,ULA_InDataNoise,derived.N.ArrayPos.');
    CSAcbf_InDataNoise = bsxfun(@times,ULA_InDataNoise,derived.CSAcbf.ArrayPos.');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Noise Power Spectral Density Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    truePSDNoise{1,a_ind} = (b*sqrt(1-alpha(a_ind)^2)./(a(1) + a(2).*exp(-1i*pi.*u))).^2*varW;
    
    % Expected AR Noise Correlation Estimates
    rULA_Noise = derived.ULA.coarray.*alpha(a_ind).^abs(gamma)/L*varW;
    derived.ULA.NoisePSD{1,a_ind} = fftshift(fft(rULA_Noise,uNfft,2));
    
    rM_Noise = derived.M.coarray.*alpha(a_ind).^abs(gamma)*M/L*varW;
    derived.M.NoisePSD{1,a_ind} = fftshift(fft(rM_Noise,uNfft,2));
    
    rN_Noise = derived.N.coarray.*alpha(a_ind).^abs(gamma)*N/L*varW;
    derived.N.NoisePSD{1,a_ind} = fftshift(fft(rN_Noise,uNfft,2));
    
    rCSAcbf_Noise = derived.CSAcbf.coarray.*alpha(a_ind).^abs(gamma)/I_csa*varW;
    derived.CSAcbf.NoisePSD{1,a_ind} = fftshift(fft(rCSAcbf_Noise,uNfft,2));
    
    rCSApp_Noise = derived.CSApp.coarray.*alpha(a_ind).^abs(gamma)/beta*varW;
    derived.CSApp.NoisePSD{1,a_ind} = fftshift(fft(rCSApp_Noise,uNfft,2));
    
    % Noise Only Covariance Matrix for each processor
    RULA_simNoise       = ULA_InDataNoise*ULA_InDataNoise'/SampleSize;
    RM_simNoise         = M_InDataNoise*M_InDataNoise'/SampleSize;
    RN_simNoise         = N_InDataNoise*N_InDataNoise'/SampleSize;
    RCSAcbf_simNoise    = CSAcbf_InDataNoise*CSAcbf_InDataNoise'/SampleSize;
    RCSApp_simNoise     = M_InDataNoise*N_InDataNoise'/SampleSize;
    
    rULA_simNoise    = zeros(size(gamma));
    rM_simNoise      = zeros(size(gamma));
    rN_simNoise      = zeros(size(gamma));
    rCSAcbf_simNoise = zeros(size(gamma));
    rCSApp_simNoise  = zeros(size(gamma));
    for ind = 1:length(gamma)
        rULA_simNoise(ind) = sum(diag(RULA_simNoise,gamma(ind)));
        rM_simNoise(ind) = sum(diag(RM_simNoise,gamma(ind)));
        rN_simNoise(ind) = sum(diag(RN_simNoise,gamma(ind)));
        rCSAcbf_simNoise(ind) = sum(diag(RCSAcbf_simNoise,gamma(ind)));
        rCSApp_simNoise(ind) = sum(diag(RCSApp_simNoise,gamma(ind)));
    end
    
    % Simulated AR Noise Only Correlation Expected Estimates
    rULA_simNoise = rULA_simNoise./L;
    rM_simNoise = rM_simNoise*M/L;
    rN_simNoise = rN_simNoise*N/L;
    rCSAcbf_simNoise = rCSAcbf_simNoise/I_csa;
    rCSApp_simNoise = rCSApp_simNoise/beta;
    
    simulated.ULA.NoisePSD{1,a_ind} = fftshift(fft(rULA_simNoise,uNfft));
    simulated.M.NoisePSD{1,a_ind} = fftshift(fft(rM_simNoise,uNfft));
    simulated.N.NoisePSD{1,a_ind} = fftshift(fft(rN_simNoise,uNfft));
    simulated.CSAcbf.NoisePSD{1,a_ind} = fftshift(fft(rCSAcbf_simNoise,uNfft));
    simulated.CSApp.NoisePSD{1,a_ind} = fftshift(fft(rCSApp_simNoise,uNfft));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      creating signal
    %%%%    simulated.*.InDataSig contains signal data in freq. domain at each sensor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = sqrt(varS/2)*randn(1,SampleSize)...
        + 1i*sqrt(varS/2)*randn(1,SampleSize);  % signal, [1 x SampleSize]
    
    varI = varS*100;
    interferer = sqrt(varI/2)*randn(1,SampleSize)...
        + 1i*sqrt(varI/2)*randn(1,SampleSize);  % signal, [1 x SampleSize]
    
    ULA_InDataSig = vS*s + vI*interferer;     % xULAsig, [L x SampleSize]
    M_InDataSig = vM*s + vI.*derived.M.ArrayPos.'*interferer;
    N_InDataSig = vN*s + vI.*derived.N.ArrayPos.'*interferer;
    CSAcbf_InDataSig = vCSA*s + vI.*derived.CSAcbf.ArrayPos.'*interferer;
    
    % Simulated Signal Only Covariance Matrix
    RULA_simSig       = ULA_InDataSig*ULA_InDataSig'/SampleSize;
    RM_simSig         = M_InDataSig*M_InDataSig'/SampleSize;
    RN_simSig         = N_InDataSig*N_InDataSig'/SampleSize;
    RCSAcbf_simSig    = CSAcbf_InDataSig*CSAcbf_InDataSig'/SampleSize;
    RCSApp_simSig     = M_InDataSig*N_InDataSig'/SampleSize;
    
    rULA_simSig    = zeros(size(gamma));
    rM_simSig      = zeros(size(gamma));
    rN_simSig      = zeros(size(gamma));
    rCSAcbf_simSig = zeros(size(gamma));
    rCSApp_simSig  = zeros(size(gamma));
    for ind = 1:length(gamma)
        rULA_simSig(ind) = sum(diag(RULA_simSig,gamma(ind)));
        rM_simSig(ind) = sum(diag(RM_simSig,gamma(ind)));
        rN_simSig(ind) = sum(diag(RN_simSig,gamma(ind)));
        rCSAcbf_simSig(ind) = sum(diag(RCSAcbf_simSig,gamma(ind)));
        rCSApp_simSig(ind) = sum(diag(RCSApp_simSig,gamma(ind)));
    end
    
    % Simulated Signal Only Correlation Function (Biased)
    rULA_simSig = rULA_simSig./L;
    rM_simSig = rM_simSig*M/L;
    rN_simSig = rN_simSig*N/L;
    rCSAcbf_simSig = rCSAcbf_simSig/I_csa;
    rCSApp_simSig = rCSApp_simSig/beta;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Input Signal + Noise
    %%%%    simulated.*.InData contains array data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ULA_InDataSigNoise      = ULA_InDataSig + ULA_InDataNoise;
    M_InDataSigNoise        = M_InDataSig + M_InDataNoise;
    N_InDataSigNoise        = N_InDataSig + N_InDataNoise;
    CSAcbf_InDataSigNoise   = CSAcbf_InDataSig + CSAcbf_InDataNoise;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Signal + Noise Power Spectral Density Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    truePSDSigNoise{1,a_ind} = SigPSD + truePSDNoise{1,a_ind};
    
    % Expected Signal+AR Noise Estimates
    rULA_SigNoise = rULA_Sig + rULA_Noise;
    derived.ULA.SigNoisePSD{1,a_ind} = fftshift(fft(rULA_SigNoise,uNfft,2));
    
    rM_SigNoise = rM_Sig + rM_Noise;
    derived.M.SigNoisePSD{1,a_ind} = fftshift(fft(rM_SigNoise,uNfft,2));
    
    rN_SigNoise = rN_Sig + rN_Noise;
    derived.N.SigNoisePSD{1,a_ind} = fftshift(fft(rN_SigNoise,uNfft,2));
    
    rCSAcbf_SigNoise = rCSAcbf_Sig + rCSAcbf_Noise;
    derived.CSAcbf.SigNoisePSD{1,a_ind} = fftshift(fft(rCSAcbf_SigNoise,uNfft,2));
    
    rCSApp_SigNoise = rCSApp_Sig + rCSApp_Noise;
    derived.CSApp.SigNoisePSD{1,a_ind} = fftshift(fft(rCSApp_SigNoise,uNfft,2));
    
    % Simulated Signal+Noise Covariance Matrix
    RULA_simSigNoise    = ULA_InDataSigNoise*ULA_InDataSigNoise'/SampleSize;
    RM_simSigNoise      = M_InDataSigNoise*M_InDataSigNoise'/SampleSize;
    RN_simSigNoise      = N_InDataSigNoise*N_InDataSigNoise'/SampleSize;
    RCSAcbf_simSigNoise = CSAcbf_InDataSigNoise*CSAcbf_InDataSigNoise'/SampleSize;
    RCSApp_simSigNoise  = M_InDataSigNoise*N_InDataSigNoise'/SampleSize;
    
    rULA_simSigNoise    = zeros(size(gamma));
    rM_simSigNoise = zeros(size(gamma));
    rN_simSigNoise = zeros(size(gamma));
    rCSAcbf_simSigNoise = zeros(size(gamma));
    rCSApp_simSigNoise = zeros(size(gamma));
    for ind = 1:length(gamma)
        rULA_simSigNoise(ind) = sum(diag(RULA_simSigNoise,gamma(ind)));
        rM_simSigNoise(ind) = sum(diag(RM_simSigNoise,gamma(ind)));
        rN_simSigNoise(ind) = sum(diag(RN_simSigNoise,gamma(ind)));
        rCSAcbf_simSigNoise(ind) = sum(diag(RCSAcbf_simSigNoise,gamma(ind)));
        rCSApp_simSigNoise(ind) = sum(diag(RCSApp_simSigNoise,gamma(ind)));
    end
    
    % Simulated Signal+Noise Correlation Functions
    rULA_simSigNoise = rULA_simSigNoise./L;
    rM_simSigNoise = rM_simSigNoise/I_M;
    rN_simSigNoise = rN_simSigNoise/I_N;
    rCSAcbf_simSigNoise = rCSAcbf_simSigNoise./I_csa;
    rCSApp_simSigNoise = rCSApp_simSigNoise/beta;
    
    simulated.ULA.SigNoisePSD{1,a_ind} = fftshift(fft(rULA_simSigNoise,uNfft));
    simulated.M.SigNoisePSD{1,a_ind} = fftshift(fft(rM_simSigNoise,uNfft));
    simulated.N.SigNoisePSD{1,a_ind} = fftshift(fft(rN_simSigNoise,uNfft));
    simulated.CSAcbf.SigNoisePSD{1,a_ind} = fftshift(fft(rCSAcbf_simSigNoise,uNfft));
    simulated.CSApp.SigNoisePSD{1,a_ind} = fftshift(fft(rCSApp_simSigNoise,uNfft));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Input Signal Power: E{s^2[l]} = varS
    %%%%    *.InSigPow is average power of signal
    %%%%    averaged over all trials, averaged over all sensors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simulated.ULA.InSigPow(1,a_ind) = trace(RULA_simSig);
    % E{(1/L)sum(rULA_simSig[gamma])}
    
    simulated.M.InSigPow(1,a_ind) = trace(RM_simSig);
    % E{(1/Me)sum(xMsig[l]^2)}
    
    simulated.N.InSigPow(1,a_ind) = trace(RN_simSig);
    % E{(1/Ne)sum(xNsig[l]^2)}
    
    simulated.CSApp.InSigPow(1,a_ind) = ...
        (simulated.M.InSigPow(1,a_ind) + simulated.N.InSigPow(1,a_ind)) / 2;
    % Avg between subarrays...might be weak definition
    
    simulated.CSAcbf.InSigPow(1,a_ind) = trace(RCSAcbf_simSig);
    % E{(1/MNe)sum(xCSAcbfsig[l]^2)}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Input Noise Power: E{c^2[l]} = varW / (1-alpha^2)
    %%%%    *.InNoisePow is average power of noise
    %%%%    averaged over all trials, averaged over all sensors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    derived.InSigPow = vS'*vS*varS;
    derived.InNoisePow(1,a_ind) = L*varW;
    
    simulated.ULA.InNoisePow(1,a_ind) = trace(ULA_InDataNoise*ULA_InDataNoise'/SampleSize);
    % E{(1/L)sum(xULAnoise[l]^2)}
    
    simulated.M.InNoisePow(1,a_ind) = trace(M_InDataNoise*M_InDataNoise'/SampleSize);
    % E{(1/Me)sum(xMnoise[l]^2)}
    
    simulated.N.InNoisePow(1,a_ind) =  trace(N_InDataNoise*N_InDataNoise'/SampleSize);
    % E{(1/Ne)sum(xNnoise[l]^2)}
    
    simulated.CSApp.InNoisePow(1,a_ind) = ...
        (simulated.M.InNoisePow(1,a_ind) + simulated.N.InNoisePow(1,a_ind)) / 2;
    % Avg between subarrays...might be weak definition
    
    simulated.CSAcbf.InNoisePow(1,a_ind) =  trace(CSAcbf_InDataNoise*CSAcbf_InDataNoise'/SampleSize);
    % E{(1/MNe)sum(xCSAcbfnoise[l]^2)}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      SNRin = varS*(1-alpha^2) / varW
    %%%%    Input SNR simulation calculation stored in *.InSNR
    %%%%    Input Power should be same at all sensors because assuming WSS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    derived.InSNR(1,a_ind) = derived.InSigPow/derived.InNoisePow(1,a_ind);
    
    simulated.ULA.InSNR(1,a_ind) = ...
        simulated.ULA.InSigPow(1,a_ind)/simulated.ULA.InNoisePow(1,a_ind);
    
    simulated.M.InSNR(1,a_ind) = ...
        simulated.M.InSigPow(1,a_ind)/simulated.M.InNoisePow(1,a_ind);
    
    simulated.N.InSNR(1,a_ind) = ...
        simulated.N.InSigPow(1,a_ind)/simulated.N.InNoisePow(1,a_ind);
    
    simulated.CSApp.InSNR(1,a_ind) = ...
        simulated.CSApp.InSigPow(1,a_ind)/simulated.CSApp.InNoisePow(1,a_ind);
    
    simulated.CSAcbf.InSNR(1,a_ind) = ...
        simulated.CSAcbf.InSigPow(1,a_ind)/simulated.CSAcbf.InNoisePow(1,a_ind);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Output Signal Power: varS
    %%%%    *.OutDataSig contains Conventionally BeamFormed (CBF) signal data
    %%%%    *.OutSigPow is average power of signal after CBF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ULA_OutDataSig = wULA'*ULA_InDataSig;
    simulated.ULA.OutSigPow(1,a_ind) = ULA_OutDataSig*ULA_OutDataSig'/SampleSize;
    % E{yULAsig^2}
    
    M_OutDataSig = wM'*M_InDataSig;
    simulated.M.OutSigPow(1,a_ind) = M_OutDataSig*M_OutDataSig'/SampleSize;
    % E{yMsig^2}
    
    N_OutDataSig = wN'*N_InDataSig;
    simulated.N.OutSigPow(1,a_ind) = N_OutDataSig*N_OutDataSig'/SampleSize;
    % E{yNsig^2}
    
    CSAcbf_OutDataSig = wCSA'*CSAcbf_InDataSig;
    simulated.CSAcbf.OutSigPow(1,a_ind) = CSAcbf_OutDataSig*CSAcbf_OutDataSig'/SampleSize;
    % E{yCSAcbfsig^2}
    
    simulated.CSApp.OutSigPow(1,a_ind) = M_OutDataSig*N_OutDataSig'/SampleSize;
    % yCSAprodsig = yMsig.*conj(yNsig)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Output Signal + Noise Power
    %%%%    *.OutData contains Conventionally BeamFormed (CBF) signal + noise data
    %%%%    *.OutDataPow is average power of signal + noise after CBF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ULA_OutDataSigNoise = wULA'*ULA_InDataSigNoise;
    ULA_OutDataSigNoisePow(1,a_ind) = ULA_OutDataSigNoise*ULA_OutDataSigNoise'/SampleSize;
    
    M_OutDataSigNoise = wM'*M_InDataSigNoise;
    simulated.M.OutDataSigNoisePow(1,a_ind) = M_OutDataSigNoise*M_OutDataSigNoise'/SampleSize;
    
    N_OutDataSigNoise = wN'*N_InDataSigNoise;
    simulated.N.OutDataSigNoisePow(1,a_ind) = N_OutDataSigNoise*N_OutDataSigNoise'/SampleSize;
    
    CSAcbf_OutDataSigNoise = wCSA'*CSAcbf_InDataSigNoise;
    simulated.CSAcbf.OutDataSigNoisePow(1,a_ind) = CSAcbf_OutDataSigNoise*CSAcbf_OutDataSigNoise'/SampleSize;
    
    simulated.CSApp.OutDataSigNoisePow(1,a_ind) = M_OutDataSigNoise*N_OutDataSigNoise'/SampleSize;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Output Noise Power simulation calculation
    %%%%    *.OutDataNoise contains Conventionally BeamFormed (CBF) noise data
    %%%%    *.OutNoisePow is average power of noise after CBF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%------------------------------------------------------------------------------------------------
    %%      ULA Output Noise Power: E{[(1/L)sum/l(c[l])]^2}
    %%%%    (1/L^2)E{sum/a(c[a])sum/b(conj(c[b]))} =
    %%%%    (1/L^2) * (varW/(1-alpha^2)) * sum(derived.ULA.coarray[k]*alpha^|k|)
    %%%%    derived.ULA.coarray[k] is auto-correlation of ULA
    %%%%    treating ULA as special case of product process, M=1 & N=1
    %%%%------------------------------------------------------------------------------------------------
    derived.ULA.OutNoisePow(1,a_ind) = ...
        (wULA'*wULA)/L*varW*sum(derived.ULA.coarray.*alpha(a_ind).^abs(gamma).*exp(1i*pi.*gamma*uS));
    
    ULA_OutDataNoise = wULA'*ULA_InDataNoise;
    simulated.ULA.OutNoisePow(1,a_ind) = ULA_OutDataNoise*ULA_OutDataNoise'/SampleSize;
    % E{yULAnoise^2}
    simulated.ULA.PowOutNoisePow(1,a_ind) = ...
        abs(mean(ULA_OutDataNoise.*conj(ULA_OutDataNoise)...
        .*conj(ULA_OutDataNoise).*ULA_OutDataNoise));
    
    %%%%------------------------------------------------------------------------------------------------
    %%      M Output Noise Power: E{[(1/Me)sum/l(c[lM])]^2}
    %%%%    (1/Me^2)E{sum/a(c[aM])sum/b(conj(c[bM]))} =
    %%%%    (1/Me^2) * (varW/(1-alpha^2)) * sum(derived.M.coarray[k]*alpha^|k|)
    %%%%    derived.M.coarray[k] is auto-correlation of M-array
    %%%%------------------------------------------------------------------------------------------------
    derived.M.OutNoisePow(1,a_ind) = ...
        (wM'*wM)/I_M*varW*sum(derived.M.coarray.*alpha(a_ind).^abs(gamma).*exp(1i*pi.*gamma*uS));
    
    M_OutDataNoise = wM'*M_InDataNoise;
    simulated.M.OutNoisePow(1,a_ind) = M_OutDataNoise*M_OutDataNoise'/SampleSize;
    % E{yMnoise^2}
    simulated.M.PowOutNoisePow(1,a_ind) = ...
        abs(mean(M_OutDataNoise.*conj(M_OutDataNoise)...
        .*conj(M_OutDataNoise).*M_OutDataNoise));
    
    %%%%------------------------------------------------------------------------------------------------
    %%      N Output Noise Power: E{[(1/Ne)sum/l(c[lN])]^2}
    %%%%    (1/Ne^2)E{sum/a(c[aN])sum/b(conj(c[bN]))} =
    %%%%    (1/Ne^2) * (varW/(1-alpha^2)) * sum(derived.N.coarray[k]*alpha^|k|)
    %%%%    derived.N.coarray[k] is auto-correlation of N-array
    %%%%------------------------------------------------------------------------------------------------
    derived.N.OutNoisePow(1,a_ind) = ...
        (wN'*wN)/I_N*varW*sum(derived.N.coarray.*alpha(a_ind).^abs(gamma).*exp(1i*pi.*gamma*uS));
    
    N_OutDataNoise = wN'*N_InDataNoise;
    simulated.N.OutNoisePow(1,a_ind) = N_OutDataNoise*N_OutDataNoise'/SampleSize;
    % E{yNnoise^2}
    simulated.N.PowOutNoisePow(1,a_ind) = ...
        abs(mean(N_OutDataNoise.*conj(N_OutDataNoise)...
        .*conj(N_OutDataNoise).*N_OutDataNoise));
    
    %%%%------------------------------------------------------------------------------------------------
    %%      CSAcbf Output Noise Power:
    %%%%    (1/MNe^2)E{sum/a(c[a(M or N)])sum/b(conj(c[b(M or N)]))} =
    %%%%    (1/MNe^2) * (varW/(1-alpha^2)) * sum(derived.CSAcbf.coarray[k]*alpha^|k|)
    %%%%    derived.CSAcbf.coarray[k] is auto-correlation of total CSA
    %%%%------------------------------------------------------------------------------------------------
    derived.CSAcbf.OutNoisePow(1,a_ind) = ...
        (wCSA'*wCSA)/I_csa*varW*sum(derived.CSAcbf.coarray.*alpha(a_ind).^abs(gamma).*exp(1i*pi.*gamma*uS));
    CSAcbf_OutDataNoise = wCSA'*CSAcbf_InDataNoise;
    simulated.CSAcbf.OutNoisePow(1,a_ind) = CSAcbf_OutDataNoise*CSAcbf_OutDataNoise'/SampleSize;
    % E{yCSAcbfnoise^2}
    
    %%%%------------------------------------------------------------------------------------------------
    %%      CSApp Output Noise Power:
    %%%%    (MN/L^2)E{sum/a(c[aM])sum/b(conj(c[bN]))} =
    %%%%    (MN/L^2) * (varW/(1-alpha^2)) * sum(derived.CSApp.coarray[k]*alpha^|k|)
    %%%%    derived.CSApp.coarray[k] is cross-correlation of M- and N- subarrays
    %%%%------------------------------------------------------------------------------------------------
    derived.CSApp.OutNoisePow(1,a_ind) = ...
        (wM'*wN)/beta*varW*sum(derived.CSApp.coarray.*alpha(a_ind).^abs(gamma).*exp(1i*pi.*gamma*uS));
    simulated.CSApp.OutNoisePow(1,a_ind) = M_OutDataNoise*N_OutDataNoise'/SampleSize;
    % E{yMnoise.*conj(yNnoise)}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%      Output SNR calculation using deflection statistics
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%------------------------------------------------------------------------------------------------
    %%      ULA: (varS*L^2*(1-alpha^2)) / (varW*sum(derived.ULA.r[k]*alpha^|k|))
    %%%%------------------------------------------------------------------------------------------------
    derived.ULA.OutSNR(1,a_ind) = varS/derived.ULA.OutNoisePow(1,a_ind);
    simulated.ULA.OutSNR(1,a_ind) = ...
        (ULA_OutDataSigNoisePow(1,a_ind) - simulated.ULA.OutNoisePow(1,a_ind)) / ...
        sqrt(simulated.ULA.PowOutNoisePow(1,a_ind) - simulated.ULA.OutNoisePow(1,a_ind)^2);
    
    %%%%------------------------------------------------------------------------------------------------
    %%      M: (varS*Me^2*(1-alpha^2)) / (varW*sum(derived.M.r[k]*alpha^|k|))
    %%%%------------------------------------------------------------------------------------------------
    derived.M.OutSNR(1,a_ind) = varS/derived.M.OutNoisePow(1,a_ind);
    simulated.M.OutSNR(1,a_ind) = ...
        (simulated.M.OutDataSigNoisePow(1,a_ind) - simulated.M.OutNoisePow(1,a_ind)) / ...
        sqrt(simulated.M.PowOutNoisePow(1,a_ind) - simulated.M.OutNoisePow(1,a_ind)^2);
    
    %%%%------------------------------------------------------------------------------------------------
    %%      N: (varS*Ne^2*(1-alpha^2)) / (varW*sum(derived.N.r[k]*alpha^|k|))
    %%%%------------------------------------------------------------------------------------------------
    derived.N.OutSNR(1,a_ind) = varS/derived.N.OutNoisePow(1,a_ind);
    simulated.N.OutSNR(1,a_ind) = ...
        (simulated.N.OutDataSigNoisePow(1,a_ind) - simulated.N.OutNoisePow(1,a_ind)) / ...
        sqrt(simulated.N.PowOutNoisePow(1,a_ind) - simulated.N.OutNoisePow(1,a_ind)^2);
    
    %%%%------------------------------------------------------------------------------------------------
    %%      CSAcbf: (varS*MNe^2*(1-alpha^2)) / (varW*sum(derived.CSAcbf.r[k]*alpha^|k|))
    %%%%------------------------------------------------------------------------------------------------
    derived.CSAcbf.OutSNR(1,a_ind) = varS/derived.CSAcbf.OutNoisePow(1,a_ind);
    simulated.CSAcbf.OutSNR(1,a_ind) = ...
        simulated.CSAcbf.OutSigPow(1,a_ind)/simulated.CSAcbf.OutNoisePow(1,a_ind);
    
    %%%%------------------------------------------------------------------------------------------------
    %%      2nd Moment of DATA y.*conj(y) = (yM.*conj(yN).*conj(yM).*yN)
    %%%%    E{y^2|H0} = M^2*N^2/L^4 * varW^2/(1-alpha^2)^2 * ...
    %%%%    SUM/(a,b,d,e)(alpha^(|aM-bN|+|dM-eN|) + alpha^(M|a-d|+N|b-e|)
    %%%%------------------------------------------------------------------------------------------------
    derived.CSApp.PowOutNoisePow(1,a_ind) = (1/I_M^2/I_N^2)*varW^2*...
        sum(alpha(a_ind).^(derived.CSApp.PTS).*SDT);
    simulated.CSApp.PowOutNoisePow(1,a_ind) = ...
        abs(mean(M_OutDataNoise.*conj(N_OutDataNoise)...
        .*conj(M_OutDataNoise).*N_OutDataNoise));
    
    %%%%------------------------------------------------------------------------------------------------
    %%%%    CSApp.OutSNR will be defined as the deflection equation given by H. Cox
    %%%%    SNRout = E{y|H1}-E{y|H0} / sqrt(E{y^2|H0}-E^2{y|H0})
    %%%%    SNRMNoutprod = L^2 / (M*N*varW*(1-alpha^2) / ...
    %%%%    sqrt( SUM/(a,b,d,e)(alpha^{|aM-bN|+|dM-eN|}+alpha^[M|a-d|+N|b-e|])-(rMN[k]*alpha^|k|)^2 )
    %%%%------------------------------------------------------------------------------------------------
    derived.CSApp.OutSNR(1,a_ind) = varS / ...
        sqrt(derived.CSApp.PowOutNoisePow(1,a_ind) - derived.CSApp.OutNoisePow(1,a_ind)^2);
    simulated.CSApp.OutSNR(1,a_ind) = ...
        (simulated.CSApp.OutDataSigNoisePow(1,a_ind) - simulated.CSApp.OutNoisePow(1,a_ind))/...
        sqrt(simulated.CSApp.PowOutNoisePow(1,a_ind) - simulated.CSApp.OutNoisePow(1,a_ind)^2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%    Noise/Signal+Noise PDF & CDF calculations:
    %%%%    Dealing with n=2 dimenional Guassian RVs for Real and Imaginary components
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for localdummy = 1
            %%%%------------------------------------------------------------------------------------------------
            %%      ULA.PDF: Power of ULA output is distributed as a central Chi-Squared random variable
            %%%%    pULA(y) = 1/(2*sigmaULA^2) * exp(-y/(2*sigmaULA^2))
            %%%%    defined in Simon p13
            %%%%------------------------------------------------------------------------------------------------
            ULA_OutDataSigNoise = wULA'*ULA_InDataSigNoise;
            ULASigNoiseData = ULA_OutDataSigNoise.*conj(ULA_OutDataSigNoise);
            
            % histogram for simulated ULA Signal + Noise PDF
            [yULASigNoisePDF,xULASigNoisePDF] = ...
                hist(ULASigNoiseData,100);
            simulated.ULA.xSigNoisePDF{1,a_ind} = xULASigNoisePDF;
            simulated.ULA.SigNoisePDF{1,a_ind} = yULASigNoisePDF./trapz(xULASigNoisePDF,yULASigNoisePDF);
            % x-axis vector for plotting derived ULA Signal + Noise PDF
            derived.ULA.xSigNoisePDF{1,a_ind} = linspace(xULASigNoisePDF(1)/100,xULASigNoisePDF(end),500);
            % derived ULA Signal + Noise variance
            varULASigNoise = varS/2 + derived.ULA.OutNoisePow(1,a_ind)/2;
            % derived ULA Signal + Noise PDF
            H_ULASigNoisePDF = @(x)Chi2PDF(x,varULASigNoise);
            derived.ULA.SigNoisePDF{1,a_ind} = H_ULASigNoisePDF(derived.ULA.xSigNoisePDF{1,a_ind});
            
            % Calculating Pd for all x values
            if doROC
                for tauind = 1:length(derived.ULA.xSigNoisePDF{1,a_ind})
                    derived.ULA.Pd{1,a_ind}(tauind) = ...
                        quadgk(H_ULASigNoisePDF,derived.ULA.xSigNoisePDF{1,a_ind}(tauind),inf);
                end
                derived.ULA.Pd{1,a_ind} = sort(derived.ULA.Pd{1,a_ind},'ascend');
            end
            
            if ~islocal
                % calculating simulated ULA Signal + Noise CDF
                [simulated.ULA.SigNoiseCDF{1,a_ind},simulated.ULA.xSigNoiseCDF{1,a_ind}] = ...
                    ecdf(ULASigNoiseData);
                % x-axis vector for plotting derived ULA CDF
                derived.ULA.xSigNoiseCDF{1,a_ind} = simulated.ULA.xSigNoiseCDF{1,a_ind};
                % derived ULA Signal + Noise CDF
                derived.ULA.SigNoiseCDF{1,a_ind} = ...
                    1 - exp(-derived.ULA.xSigNoiseCDF{1,a_ind}/(2*varULASigNoise));
                % ULA Signal + Noise KS Test
                derived.ULA.SigNoiseKSTest(1,a_ind) = ...
                    kstest(ULASigNoiseData,[derived.ULA.xSigNoiseCDF{1,a_ind},derived.ULA.SigNoiseCDF{1,a_ind}]);
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      ULA.NoisePDF: Power of ULA output is distributed as a central Chi-Squared random variable
            %%%%    pULA(y) = 1/(2*sigmaULA^2) * exp(-y/(2*sigmaULA^2))
            %%%%    defined in Simon p13
            %%%%------------------------------------------------------------------------------------------------
            % simulated Noise Only sample data
            ULA_OutDataNoise = wULA'*ULA_InDataNoise;
            ULANoiseData = ULA_OutDataNoise.*conj(ULA_OutDataNoise);
            
            % histogram for simulated ULA Noise Only PDF
            [yULANoisePDF,xULANoisePDF] = hist(ULANoiseData,100);
            % storing x- and y- values
            simulated.ULA.xNoisePDF{1,a_ind} = xULANoisePDF;
            simulated.ULA.NoisePDF{1,a_ind} = yULANoisePDF./trapz(xULANoisePDF,yULANoisePDF);
            % x-axis vector for plotting derived ULA Noise Only PDF
            derived.ULA.xNoisePDF{1,a_ind} = linspace(xULANoisePDF(1)/100,xULANoisePDF(end),500);
            % derived ULA Noise Only variance
            varULANoise = derived.ULA.OutNoisePow(1,a_ind)/2;
            % derived ULA PDF
            H_ULANoisePDF = @(x)Chi2PDF(x,varULANoise);
            derived.ULA.NoisePDF{1,a_ind} = H_ULANoisePDF(derived.ULA.xNoisePDF{1,a_ind});
            
            
            if ~islocal
                % calculating simulated ULA Noise Only CDF
                [simulated.ULA.NoiseCDF{1,a_ind},simulated.ULA.xNoiseCDF{1,a_ind}] = ...
                    ecdf(ULANoiseData);
                % x-axis vector for plotting derived ULA Noise Only CDF
                derived.ULA.xNoiseCDF{1,a_ind} = simulated.ULA.xNoiseCDF{1,a_ind};
                % derived ULA Noise Only CDF
                derived.ULA.NoiseCDF{1,a_ind} = 1 - exp(-derived.ULA.xNoiseCDF{1,a_ind}/(2*varULANoise));
                % ULA Noise Only KS Test
                derived.ULA.NoiseKSTest(1,a_ind) = ...
                    kstest(ULANoiseData,[derived.ULA.xNoiseCDF{1,a_ind},derived.ULA.NoiseCDF{1,a_ind}]);
                % Calculating Pfa for all x values
            end
            if doROC
                for tauind = 1:length(derived.ULA.xSigNoisePDF{1,a_ind})
                    derived.ULA.Pfa{1,a_ind}(tauind) = ...
                        quadgk(H_ULANoisePDF,derived.ULA.xSigNoisePDF{1,a_ind}(tauind),inf);
                end
                derived.ULA.Pfa{1,a_ind} = sort(derived.ULA.Pfa{1,a_ind},'ascend');
                for tauind = 1:length(xULANoisePDF)
                    Detection = ULASigNoiseData > xULANoisePDF(tauind);
                    simulated.ULA.Pd{1,a_ind}(tauind) = sum(Detection)/SampleSize;
                    FalseAlarm = ULANoiseData > xULANoisePDF(tauind);
                    simulated.ULA.Pfa{1,a_ind}(tauind) = sum(FalseAlarm)/SampleSize;
                end
            end
            
            
            %%%%------------------------------------------------------------------------------------------------
            %%      M.NoisePDF: Magnitude of M output is distributed as a Rayleigh random variable
            %%%%    pN(r) = r/sigmaN^2 * exp(-r^2/(2*sigmaM^2))
            %%%%    defined in Simon p10
            %%%%------------------------------------------------------------------------------------------------
            MNoiseData = abs(M_OutDataNoise);
            
            % histogram for simulated M Noise Only PDF
            [simulated.M.NoisePDF{1,a_ind},simulated.M.xNoisePDF{1,a_ind}] = ...
                hist(MNoiseData,100);
            simulated.M.NoisePDF{1,a_ind} = simulated.M.NoisePDF{1,a_ind}...
                ./trapz(simulated.M.xNoisePDF{1,a_ind},simulated.M.NoisePDF{1,a_ind});
            % x-axis for plotting derived M Noise Only PDF
            derived.M.xNoisePDF{1,a_ind} = linspace(min(simulated.M.xNoisePDF{1,a_ind})/100,...
                max(simulated.M.xNoisePDF{1,a_ind}),500);
            % derived M Noise Only standard deviation
            sigmaM_Noise = sqrt(derived.M.OutNoisePow(1,a_ind)/2);
            % derived M Noise Only PDF
            derived.M.NoisePDF{1,a_ind} = derived.M.xNoisePDF{1,a_ind}./sigmaM_Noise^2 .* ...
                exp(-(derived.M.xNoisePDF{1,a_ind}.^2./(2*sigmaM_Noise^2)));
            
            if ~islocal
            % calculating simulated M Noise Only CDF
            [simulated.M.NoiseCDF{1,a_ind},simulated.M.xNoiseCDF{1,a_ind}] = ...
                ecdf(MNoiseData);
            % x-axis vector for plotting derived M Noise Only CDF
            derived.M.xNoiseCDF{1,a_ind} = simulated.M.xNoiseCDF{1,a_ind};
            % derived M Noise Only CDF
            derived.M.NoiseCDF{1,a_ind} = ...
                1 - exp(-derived.M.xNoiseCDF{1,a_ind}.^2/(2*sigmaM_Noise^2));
            % M KS Noise Only Test
            derived.M.NoiseKSTest(1,a_ind) = ...
                kstest(MNoiseData,[derived.M.xNoiseCDF{1,a_ind},derived.M.NoiseCDF{1,a_ind}]);
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      M.SigNoisePDF: Magnitude of M output is distributed as a Rayleigh random variable
            %%%%    pN(r) = r/sigmaN^2 * exp(-r^2/(2*sigmaM^2))
            %%%%    defined in Simon p10
            %%%%------------------------------------------------------------------------------------------------
            simulated.M.OutDataSigNoise = wM'*M_InDataSigNoise;
            MSigNoiseData = abs(simulated.M.OutDataSigNoise);
            
            % histogram for simulated M Signal + Noise PDF
            [simulated.M.SigNoisePDF{1,a_ind},simulated.M.xSigNoisePDF{1,a_ind}] = ...
                hist(MSigNoiseData,100);
            simulated.M.SigNoisePDF{1,a_ind} = simulated.M.SigNoisePDF{1,a_ind}...
                ./trapz(simulated.M.xSigNoisePDF{1,a_ind},simulated.M.SigNoisePDF{1,a_ind});
            % x-axis for plotting derived M Signal + Noise PDF
            derived.M.xSigNoisePDF{1,a_ind} = linspace(min(simulated.M.xSigNoisePDF{1,a_ind})/100,...
                max(simulated.M.xSigNoisePDF{1,a_ind}),500);
            % derived M Signal + Noise standard deviation
            sigmaM_SigNoise = sqrt(varS/2 + derived.M.OutNoisePow(1,a_ind)/2);
            % derived M Signal + Noise PDF
            derived.M.SigNoisePDF{1,a_ind} = derived.M.xSigNoisePDF{1,a_ind}./sigmaM_SigNoise^2 ...
                .* exp(-(derived.M.xSigNoisePDF{1,a_ind}.^2./(2*sigmaM_SigNoise^2)));
            
            if ~islocal
            % calculating simulated M Signal + Noise CDF
            [simulated.M.SigNoiseCDF{1,a_ind},simulated.M.xSigNoiseCDF{1,a_ind}] = ...
                ecdf(MSigNoiseData);
            % x-axis vector for plotting derived M Signal + Noise CDF
            derived.M.xSigNoiseCDF{1,a_ind} = simulated.M.xSigNoiseCDF{1,a_ind};
            % derived M Signal + Noise CDF
            derived.M.SigNoiseCDF{1,a_ind} = ...
                1 - exp(-derived.M.xSigNoiseCDF{1,a_ind}.^2/(2*sigmaM_SigNoise^2));
            % M Signal + Noise KS Test
            derived.M.SigNoiseKSTest(1,a_ind) = ...
                kstest(MSigNoiseData,[derived.M.xSigNoiseCDF{1,a_ind},derived.M.SigNoiseCDF{1,a_ind}]);
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      N.PDF: Magnitude of N output is distributed as a Rayleigh random variable
            %%%%    pN(r) = r/sigmaN^2 * exp(-r^2/(2*sigmaN^2))
            %%%%    defined in Simon p10
            %%%%------------------------------------------------------------------------------------------------
            NNoiseData = abs(N_OutDataNoise);
            
            % histogram for simulated N Noise Only PDF
            [simulated.N.NoisePDF{1,a_ind},simulated.N.xNoisePDF{1,a_ind}] = ...
                hist(NNoiseData,100);
            simulated.N.NoisePDF{1,a_ind} = simulated.N.NoisePDF{1,a_ind}...
                ./trapz(simulated.N.xNoisePDF{1,a_ind},simulated.N.NoisePDF{1,a_ind});
            % x-axis for plotting derived N Noise Only PDF
            derived.N.xNoisePDF{1,a_ind} = linspace(min(simulated.N.xNoisePDF{1,a_ind})/100,...
                max(simulated.N.xNoisePDF{1,a_ind}),500);
            % derived N Noise Only standard deviation
            sigmaN_Noise = sqrt(derived.N.OutNoisePow(1,a_ind)/2);
            % derived N Noise Only PDF
            derived.N.NoisePDF{1,a_ind} = derived.N.xNoisePDF{1,a_ind}./sigmaN_Noise^2 .* ...
                exp(-(derived.N.xNoisePDF{1,a_ind}.^2./(2*sigmaN_Noise^2)));
            
            if ~islocal
            % calculating simulated N Noise Only CDF
            [simulated.N.NoiseCDF{1,a_ind},simulated.N.xNoiseCDF{1,a_ind}] = ...
                ecdf(NNoiseData);
            % x-axis vector for plotting derived N Noise Only CDF
            derived.N.xNoiseCDF{1,a_ind} = simulated.N.xNoiseCDF{1,a_ind};
            % derived N Noise Only CDF
            derived.N.NoiseCDF{1,a_ind} = ...
                1 - exp(-derived.N.xNoiseCDF{1,a_ind}.^2/(2*sigmaN_Noise^2));
            % N KS Noise Only Test
            derived.N.NoiseKSTest(1,a_ind) = ...
                kstest(NNoiseData,[derived.N.xNoiseCDF{1,a_ind},derived.N.NoiseCDF{1,a_ind}]);
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      N.PDF: Magnitude of N output is distributed as a Rayleigh random variable
            %%%%    pN(r) = r/sigmaN^2 * exp(-r^2/(2*sigmaN^2))
            %%%%    defined in Simon p10
            %%%%------------------------------------------------------------------------------------------------
            simulated.N.OutDataSigNoise = wN'*N_InDataSigNoise;
            NSigNoiseData = abs(simulated.N.OutDataSigNoise);
            
            % histogram for simulated N PDF
            [simulated.N.SigNoisePDF{1,a_ind},simulated.N.xSigNoisePDF{1,a_ind}] = ...
                hist(NSigNoiseData,100);
            simulated.N.SigNoisePDF{1,a_ind} = simulated.N.SigNoisePDF{1,a_ind}...
                ./trapz(simulated.N.xSigNoisePDF{1,a_ind},simulated.N.SigNoisePDF{1,a_ind});
            % x-axis for plotting derived N PDF
            derived.N.xSigNoisePDF{1,a_ind} = linspace(min(simulated.N.xSigNoisePDF{1,a_ind})/100,...
                max(simulated.N.xSigNoisePDF{1,a_ind}),500);
            % derived N standard deviation
            sigmaN_SigNoise = sqrt(varS/2 + derived.N.OutNoisePow(1,a_ind)/2);
            % derived N PDF
            derived.N.SigNoisePDF{1,a_ind} = derived.N.xSigNoisePDF{1,a_ind}./sigmaN_SigNoise^2 ...
                .* exp(-(derived.N.xSigNoisePDF{1,a_ind}.^2./(2*sigmaN_SigNoise^2)));
            
            if ~islocal
            % calculating simulated N CDF
            [simulated.N.SigNoiseCDF{1,a_ind},simulated.N.xSigNoiseCDF{1,a_ind}] = ...
                ecdf(NSigNoiseData);
            % x-axis vector for plotting derived N CDF
            derived.N.xSigNoiseCDF{1,a_ind} = simulated.N.xSigNoiseCDF{1,a_ind};
            % derived N CDF
            derived.N.SigNoiseCDF{1,a_ind} = ...
                1 - exp(-derived.N.xSigNoiseCDF{1,a_ind}.^2/(2*sigmaN_SigNoise^2));
            % N Signal + Noise KS Test
            derived.N.SigNoiseKSTest(1,a_ind) = ...
                kstest(NSigNoiseData,[derived.N.xSigNoiseCDF{1,a_ind},derived.N.SigNoiseCDF{1,a_ind}]);
            end
            
            %%%%------------------------------------------------------------------------------------------------
            %%      CSAcbf.PDF: Power of CSAcbf output is distributed as a central Chi-Squared random variable
            %%%%    pCSAcbf(y) = 1/(2*sigmaCSAcbf^2) * exp(-y/(2*sigmaCSAcbf^2))
            %%%%    defined in Simon p13
            %%%%------------------------------------------------------------------------------------------------
            simulated.CSAcbf.OutDataSigNoise = wCSA'*CSAcbf_InDataSigNoise;
            CSAcbfSigNoiseData = simulated.CSAcbf.OutDataSigNoise.*conj(simulated.CSAcbf.OutDataSigNoise);
            
            % histogram for simulated CSAcbf PDF
            [simulated.CSAcbf.SigNoisePDF{1,a_ind},simulated.CSAcbf.xSigNoisePDF{1,a_ind}] = ...
                hist(CSAcbfSigNoiseData,100);
            simulated.CSAcbf.SigNoisePDF{1,a_ind} = simulated.CSAcbf.SigNoisePDF{1,a_ind}...
                ./trapz(simulated.CSAcbf.xSigNoisePDF{1,a_ind},simulated.CSAcbf.SigNoisePDF{1,a_ind});
            % x-axis for plotting derived CSAcbf PDF
            derived.CSAcbf.xSigNoisePDF{1,a_ind} = linspace(min(simulated.CSAcbf.xSigNoisePDF{1,a_ind})/100,...
                max(simulated.CSAcbf.xSigNoisePDF{1,a_ind}),1000);
            % derived CSAcbf variance
            varCSAcbfSigNoise = varS/2 + derived.CSAcbf.OutNoisePow(1,a_ind)/2;
            % derived CSAcbf PDF
            H_CSAcbfSigNoisePDF = @(x)Chi2PDF(x,varCSAcbfSigNoise);
            derived.CSAcbf.SigNoisePDF{1,a_ind} = H_CSAcbfSigNoisePDF(derived.CSAcbf.xSigNoisePDF{1,a_ind});
            %             1/(2*varCSAcbfSigNoise).*exp(-(derived.CSAcbf.xSigNoisePDF{1,a_ind}/(2*varCSAcbfSigNoise)));
            
            if ~islocal
            % calculating simulated CSAcbf CDF
            [simulated.CSAcbf.SigNoiseCDF{1,a_ind},simulated.CSAcbf.xSigNoiseCDF{1,a_ind}] = ...
                ecdf(CSAcbfSigNoiseData);
            % x-axis vector for plotting derived CSAcbf CDF
            derived.CSAcbf.xSigNoiseCDF{1,a_ind} = simulated.CSAcbf.xSigNoiseCDF{1,a_ind};
            % derived CSAcbf CDF
            derived.CSAcbf.SigNoiseCDF{1,a_ind} = ...
                1 - exp(-derived.CSAcbf.xSigNoiseCDF{1,a_ind}/(2*varCSAcbfSigNoise));
            % CSAcbf Signal + Noise KS Test
            derived.CSAcbf.SigNoiseKSTest(1,a_ind) = ...
                kstest(CSAcbfSigNoiseData,[derived.CSAcbf.xSigNoiseCDF{1,a_ind},derived.CSAcbf.SigNoiseCDF{1,a_ind}]);
            
                if doROC
                    for tauind = 1:length(derived.CSAcbf.xSigNoisePDF{1,a_ind})
                        derived.CSAcbf.Pd{1,a_ind}(tauind) = ...
                            quadgk(H_CSAcbfSigNoisePDF,derived.CSAcbf.xSigNoisePDF{1,a_ind}(tauind),inf);
                    end
                    derived.CSAcbf.Pd{1,a_ind} = sort(derived.CSAcbf.Pd{1,a_ind},'ascend');

                end
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      CSAcbf.PDF: Power of CSAcbf output is distributed as a central Chi-Squared random variable
            %%%%    pCSAcbf(y) = 1/(2*sigmaCSAcbf^2) * exp(-y/(2*sigmaCSAcbf^2))
            %%%%    defined in Simon p13
            %%%%------------------------------------------------------------------------------------------------
            simulated.CSAcbf.OutDataNoise = wCSA'*CSAcbf_InDataNoise;
            CSAcbfNoiseData = simulated.CSAcbf.OutDataNoise.*conj(simulated.CSAcbf.OutDataNoise);
            
            % histogram for simulated CSAcbf Noise Only PDF
            [simulated.CSAcbf.NoisePDF{1,a_ind},simulated.CSAcbf.xNoisePDF{1,a_ind}] = ...
                hist(CSAcbfNoiseData,100);
            simulated.CSAcbf.NoisePDF{1,a_ind} = simulated.CSAcbf.NoisePDF{1,a_ind}...
                ./trapz(simulated.CSAcbf.xNoisePDF{1,a_ind},simulated.CSAcbf.NoisePDF{1,a_ind});
            % x-axis for plotting derived CSAcbf Noise Only PDF
            derived.CSAcbf.xNoisePDF{1,a_ind} = linspace(min(simulated.CSAcbf.xNoisePDF{1,a_ind})/100,...
                max(simulated.CSAcbf.xNoisePDF{1,a_ind}),1000);
            % derived CSAcbf Noise Only variance
            varCSAcbfNoise = derived.CSAcbf.OutNoisePow(1,a_ind)/2;
            % derived CSAcbf PDF
            H_CSAcbfNoisePDF = @(x)Chi2PDF(x,varCSAcbfNoise);
            derived.CSAcbf.NoisePDF{1,a_ind} = H_CSAcbfNoisePDF(derived.CSAcbf.xNoisePDF{1,a_ind});
            
            if ~islocal
            % calculating simulated CSAcbf Noise Only CDF
            [simulated.CSAcbf.NoiseCDF{1,a_ind},simulated.CSAcbf.xNoiseCDF{1,a_ind}] = ...
                ecdf(CSAcbfNoiseData);
            % x-axis vector for plotting derived CSAcbf Noise Only CDF
            derived.CSAcbf.xNoiseCDF{1,a_ind} = simulated.CSAcbf.xNoiseCDF{1,a_ind};
            % derived CSAcbf Noise Only CDF
            derived.CSAcbf.NoiseCDF{1,a_ind} = ...
                1 - exp(-derived.CSAcbf.xNoiseCDF{1,a_ind}/(2*varCSAcbfNoise));
            derived.CSAcbf.NoiseKSTest(1,a_ind) = ...
                kstest(CSAcbfNoiseData,[derived.CSAcbf.xNoiseCDF{1,a_ind},derived.CSAcbf.NoiseCDF{1,a_ind}]);
            if doROC
                for tauind = 1:length(derived.CSAcbf.xSigNoisePDF{1,a_ind})
                    derived.CSAcbf.Pfa{1,a_ind}(tauind) = ...
                        quadgk(H_CSAcbfNoisePDF,derived.CSAcbf.xSigNoisePDF{1,a_ind}(tauind),inf);
                end
                derived.CSAcbf.Pfa{1,a_ind} = sort(derived.CSAcbf.Pfa{1,a_ind},'ascend');
                for tauind = 1:length(simulated.CSAcbf.xNoisePDF{1,a_ind})
                    Detection = CSAcbfSigNoiseData > simulated.CSAcbf.xNoisePDF{1,a_ind}(tauind);
                    simulated.CSAcbf.Pd{1,a_ind}(tauind) = sum(Detection)/SampleSize;
                    FalseAlarm = CSAcbfNoiseData > simulated.CSAcbf.xNoisePDF{1,a_ind}(tauind);
                    simulated.CSAcbf.Pfa{1,a_ind}(tauind) = sum(FalseAlarm)/SampleSize;
                end
            end
            end
            
            %%%%------------------------------------------------------------------------------------------------
            %%      CSApp.PDF: Power of the product of the subarray outputs is distributed as a product
            %%%%    of two Rayleigh random variables: R = R1*R2
            %%%%    with x = r/(sigmaM*sigmaN*(1-rho^2))
            %%%%    pCSAprod(x) = x/(sigmaM*sigmaN) * I0(x|rho|) * K0(x)
            %%%%    defined in Simon p57
            %%%%------------------------------------------------------------------------------------------------
            CSAprodSigNoiseData = abs(simulated.M.OutDataSigNoise.*conj(simulated.N.OutDataSigNoise));
            
            % histogram for simulated CSApp PDF
            [simulated.CSApp.SigNoisePDF{1,a_ind},simulated.CSApp.xSigNoisePDF{1,a_ind}] = ...
                hist(CSAprodSigNoiseData,100);
            simulated.CSApp.SigNoisePDF{1,a_ind} = simulated.CSApp.SigNoisePDF{1,a_ind}...
                ./trapz(simulated.CSApp.xSigNoisePDF{1,a_ind},simulated.CSApp.SigNoisePDF{1,a_ind});
            % calculating simulated CSApp CDF
            % derived correlation coefficient
            derived.CSApp.rhoSigNoise(1,a_ind) = (varS + derived.CSApp.OutNoisePow(1,a_ind))/2/(sigmaM_SigNoise*sigmaN_SigNoise);
            rhoSigNoiseMat = corrcoef(simulated.M.OutDataSigNoise,simulated.N.OutDataSigNoise);
            simulated.CSApp.rhoSigNoise(1,a_ind) = rhoSigNoiseMat(1,2);
            if ~islocal
            [simulated.CSApp.SigNoiseCDF{1,a_ind},simulated.CSApp.xSigNoiseCDF{1,a_ind}] = ...
                ecdf(CSAprodSigNoiseData);
            % x-axis for plotting derived CSApp PDF
            derived.CSApp.xSigNoisePDF{1,a_ind} = simulated.CSApp.xSigNoiseCDF{1,a_ind};
            derived.CSApp.xSigNoisePDF_forROC = linspace(simulated.CSApp.xSigNoiseCDF{1,a_ind}(1),simulated.CSApp.xSigNoiseCDF{1,a_ind}(end),length(simulated.CSApp.xSigNoiseCDF{1,a_ind}));
            
            H_CSAprodSigNoisePDF = @(x)RxRPDF(x,sigmaM_SigNoise,sigmaN_SigNoise,derived.CSApp.rhoSigNoise(1,a_ind));
            derived.CSApp.SigNoisePDF{1,a_ind} = H_CSAprodSigNoisePDF(derived.CSApp.xSigNoisePDF{1,a_ind});
            
            % x-axis for plotting derived CSApp CDF
            derived.CSApp.xSigNoiseCDF{1,a_ind} = simulated.CSApp.xSigNoiseCDF{1,a_ind};
            % numerically integrated CSApp CDF
            
            end
            if ~islocal
                derived.CSApp.SigNoiseCDF{1,a_ind} = zeros(length(derived.CSApp.xSigNoiseCDF{1,a_ind}),1);
                for xind = 1:length(derived.CSApp.xSigNoiseCDF{1,a_ind})
                    derived.CSApp.SigNoiseCDF{1,a_ind}(xind) = quadgk(H_CSAprodSigNoisePDF,0,derived.CSApp.xSigNoiseCDF{1,a_ind}(xind));
                end
            derived.CSApp.SigNoiseKSTest(1,a_ind) = ...
                kstest(CSAprodSigNoiseData,[derived.CSApp.xSigNoiseCDF{1,a_ind},derived.CSApp.SigNoiseCDF{1,a_ind}]);
            if doROC
                for tauind = 1:length(derived.CSApp.xSigNoisePDF_forROC)
                    derived.CSApp.Pd{1,a_ind}(tauind) = ...
                        quadgk(H_CSAprodSigNoisePDF,derived.CSApp.xSigNoisePDF_forROC(tauind),inf);
                end
                derived.CSApp.Pd{1,a_ind} = sort(derived.CSApp.Pd{1,a_ind},'ascend');
            end
            end
            %%%%------------------------------------------------------------------------------------------------
            %%      CSApp.PDF: Power of the product of the subarray outputs is distributed as a product
            %%%%    of two Rayleigh random variables: R = R1*R2
            %%%%    with x = r/(sigmaM*sigmaN*(1-rho^2))
            %%%%    pCSAprod(x) = x/(sigmaM*sigmaN) * I0(x|rho|) * K0(x)
            %%%%    defined in Simon p57
            %%%%------------------------------------------------------------------------------------------------
            CSAprodNoiseData = abs(M_OutDataNoise.*conj(N_OutDataNoise));
            
            % histogram for simulated CSApp Noise Only PDF
            [simulated.CSApp.NoisePDF{1,a_ind},simulated.CSApp.xNoisePDF{1,a_ind}] = ...
                hist(CSAprodNoiseData,100);
            simulated.CSApp.NoisePDF{1,a_ind} = simulated.CSApp.NoisePDF{1,a_ind}...
                ./trapz(simulated.CSApp.xNoisePDF{1,a_ind},simulated.CSApp.NoisePDF{1,a_ind});
            % calculating simulated CSApp Noise Only CDF
            
            % derived Noise Only correlation coefficient
            derived.CSApp.rhoNoise(1,a_ind) = derived.CSApp.OutNoisePow(1,a_ind)/2/(sigmaM_Noise*sigmaN_Noise);
            rhoNoiseMat = corrcoef(M_OutDataNoise,N_OutDataNoise);
            simulated.CSApp.rhoNoise(1,a_ind) = rhoNoiseMat(1,2);
            if ~islocal
            [simulated.CSApp.NoiseCDF{1,a_ind},simulated.CSApp.xNoiseCDF{1,a_ind}] = ...
                ecdf(CSAprodNoiseData);
            
            % x-axis for plotting derived CSApp Noise Only PDF & CDF
            derived.CSApp.xNoisePDF{1,a_ind} = simulated.CSApp.xNoiseCDF{1,a_ind};
            H_CSAprodNoisePDF = @(x)RxRPDF(x,sigmaM_Noise,sigmaN_Noise,derived.CSApp.rhoNoise(1,a_ind));
            derived.CSApp.NoisePDF{1,a_ind} = H_CSAprodNoisePDF(derived.CSApp.xNoisePDF{1,a_ind});
            %             CSAproductPDF(derived.CSApp.xNoisePDF{1,a_ind},sigmaM_Noise,sigmaN_Noise,derived.CSApp.rhoNoise(1,a_ind));
            
            % x-axis for plotting derived CSApp Noise Only CDF
            derived.CSApp.xNoiseCDF{1,a_ind} = simulated.CSApp.xNoiseCDF{1,a_ind};
            end
            % numerically integrated CSApp Noise Only CDF
            if ~islocal
                derived.CSApp.NoiseCDF{1,a_ind} = zeros(length(derived.CSApp.xNoiseCDF{1,a_ind}),1);
                for xind = 1:length(derived.CSApp.xNoiseCDF{1,a_ind})
                    derived.CSApp.NoiseCDF{1,a_ind}(xind) = quadgk(H_CSAprodNoisePDF,0,derived.CSApp.xNoiseCDF{1,a_ind}(xind));
                end
%             KSTestCSAppNoiseCDF = cumtrapz(derived.CSApp.xNoiseCDF{1,a_ind},derived.CSApp.NoisePDF{1,a_ind});
            derived.CSApp.NoiseKSTest(1,a_ind) = ...
                kstest(CSAprodNoiseData,[derived.CSApp.xNoiseCDF{1,a_ind},derived.CSApp.NoiseCDF{1,a_ind}]);
            end
            % Calculating Pfa for all x values
            if doROC
                for tauind = 1:length(derived.CSApp.xSigNoisePDF_forROC)
                    derived.CSApp.Pfa{1,a_ind}(tauind) = ...
                        quadgk(H_CSAprodNoisePDF,derived.CSApp.xSigNoisePDF_forROC(tauind),inf);
                end
                derived.CSApp.Pfa{1,a_ind} = sort(derived.CSApp.Pfa{1,a_ind},'ascend');
                for tauind = 1:length(simulated.CSApp.xNoisePDF{1,a_ind})
                    Detection = CSAprodSigNoiseData > simulated.CSApp.xNoisePDF{1,a_ind}(tauind);
                    simulated.CSApp.Pd{1,a_ind}(tauind) = sum(Detection)/SampleSize;
                    FalseAlarm = CSAprodNoiseData > simulated.CSApp.xNoisePDF{1,a_ind}(tauind);
                    simulated.CSApp.Pfa{1,a_ind}(tauind) = sum(FalseAlarm)/SampleSize;
                end
            end
    end %for dummylocal
    
end % alpha loop


%% End of simulation. Calculating AG = SNRout/SNRin
derived.ULA.AG = derived.ULA.OutSNR./derived.InSNR;
simulated.ULA.AG = simulated.ULA.OutSNR./simulated.ULA.InSNR;

derived.M.AG = derived.M.OutSNR./derived.InSNR;
simulated.M.AG = simulated.M.OutSNR./simulated.M.InSNR;

derived.N.AG = derived.N.OutSNR./derived.InSNR;
simulated.N.AG = simulated.N.OutSNR./simulated.N.InSNR;

derived.CSAcbf.AG = derived.CSAcbf.OutSNR./derived.InSNR;
simulated.CSAcbf.AG = simulated.CSAcbf.OutSNR./simulated.CSAcbf.InSNR;

derived.CSApp.AG = derived.CSApp.OutSNR./derived.InSNR;
simulated.CSApp.AG = simulated.CSApp.OutSNR./simulated.CSApp.InSNR;
derived.CSApp.FalseAG = (varS./derived.CSApp.OutNoisePow)./(varS./derived.InNoisePow);
simulated.CSApp.FalseAG = 1./simulated.CSApp.OutNoisePow;
profile report

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Array Gain plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
present(0)
load CSAsim2June.mat;
f1 = figure(1);
set(f1,'DefaultAxesFontSize',24)
for b = 1:1
    clf
    semilogx(1-alpha,10*log10(abs(derived.ULA.AG(b,:))),'color',plotcolor.ULA)
    hold on % needs to be after first semilogx command
    semilogx(1-alpha,10*log10(abs(derived.CSAcbf.AG(b,:))),'color',plotcolor.CSAcbf)
    semilogx(1-alpha,10*log10(abs(derived.CSApp.AG(b,:))),'color',plotcolor.CSApp)
    semilogx(1-alpha,10*log10(abs(simulated.ULA.AG(b,:))),':^','MarkerSize',12,'Color',plotcolor.ULA)
    semilogx(1-alpha,10*log10(abs(simulated.CSAcbf.AG(b,:))),':+','MarkerSize',15,'color',plotcolor.CSAcbf)
    semilogx(1-alpha,10*log10(abs(simulated.CSApp.AG(b,:))),':x','MarkerSize',18,'color',plotcolor.CSApp)
    grid on
    ax = gca;
    xtick = linspace(.1,1,10); set(ax,'XTick',xtick);
    xticklabel = abs(round(xtick.*100)./100 - 1);
    set(ax,'Xdir','reverse','XTickLabel',xticklabel)
    xlim([1-.9 1]);
    ylim([7 inf]);
    xlabel('$\alpha$'); ylabel('[dB]')
    legend('','','','ULA','CSA$_{\mathrm{cbf}}$','CSA$_{\mathrm{pp}}$','location','northeast')
    title({'Array Gain vs. $\alpha$ - Fixed Aperture Constraint ($\approx 50\lambda_o$)',['CSA Parameters: $M = 2$, $N = 5$, $\beta = 10$']})
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Noise Only Biased PSD estimate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiasedPSDanimation(length(alpha)) = struct('cdata',[],'colormap',[]);
plotcolor.truePSD = [0 1 1];
zeta = {'ULA';'CSAcbf';'CSApp'};
whitenoise = 1;
fullcorr = length(alpha);
f2 = figure(2);
set(f2,'DefaultAxesFontSize',24)
for aplot = 1
    clf
    hold on
%     box on
    for plotind = 1:length(zeta)
        plot(u,10*log10(abs(derived.(zeta{plotind}).NoisePSD{1,aplot})),'color',plotcolor.(zeta{plotind}))
        plot(u(1:6:end),10*log10(abs(simulated.(zeta{plotind}).NoisePSD{1,aplot}(1:6:end))),':d','MarkerSize',7,'LineWidth',2,'color',plotcolor.(zeta{plotind}))
    end% create a new pair of axes inside current figure
    plot(u,10*log10(abs(truePSDNoise{1,aplot})),'-.','LineWidth',3,'color',plotcolor.truePSD)
    title(['Noise Only PSD Estimates: $\alpha = $ ' num2str(round(alpha(aplot)*1000)/1000) ', $\sigma^2_W = $ ' num2str(varW)])
    legend({'ULA','','CSA$_{\mathrm{cbf}}$','','CSA$_{\mathrm{pp}}$','','true PSD'},'location','northeast')
    xlabel('u'); ylabel('[dB]')
%     ylim([10*log10(varW)-6 6])
    
    BiasedPSDanimation(aplot) = getframe(gcf);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Noise Only Probability/Cumulative Distribution Function plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NoisePDFanimation(length(alpha)) = struct('cdata',[],'colormap',[]);
f3 = figure(3);
set(f3,'DefaultAxesFontSize',24)
    for aplot = 1
        if ~islocal
%         clf
        skip_data = 4;
        hold on
        plot(derived.ULA.xNoisePDF{1,aplot},derived.ULA.NoisePDF{1,aplot}...
            /trapz(derived.ULA.xNoisePDF{1,aplot},derived.ULA.NoisePDF{1,aplot})...
            ,'LineWidth',4,'color',plotcolor.ULA);
        plot(simulated.ULA.xNoisePDF{1,aplot}(1:end-80),simulated.ULA.NoisePDF{1,aplot}(1:end-80)...
            ,':d','color',plotcolor.ULA);
        plot(simulated.ULA.xNoisePDF{1,aplot}(end-80:skip_data:end),simulated.ULA.NoisePDF{1,aplot}(end-80:skip_data:end)...
            ,':d','color',plotcolor.ULA);
        
        hold on
        plot(derived.CSAcbf.xNoisePDF{1,aplot},derived.CSAcbf.NoisePDF{1,aplot}...
            ,'LineWidth',4,'color',plotcolor.CSAcbf);
        plot(simulated.CSAcbf.xNoisePDF{1,aplot}(1:end-80),simulated.CSAcbf.NoisePDF{1,aplot}(1:end-80)...
            ,':+','MarkerSize',15,'color',plotcolor.CSAcbf);
        plot(simulated.CSAcbf.xNoisePDF{1,aplot}(end-80:skip_data:end),simulated.CSAcbf.NoisePDF{1,aplot}(end-80:skip_data:end)...
            ,':+','MarkerSize',15,'color',plotcolor.CSAcbf);
        

        plot(derived.CSApp.xNoisePDF{1,aplot},abs(derived.CSApp.NoisePDF{1,aplot})...
            ,'LineWidth',4,'color',plotcolor.CSApp);
        plot(simulated.CSApp.xNoisePDF{1,aplot}(1:end-80),simulated.CSApp.NoisePDF{1,aplot}(1:end-80)...
            ,':x','MarkerSize',18,'color',plotcolor.CSApp);
        plot(simulated.CSApp.xNoisePDF{1,aplot}(end-80:skip_data:end),simulated.CSApp.NoisePDF{1,aplot}(end-80:skip_data:end)...
            ,':x','MarkerSize',18,'color',plotcolor.CSApp);
        
        
        hold on
        axis tight
        xlabel('$z$')
        title('$H_0$')
        ylabel('$M=5$, $N=6$, $\beta=10$')
        
        end
    end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Cumulative Distribution Function plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoiseCDFanimation(length(alpha)) = struct('cdata',[],'colormap',[]);

for aplot = 1
    if ~islocal
    f4 = figure(4);
    set(f4,'DefaultAxesFontSize',24)
    clf
    
    subplot(2,3,[2 5])
    hold on
    plot(derived.CSAcbf.xNoiseCDF{1,aplot},derived.CSAcbf.NoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xNoiseCDF{1,aplot}(1:100:end-20),...
        simulated.CSAcbf.NoiseCDF{1,aplot}(1:100:end-20),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xNoiseCDF{1,aplot}(end-50:5:end),...
        simulated.CSAcbf.NoiseCDF{1,aplot}(end-50:5:end),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    title(['Noise Only CDFs, $\alpha$ = ' num2str(round(alpha(aplot)*100)/100)])
    xlabel('$z$')
    
    subplot(2,3,[1 4])
    hold on
    plot(derived.ULA.xNoiseCDF{1,aplot},derived.ULA.NoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.ULA);
    plot(simulated.ULA.xNoiseCDF{1,aplot}(1:100:end-20),...
        simulated.ULA.NoiseCDF{1,aplot}(1:100:end-20),':d','MarkerSize',9,'color',plotcolor.ULA);
    plot(simulated.ULA.xNoiseCDF{1,aplot}(end-50:5:end),...
        simulated.ULA.NoiseCDF{1,aplot}(end-50:5:end),':d','MarkerSize',9,'color',plotcolor.ULA);
    
    subplot(2,3,[3 6])
    hold on
    plot(derived.CSApp.xNoiseCDF{1,aplot},derived.CSApp.NoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xNoiseCDF{1,aplot}(1:150:end-20),...
        simulated.CSApp.NoiseCDF{1,aplot}(1:150:end-20),':d','MarkerSize',9,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xNoiseCDF{1,aplot}(end-50:5:end),...
        simulated.CSApp.NoiseCDF{1,aplot}(end-50:5:end),':d','MarkerSize',9,'color',plotcolor.CSApp);
    end
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Signal + Noise Biased PSD estimate plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta = {'CSAcbf';'CSApp'};
f5 = figure(5);
set(f5,'DefaultAxesFontSize',24)
    SigindOfInterest = (u < uS+.175 & u > uS-.175);
    uSig = u(SigindOfInterest);
    RightOfSig = (u > uS+.175);
    uRofS = u(RightOfSig);
    LeftOfSig = (u < uS-.175);
    uLofS = u(LeftOfSig);
zsigsim = zeros(1,length(zeta));  
for aplot = 1
    
    clf
    hold on
    sigline = plot([uS uS],[-5 10*log10(2*L)].','-.m');
    zsigder = zeros(1,length(zeta));
    for plotind = 1:length(zeta)
    simSigNoisePSD = 10*log10(abs(simulated.(zeta{plotind}).SigNoisePSD{1,aplot}(SigindOfInterest)));
    simRoS = 10*log10(abs(simulated.(zeta{plotind}).SigNoisePSD{1,aplot}(RightOfSig)));
    simLoS = 10*log10(abs(simulated.(zeta{plotind}).SigNoisePSD{1,aplot}(LeftOfSig)));
    
    zsigder(plotind) = plot(u,10*log10(abs(derived.(zeta{plotind}).SigNoisePSD{1,aplot})),'color',plotcolor.(zeta{plotind}));
    
%     plot(u,10*log10(abs(simulated.(zeta{plotind}).SigNoisePSD{1,aplot})),':','color',plotcolor.sim)
    zsigsim(plotind) = plot(uSig(1:3:end),simSigNoisePSD(1:3:end),':d','MarkerSize',9,'LineWidth',2,'color',plotcolor.(zeta{plotind}));
    plot(uRofS(1:10:end),simRoS(1:10:end),':d','MarkerSize',9,'LineWidth',2,'color',plotcolor.(zeta{plotind}))
    plot(uLofS(1:10:end),simLoS(1:10:end),':d','MarkerSize',9,'LineWidth',2,'color',plotcolor.(zeta{plotind}))
    end
    tpsd = plot(u,10*log10(abs(derived.ULA.SigNoisePSD{1,aplot})),'-.','LineWidth',3,'color',plotcolor.truePSD);
    legend([sigline zsigder(1) zsigsim(1) zsigder(2) zsigsim(2) tpsd ],...
        {'Signal' 'CSA$_{\mathrm{cbf}}$' '' 'CSA$_{\mathrm{pp}}$' '' 'ULA'},'location','northeast')
%     legend('true PSD','Signal',zeta{1},'',zeta{2},'',zeta{3},'')
    xlabel('u'); ylabel('[dB]')
    title({'Signal arriving from $u_s=-0.4$: Fixed Aperture ($\approx$50$\lambda_o$)',['SNR$_{in}$ = ' num2str(10*log10(varS/varW)) 'dB, $\alpha$ = ' num2str(round(alpha(aplot)*1000)/1000) ', CSA of M = ' num2str(M) ', N = ' num2str(N) ', $\beta$ = ' num2str(beta) ]})
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Signal + Noise Probability/Cumulative Distribution Function plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotcolor.sim = [1 1 1];
SigNoisePDFanimation(length(alpha)) = struct('cdata',[],'colormap',[]);

if ~islocal
    f6 = figure(6);
    set(f6,'DefaultAxesFontSize',24);
for aplot = aplot
    clf
    hold on

    plot(derived.CSAcbf.xSigNoisePDF{1,aplot},derived.CSAcbf.SigNoisePDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.CSAcbf);
    plot(derived.ULA.xSigNoisePDF{1,aplot}(1:end),derived.ULA.SigNoisePDF{1,aplot}(1:end),'color',plotcolor.ULA);
    plot(simulated.ULA.xSigNoisePDF{1,aplot}(1:end),simulated.ULA.SigNoisePDF{1,aplot}(1:end),':d','color',plotcolor.ULA);

    
    plot(simulated.CSAcbf.xSigNoisePDF{1,aplot}(1:end-80),simulated.CSAcbf.SigNoisePDF{1,aplot}(1:end-80)...
        ,':+','MarkerSize',15,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xSigNoisePDF{1,aplot}(end-80:skip_data:end),simulated.CSAcbf.SigNoisePDF{1,aplot}(end-80:skip_data:end)...
        ,':+','MarkerSize',15,'color',plotcolor.CSAcbf);
    title(['Signal + Noise PDFs, L = ' num2str(L) ', $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ', $\alpha$ = ' num2str(round(alpha(aplot)*1000)/1000)])

    
    plot(derived.CSApp.xSigNoisePDF{1,aplot},abs(derived.CSApp.SigNoisePDF{1,aplot})...
        ,'LineWidth',5,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoisePDF{1,aplot}(1:end-80),simulated.CSApp.SigNoisePDF{1,aplot}(1:end-80)...
        ,':x','MarkerSize',18,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoisePDF{1,aplot}(end-80:skip_data:end),simulated.CSApp.SigNoisePDF{1,aplot}(end-80:skip_data:end)...
        ,':x','MarkerSize',18,'color',plotcolor.CSApp);
    
    xlabel('$z$')
    title('$H_1$')
    ylabel('Fixed No. Sensors Constraint')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Cumulative Distribution Function plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SigNoisePDFanimation(length(alpha)) = struct('cdata',[],'colormap',[]);
f7 = figure(7);
set(f7,'DefaultAxesFontSize',24);
for aplot = 1
    clf

    subplot(2,3,[2 5])
    hold on
    plot(derived.CSAcbf.xSigNoiseCDF{1,aplot},derived.CSAcbf.SigNoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xSigNoiseCDF{1,aplot}(1:200:end-200),...
        simulated.CSAcbf.SigNoiseCDF{1,aplot}(1:200:end-200),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xSigNoiseCDF{1,aplot}(end-200:50:end),...
        simulated.CSAcbf.SigNoiseCDF{1,aplot}(end-200:50:end),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xSigNoiseCDF{1,aplot}(end-60:15:end),...
        simulated.CSAcbf.SigNoiseCDF{1,aplot}(end-60:15:end),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    plot(simulated.CSAcbf.xSigNoiseCDF{1,aplot}(end-10:end),...
        simulated.CSAcbf.SigNoiseCDF{1,aplot}(end-10:end),':d','MarkerSize',9,'color',plotcolor.CSAcbf);
    title(['Signal + Noise CDFs, $\alpha$ = ' num2str(round(alpha(aplot)*1000)/1000)])
%     ylabel(['CBF CSA: \chi^2 Distributed\newline \newline \Lambda = ' num2str(Lmn)])
    xlabel('$z$')

    
    subplot(2,3,[1 4])
    hold on
    plot(derived.ULA.xSigNoiseCDF{1,aplot},derived.ULA.SigNoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.ULA);
    plot(simulated.ULA.xSigNoiseCDF{1,aplot}(1:200:end-200),...
        simulated.ULA.SigNoiseCDF{1,aplot}(1:200:end-200),':d','MarkerSize',9,'color',plotcolor.ULA);
    plot(simulated.ULA.xSigNoiseCDF{1,aplot}(end-200:60:end),...
        simulated.ULA.SigNoiseCDF{1,aplot}(end-200:60:end),':d','MarkerSize',9,'color',plotcolor.ULA);
    plot(simulated.ULA.xSigNoiseCDF{1,aplot}(end-60:15:end),...
        simulated.ULA.SigNoiseCDF{1,aplot}(end-60:15:end),':d','MarkerSize',9,'color',plotcolor.ULA);
    plot(simulated.ULA.xSigNoiseCDF{1,aplot}(end-10:end),...
        simulated.ULA.SigNoiseCDF{1,aplot}(end-10:end),':d','MarkerSize',9,'color',plotcolor.ULA);

    
    subplot(2,3,[3 6])
    hold on
    plot(derived.CSApp.xSigNoiseCDF{1,aplot},derived.CSApp.SigNoiseCDF{1,aplot}...
        ,'LineWidth',4,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoiseCDF{1,aplot}(1:200:end-200),...
        simulated.CSApp.SigNoiseCDF{1,aplot}(1:200:end-200),':d','MarkerSize',9,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoiseCDF{1,aplot}(end-200:50:end),...
        simulated.CSApp.SigNoiseCDF{1,aplot}(end-200:50:end),':d','MarkerSize',9,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoiseCDF{1,aplot}(end-60:15:end),...
        simulated.CSApp.SigNoiseCDF{1,aplot}(end-60:15:end),':d','MarkerSize',9,'color',plotcolor.CSApp);
    plot(simulated.CSApp.xSigNoiseCDF{1,aplot}(end-10:end),...
        simulated.CSApp.SigNoiseCDF{1,aplot}(end-10:end),':d','MarkerSize',9,'color',plotcolor.CSApp);

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    rho vs alpha plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load  rho.mat
f10 = figure(10);
set(f10,'DefaultAxesFontSize',24)
for b = 1:1
    clf
    semilogx(1-alpha,abs(derived.CSApp.rhoNoise(b,:)),'r','LineWidth',4)
    hold on
    semilogx(1-alpha(1:2:end),abs(simulated.CSApp.rhoNoise(1:2:end)),':dr','MarkerSize',9)
%     semilogx(1-alpha,(rhoSigNoiseder(b,:)),'LineWidth',4,'color',[0 0 1])
%     semilogx(1-alpha(1:2:end),abs(rhoSigNoisesim(1:2:end)),':d','MarkerSize',9,'color',[0 0 1])
    semilogx(1-alpha,abs(derived.CSApp.rhoSigNoise(b,:)),'LineWidth',4,'color',[0 .6 .1])
    semilogx(1-alpha(1:2:end),abs(simulated.CSApp.rhoSigNoise(1:2:end)),':d','MarkerSize',9,'color',[0 .6 .1])
    grid on
    ax = gca;
    xtick = linspace(.1,1,10);
    set(ax,'XTick',xtick);
    xticklabel = abs(round(xtick.*100)./100 - 1);
    set(ax,'XTickLabel',xticklabel);
    set(ax,'Xdir','reverse')
%     ax.XTickLabelRotation=45;
    xlim([0.01 1]); ylim([0 1.1]);
    xlabel('$\alpha$')
    ylabel('$\rho$')
    legend('$\rho_{0}$','','$\rho_{1}$','','Location','SouthEast')
    title(['Correlation Coefficient $\rho$ vs. $\alpha$, M = ' num2str(M) ', N = ' num2str(N)])
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    KSTest plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f11 = figure(11);
    set(f11,'DefaultAxesFontSize',24)
    clf
for dummy = 1
%     subplot(4,3,1)
%     stem(alpha,derived.M.NoiseKSTest,'color',plotcolor.M);
%     subplot(4,3,2)
%     stem(alpha,derived.N.NoiseKSTest,'color',plotcolor.N);
   
    subplot(2,3,1)
    stem(alpha,derived.ULA.NoiseKSTest,'MarkerSize',10,'color',plotcolor.ULA);
    title('ULA')
    subplot(2,3,2)
    stem(alpha,derived.CSAcbf.NoiseKSTest,'MarkerSize',10,'color',plotcolor.CSAcbf);
    xlabel('Noise Only KS Test (above)')
    title('CSA$_{\mathrm{cbf}}$')
    subplot(2,3,3)
    stem(alpha,derived.CSApp.NoiseKSTest,'MarkerSize',10,'color',plotcolor.CSApp);
    title('CSA$_{\mathrm{pp}}$')

%     subplot(4,3,7)
%     stem(alpha,derived.M.SigNoiseKSTest,'color',plotcolor.M);
%     subplot(4,3,8)
%     stem(alpha,derived.N.SigNoiseKSTest,'color',plotcolor.N);
   
    subplot(2,3,5)
    stem(alpha,derived.CSAcbf.SigNoiseKSTest,'MarkerSize',10,'color',plotcolor.CSAcbf);
    title('Signal + Noise KS Test (below)');xlabel('$\alpha$')
     subplot(2,3,4)
    stem(alpha,derived.ULA.SigNoiseKSTest,'MarkerSize',10,'color',plotcolor.ULA);
    subplot(2,3,6)
    stem(alpha,derived.CSApp.SigNoiseKSTest,'MarkerSize',10,'color',plotcolor.CSApp);
end

end
%% ROC plots
if doROC
f12 = figure(12);
set(f12,'DefaultAxesFontSize',24);
    for dummy = 1
        ROC1 = 7;
        skip1 = 3;
        endskip = 8;
        clf
        hold on
        plot(log10(abs(derived.ULA.Pfa{1,ROC1})),log10(abs(derived.ULA.Pd{1,ROC1})),'color',plotcolor.ULA)
        plot(log10(simulated.ULA.Pfa{1,ROC1}(1:skip1:end-endskip)),log10(simulated.ULA.Pd{1,ROC1}(1:skip1:end-endskip)),':s','MarkerSize',10,'color',plotcolor.ULA)
        plot(log10(abs(derived.CSAcbf.Pfa{1,ROC1})),log10(abs(derived.CSAcbf.Pd{1,ROC1})),'color',plotcolor.CSAcbf)
        plot(log10(simulated.CSAcbf.Pfa{1,ROC1}(1:skip1:end-endskip)),log10(simulated.CSAcbf.Pd{1,ROC1}(1:skip1:end-endskip)),':s','MarkerSize',10,'color',plotcolor.CSAcbf)
        plot(log10(abs(derived.CSApp.Pfa{1,ROC1})),log10(abs(derived.CSApp.Pd{1,ROC1})),'color',plotcolor.CSApp)
        plot(log10(simulated.CSApp.Pfa{1,ROC1}(1:skip1:end-endskip)),log10(simulated.CSApp.Pd{1,ROC1}(1:skip1:end-endskip)),':s','MarkerSize',10,'color',plotcolor.CSApp)
          
        legend('ULA','','CSA$_{\mathrm{cbf}}$','','CSA$_{\mathrm{pp}}$','','location','southeast')
        ylabel('log$_{10}(P_d)$')
        xlabel('log$_{10}(P_{fa})$')
        title(['Monte-Carlo ROC for $\alpha$ = ' num2str(alpha(ROC1)) ', SNR$_{in}=$ ' num2str(10*log10(varS/varW)) 'dB, L = ' num2str(L) ', $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ])
        xlim([-3.5 0])
        ylim([log10(min(simulated.CSAcbf.Pd{1,ROC1}))-.01 0])
%         axis tight
    end
end



% %% ROC plots
% if doROC
%     load ULAROCs.mat
% f12 = figure(12);
% set(f12,'DefaultAxesFontSize',24);
%     for dummy = 1
%         ROC1 = 2;
%         ROC2 = 5;
%         ROC3 = 7;
%         skip1 = 3;
%         endskip = 8;
%         clf
%         hold on
%         plot(log10(derived.ULA.Pfa{1,ROC1}(1:5:end)),log10(derived.ULA.Pd{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
% %         plot(log10(ULAPfaDerived{1,ROC1}(1:5:end)),log10(ULAPdDerived{1,ROC1}(1:5:end)),':o','MarkerSize',10,'color',plotcolor.ULA)
%         plot(log10(derived.CSAcbf.Pfa{1,ROC1}(1:30:end)),log10(derived.CSAcbf.Pd{1,ROC1}(1:30:end)),':o','MarkerSize',10,'color',plotcolor.CSAcbf)
%         plot(log10(derived.CSApp.Pfa{1,ROC1}(1:40:end)),log10(derived.CSApp.Pd{1,ROC1}(1:40:end)),':o','MarkerSize',10,'color',plotcolor.CSApp)
%         
% %         plot(log10(derived.ULA.Pfa{1,ROC2}(1:8:end)),log10(derived.ULA.Pd{1,ROC2}(1:8:end)),':d','MarkerSize',10,'color',plotcolor.ULA)
% %         plot(log10(derived.CSAcbf.Pfa{1,ROC2}(1:15:end)),log10(derived.CSAcbf.Pd{1,ROC2}(1:15:end)),':d','MarkerSize',10,'color',plotcolor.CSAcbf)
% %         plot(log10(derived.CSApp.Pfa{1,ROC2}(1:60:end)),log10(derived.CSApp.Pd{1,ROC2}(1:60:end)),':d','MarkerSize',10,'color',plotcolor.CSApp)
% %         
% %         plot(log10(derived.ULA.Pfa{1,ROC3}(1:20:end)),log10(derived.ULA.Pd{1,ROC3}(1:20:end)),':s','MarkerSize',10,'color',plotcolor.ULA)
% %         plot(log10(derived.CSAcbf.Pfa{1,ROC3}(1:50:end)),log10(derived.CSAcbf.Pd{1,ROC3}(1:50:end)),':s','MarkerSize',10,'color',plotcolor.CSAcbf)
% %         plot(log10(derived.CSApp.Pfa{1,ROC3}(1:100:end)),log10(derived.CSApp.Pd{1,ROC3}(1:100:end)),':s','MarkerSize',10,'color',plotcolor.CSApp)
%         
%         legend('',['$\alpha$ = ' num2str(alpha(ROC1))],'','',['$\alpha$ = ' num2str(alpha(ROC2))],'','',['$\alpha$ = ' num2str(alpha(ROC3))],'','location','southeast')
%         ylabel('log$_{10}(P_d)$')
%         xlabel('log$_{10}(P_{fa})$')
%         title({'ROCs: ULA (red), CSA$_{\mathrm{cbf}}$ (orange), and CSA$_{\mathrm{pp}}$ (black)',['SNR$_{in}=$ ' num2str(10*log10(varS/varW)) 'dB, L = ' num2str(L) ', $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ]})
% %         xlim([log10(min(simulated.ULA.Pfa{1,ROC1}))-.1 0])
%         title({['SNR$_{in}=$ ' num2str(10*log10(varS/varW)) 'dB ROCs: ULA (red), CSA$_{\mathrm{cbf}}$ (orange), and CSA$_{\mathrm{pp}}$ (black)'],['ULA Aperture = 50$\lambda_o$, CSA Aperture = 150$\lambda_o$, $\beta$ = ' num2str(beta) ', M = ' num2str(M) ', N = ' num2str(N) ]})
% 
% %         ylim([log10(min(simulated.CSAcbf.Pd{1,ROC1}))-.01 0])
%         ylim([-3.5 0])
%         xlim([-150 0])
%         axis tight
%     end
% end