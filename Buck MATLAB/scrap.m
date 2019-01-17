clearvars
close all
clc
present(0)
profile on

%%
L = 11; % Number of Sensors
SampleSize = 3*L;

% Creating Signal of Interest (SoI)
varS    = 10;    % variance of SoI
uS      = 0;    % SoI direction of arrival (DoA)
vS      = exp(1i*pi*(0:L-1)'*uS);  % SoI replica vector
s       = sqrt(varS/2)*randn(1,SampleSize)...   % SoI Samples, [1 x SampleSize]
          + 1i*sqrt(varS/2)*randn(1,SampleSize);
SoI_Data    = vS*s;     % SoI Data, [L x SampleSize]
Ens_SoI     = varS*(vS*vS');
R_SoI       = (SoI_Data*SoI_Data')/SampleSize;  % SoI Sample Covariance Matrix (SCM)

% Creating Interferer (Int)
varI    = varS*100; % variance of Int
uI      = 3/L;      % Int DoA
vI      = exp(1i*pi*(0:L-1)'*uI);  % Int replica vector
Int     = sqrt(varI/2)*randn(1,SampleSize)...   % Int Samples, [1 x SampleSize]
          + 1i*sqrt(varI/2)*randn(1,SampleSize);
Int_Data    = vI*Int;     % Int Data, [L x SampleSize]
Ens_Int     = varI*(vI*vI');
R_Int       = (Int_Data*Int_Data')/SampleSize;  % Int SCM

% Creating uncorrelated Gaussian Noise
varW = 10;
Noise_Data = sqrt(varW/2)*randn(L,SampleSize)...
    + 1i*sqrt(varW/2)*randn(L,SampleSize);
Ens_Noise = varW*eye(L);
R_Noise = (Noise_Data*Noise_Data')/SampleSize;

x = SoI_Data + Int_Data + Noise_Data;   % Input vector
S = (x*x')/SampleSize;      % Data SCM

u0      = 0;    % Look direction
v0      = exp(1i*pi*(0:L-1)'*u0);	% manifold vector
DL      = varW*eye(L);
EMI     = (Ens_SoI + Ens_Int + Ens_Noise)\v0;
Ens_wMVDR   = EMI/(v0.'*EMI);
SMI     = (R_Int+R_Noise+DL)\v0;
wMVDR = SMI/(v0.'*SMI);

uNfft = 2048;
ur = linspace(0,1,uNfft/2);
ul = linspace(-1,0-1/uNfft,uNfft/2);
u = [ul ur];
B_MVDR = fftshift(fft(Ens_wMVDR,uNfft));
SB_MVDR = fftshift(fft(wMVDR,uNfft));

subplot(211)
stem(uI,10*log10(varI));
hold on
stem(uS,10*log10(varS))
xlim([-1 1])

subplot(212)
plot(u,10*log10(abs(B_MVDR)))
hold on
plot(u,10*log10(abs(SB_MVDR)),'k')
