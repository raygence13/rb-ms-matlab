% Example in Burdic, p418
clear all
close all
present(0)

% Use example to create signal and noise spectrum for testing broadband
f =  [100   150   200 250   300   350   400   450   500    550   650   750   850   950   1000]; % [Hz]
SL = [106.5 106.8 107 107.2 107.4 107.5 107.6 107.6 107.45 107.3 106.8 106.2 105.7 105.2 105];  % [dB]
%SL = [105.2 106 106.55 107 107.3 107.45 107.5 107.5 107.45 107.3 106.8 106.2 105.7 105.2 105];
NL = [73.3  71.5 70  68  66  64  62  60  58  55  48  42  38  35  34.2]; % [dB]
sl = 10 .^(SL/10);  % linear
nl = 10 .^(NL/10);  % linear
f_i = 100 : 1000;   % [Hz]
sl_i = interp1(f, sl, f_i); 
nl_i = interp1(f, nl, f_i, 'cubic');
SL_i = 10*log10(sl_i);
NL_i = 10*log10(nl_i);

subplot(121)
semilogx(f_i, NL_i, 'r',  f, NL, 'rx', 'Linewidth', 1.5);
grid on
axis([100 1000 30 90])
ylabel('Spectrum Level db re uPa per Hz')
xlabel('Frequency Hz')
title('Noise Spectrum Level')

subplot(122)
semilogx(f_i, SL_i, 'b',f, SL, 'bx', 'Linewidth', 1.5)
grid on
axis([100 1000 100 120])
title('Signal Spectrum Level')
xlabel('Frequency Hz')


f_u = 1000;
f_l = 500;
df = 10;
[SL_T,NL_T,TW] = PBD( f_i,f_u,f_l,df,sl_i, nl_i )