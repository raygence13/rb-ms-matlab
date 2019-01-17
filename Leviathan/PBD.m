function [ SL_T,NL_T,TW ] = PBD( f,f_u,f_l,df,Signal_Spectrum, Noise_Spectrum )
% Passive Broadband Detect function
% Inputs
% f: frequency vector [Hz]
% f_u: upper corner frequency value of filter [Hz]
% f_l: lower corner frequency value of filter [Hz]
% df: delta f for integration [Hz]
% Signal_Spectrum: signal level values across frequency (linear)
% Noise_Spectrum: noise level values across frequency (linear)

% Outputs
% SL_T: Signal power in band between f_u and f_l [dB]
% NL_T: Noise power in band between f_u and f_l [dB]
% TW: Effective noise bandwidth

if f_l > f_u
    f_l = temp;
    f_l = f_u;
    f_u = temp;
end


if abs(length(f)-length(Signal_Spectrum))>0
    warning('Frequency vector and signal spectrum are not of the same length')
elseif abs(length(f)-length(Noise_Spectrum))>0
    warning('Frequency vector and noise spectrum are not of the same length')
end

if f_l < f(1)
    warning(['Lower frequency not within frequency vector. Changing f_l to ' num2str(f(1))])
    f_l = f(1);
elseif f_u > f(end)
    warning(['Upper frequency not within frequency vector. Changing f_u to ' num2str(f(end))])
end

if abs(f_u-f_l) < df
    warning('chosen df is greater than the band of integration')
    df = abs(f_u-f_l);
end
fl_index = find(abs(f-f_l) < .001);
fu_index = find(abs(f-f_u) <.001);


num_int = floor((f_u-f_l)/df);  % number of bands to integrate
SL_of_interest = Signal_Spectrum(fl_index:fu_index);
NL_of_interest = Noise_Spectrum(fl_index:fu_index);

SLs = zeros(1,num_int);
NLs = zeros(1,num_int);
base_mid = round(length(SL_of_interest)/(num_int+1));   % base center frequency index
for ii = 1:num_int
   SLs(ii) = SL_of_interest(base_mid*ii);
   NLs(ii) = NL_of_interest(base_mid*ii);
end

SL_T = 10*log10(df) + 10*log10(sum(SLs));
NL_T = 10*log10(df) + 10*log10(sum(NLs));
TW = df*sum(NLs)^2/sum(NLs.^2);
% 
% % SSI = sum(Signal_Spectrum(fl_index:fu_index));
% SSI = trapz(f(fl_index:fu_index),Signal_Spectrum(fl_index:fu_index));
% SL_T = 10*log10(df) + 10*log10(SSI);
% 
% % Ai = 104.312;
% % NSI = sum(Noise_Spectrum(fl_index:fu_index)/Ai);
% % NSI2 = sum((Noise_Spectrum(fl_index:fu_index)/Ai).^2);
% 
% NSI = sum(Noise_Spectrum(fl_index:fu_index));
% NSI2 = sum((Noise_Spectrum(fl_index:fu_index)).^2);
% 
% % NSI = trapz(f(fl_index:fu_index),Noise_Spectrum(fl_index:fu_index)/Ai);
% % NSI2 = trapz(f(fl_index:fu_index),(Noise_Spectrum(fl_index:fu_index)/Ai).^2);
% 
% NL_T = 10*log10(df*NSI);
% TW = df*(NSI*sqrt(Ai))^2/(NSI2*Ai^2);
% TW = df*(NSI)^2/(NSI2);
end

