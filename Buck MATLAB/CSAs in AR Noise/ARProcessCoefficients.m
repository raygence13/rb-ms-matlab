function [ b,a,extrasamples ] = ARProcessCoefficients( alpha,USF,threshold )
%[ b,a,extrasamples ] = ARProcessCoefficients28Feb2016( alpha,USF,threshold )
%   Function to generate LCCDE to model colored noise as AR system based on
%   undersampling factor in a ULA
%   y[n] - alpha^USF*y[n-USF] = sqrt(1-alpha^(2*USF))*x[n]

%   Inputs
% alpha is correlation coefficient
% USF is undersampling factor
% threshold is used to set extrasamples

%   Outputs
% b is the coefficient for the x[n] terms normalized to achieve varW noise
% a is a vector for the y[n] terms
% extrasamples is the number of extra samples to oversample the process and
% reduce transient effects
%%% alpha^(extrasamples) = 10^(threshold) sets value for extrasamples

% b = sqrt(1 - alpha^(2*USF));
b = 1;
a = [1 zeros(1,USF-1) -alpha^(USF)];
extrasamples = ceil(threshold/log10(alpha));

end

