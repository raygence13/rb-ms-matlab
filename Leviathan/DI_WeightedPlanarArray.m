function [DI] = DI_WeightedPlanarArray(N, M, f, c, dx, dy, phi_s, W);
%
% [DI] = DI_WeightedPlanarArray(N, M, f, c, dx, dy, phi_s, W);
% This function computes the Directivity Index of a plane rectangular array.
% It uses the approximation given in "Approximations to Directivity for 
% Linear, Planar and Volumetric Apertures and Arrays," A. Nuttall, B. Cray,
% NUWC-NPT Technical Report 10,798.
%
% INPUTS:
% N == Number of elements along one dimension in array
% M == Number of elements along other dimension in array
% f == Frequencies of interest (vector input if desired), Hz
% c == Sound speed, feet/sec
% dx == Element spacing along one dimension, units consistent with element
%       spacing dy and sound speed c
% dy == Element spacing along other dimension, units consistent with element
%       spacing dx and sound speed c
% phi_s == Steering direction, degrees;  phi_s = 0 is broadside
% W == N x M array of weights
% OUTPUTS:
% DI == Directivity Index for a Baffled, Weighted planar array, dB
% Note: For unbaffled array, subtract 3dB from DI

% F. Khan NUWCNPT COde 1513

lambda = c ./ f;
n0 = (N+1)/2;
m0 = (M+1)/2;
sum_n = [];
for n = 1 : N-1;
    sum_n(n) = sum(W(n, :) *(n - n0)^2);
end
sum_m = [];
for m = 1 : M-1;
    sum_m(m) = sum(W(:, m) *(m - m0)^2);
end
eff_Area = dx * dy * sqrt(sum(sum_n)) * sqrt(sum(sum_m)) / sum(sum(W));
DI = 2 * 8 * pi * pi * cosd(phi_s) * eff_Area ./(lambda.^2);
DI = 10 * log10(DI);

% figure
% plot(f/1000, DI, '-x', 'Linewidth', 1.5);grid
% xlabel('Frequency, KHz')
% ylabel('DI, dB')
% title('Directivity Index For Planar Array')