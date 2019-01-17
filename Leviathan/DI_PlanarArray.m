function [DI,varargout] = DI_PlanarArray(N, M, f, f0, c, phi_s,flag)
%
% [DI] = DI_PlanarArray(N, M, f, c, dx, dy, phi_s);
% This function computes the Directivity Index of a plane rectangular array.
% It uses the approximation given in Approximations to Directivity for 
% Linear, Planar and Volumetric Apertures and Arrays, A. Nuttall, B. Cray,
% NUWC-NPT Technical Report 10,798.
%
% INPUTS:
% flag == 'S' for fixed sensors constraint N (x dir) & M (y dir) 
% N == Number of elements along one dimension in array
% M == Number of elements along other dimension in array
% flag == 'A' for fixed aperture constraint, N & M are in ft.
% N == Length [ft] along x dimension in array
% M == Length [ft] along y dimension in array
% f == Frequencies of interest (vector input if desired), Hz
% f0 == Design frequency of array
% c == Sound speed, feet/sec
% phi_s == Steering direction, degrees;  phi_s = 0 is broadside

% OUTPUTS:
% DI == Directivity Index for a BAFFLED planar array, dB
% For unbaffled array, subtract 3dB from DI

% F. Khan NUWCNPT Code 1513
lambda = c ./ f;    % (N1 * DX)*(M1 * dy)=Lx * Ly is the effective area
dx = c/f0./2;
dy = c/f0./2;

if strcmp(flag,'A')
    varargout{1} = floor(N./dx);   % Number of sensors in x direction
    varargout{2} = floor(M./dy);   % Number of sensors in y direction
    N1 = sqrt(varargout{1}.^2-1);
    M1 = sqrt(varargout{2}.^2-1);
elseif strcmp(flag,'S')
    N1 = sqrt(N.^2-1);
    M1 = sqrt(M.^2-1);
end

DI = (pi/3) * 4 * pi * cosd(phi_s) * dx * dy * N1 * M1 ./(lambda.^2);
DI = 10 * log10(DI).';
DI((DI <=0)) = 0;

% figure
% plot(f/1000, DI, '-x', 'Linewidth', 1.5);grid
% xlabel('Frequency, KHz')
% ylabel('DI, dB')
% title('Directivity Index For Planar Array')