function [TL] = TransmissionLoss(f, R, flag);
% function [TL] = TransmissionLoss(f, R)
% Function computes the transmission loss in dB for specified frequencies
% out to a maximum range.  This function uses the frequency attenuation
% constant due to Thorp (see Urick, page 108, 3rd edition)
% INPUTS:
% f = [fmin : delta f : fmax], a vector of frequencies in KHz, 1xNf
% R = [Rmin : delta R : Rmax], a vector of Ranges in Yds; 1xNr
% flag; If flag ==  's', use Spherical spreading, 'c', cylindrical
% OUTPUTS:
% TL in db as function of frequency and range, NfxNr

% F. Khan - Code 1513

f2 = f .^ 2;
a1 = 0.1 * f2 ./ (1 + f2);
a2 = 40 * f2 ./ (4100 + f2);
a3 = .000275 * f2 + 0.003;
a = a1 + a2 + a3;   % attenuation coefficient, dB/Kyd

if flag == 's'
    g = 20;
end
if flag == 'c'
    g = 10;
end

TL = [];
for k = 1 : length(f),
    %TL(k, :) = (20 *log10(R) + a(k) * R*1e-3);
    TL(k, :) = (g *log10(R) + a(k) * R*1e-3);
end

% minTL = 10*floor(min(min(TL))/10);
% figure
% %set(1,'Position', [15 490 825 500])
% plot(R, TL, 'LineWidth', 1.5);grid
% xlabel('Range Yds')
% ylabel('Transmission Loss dB')
% axis([0 max(R) minTL -40])
% title('Transmisson Loss, dB - see Urick, p108, 3rd Ed')
% 
% figure
% imagesc(R/1000, f/8, TL, [60 200]);colorbar
% ylabel('f/f{_0}','Fontsize',12)
% xlabel('Range Kyds','Fontsize',12)
% axis([.95 20.05 .1 3.25])
% set(gca,'ytick',[ .125 .5:.5:3 3.25])
% set(gca,'xtick',[1 5 10 15 20])

% min(min(TL))
% max(max(TL))

return