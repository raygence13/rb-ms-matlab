function PDF = CSAproductPDF(z,sigmaM,sigmaN,rho)
rho = abs(rho);
x = z./(sigmaM*sigmaN*(1-rho^2));
BesselIterm = besseli(0,x*abs(rho));
BesselKterm = besselk(0,x);

PDF = x./(sigmaM*sigmaN).*BesselIterm.*BesselKterm;

isPDFNaN = isnan(PDF);
isPDFzero = PDF==0;
ABind = or(isPDFNaN,isPDFzero);

if any(ABind)
%     newBesselterm = exp(x(ABind).*(rho-1))./(2.*x(ABind).*sqrt(rho));
    PDF(ABind) = exp(x(ABind).*(rho-1))./(2*sigmaM*sigmaN*sqrt(rho));
end
end

% kterms = 100;
% 
% BesselIgammaseries = zeros(length(z),kterms);
% BesselKgammaseries = zeros(length(z),kterms);
% sumBesselIgammaseries = zeros(length(z),1);
% sumBesselKgammaseries = zeros(length(z),1);
% for zind = 1:length(z)
%     for k = 0:kterms
%         BesselIgammaseries(zind,k+1) = (-1)^k./((2.*z(zind)*rho).^k).*gamma(k+1/2)./factorial(k)/gamma(1/2-k);
%         BesselKgammaseries(zind,k+1) = 1./((2.*z(zind)).^k).*gamma(k+1/2)./factorial(k)/gamma(1/2-k);
%     end
%     sumBesselIgammaseries(zind,1) = sum(BesselIgammaseries(zind,:),2);
%     sumBesselKgammaseries(zind,1) = sum(BesselKgammaseries(zind,:),2);
% end
% 
% y = scaleterm.*sumBesselIgammaseries.*sumBesselKgammaseries;