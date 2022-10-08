function [svv,shh] = SPM3D(kinc, thetai, ktran, h, ratio)
% calculate rough surface backscattering by 1st order SPM.
% h is rms height in meter
% ratio is correlation length / rms height
% thetai is the angle in the incidence medium
% 

l = ratio*h;
k = kinc;
k1 = ktran;
% k1 = k.*sqrt(epsr_soil);
kzi = k.*cos(thetai);
kxi = k.*sin(thetai);
k1zi = sqrt(k1.^2 - kxi.^2);

kz = k.*cos(thetai);
kx = -k.*sin(thetai);
k1z = sqrt(k1.^2-kx.^2);

termh = (abs((kzi - k1zi)./(kzi + k1zi))).^2;
termv = (abs((k1.^2 - k.^2).*(k1.^2.*k.^2.*(sin(thetai)).^2 ...
         + k.^2.*k1z.*k1z)./(k1.^2.*kz+k.^2.*k1z).^2)).^2;

% krms = k.*h;
shh = 8.*k.^4.*h.^2.*l.^2./(4.*k.^2.*l.^2.*(sin(thetai)).^2+1).^(3./2).*(cos(thetai)).^4.*termh;
svv = 8.*k.^4.*h.^2.*l.^2./(4.*k.^2.*l.^2.*(sin(thetai)).^2+1).^(3./2).*(cos(thetai)).^4.*termv;
% cross pol - OH model
% shv = 0.23.*sqrt((sqrt(real(epsr_soil))-1)./(sqrt(real(epsr_soil))+1)).*(1-exp(-krms)).*svv; 

end