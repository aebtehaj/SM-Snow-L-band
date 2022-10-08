function epsr = soil_perm_MBSDM_Mironov(mv,clayfrac,fGHz)
% calculate soil permittivity using the mineralogically based spectroscopic
% dielectric model (MBSDM) for Moist Soils using the generalized refractive
% mixing dielectric moel (GRMDM).
% Scalar Inputs:
%   mv: soil moisture, (volumetric fraction)
%   clayfrac: clay content weight fraction
%   fGHz:  wave frequency in GHz
% Outputs
%   epsr: soil complex dielectric constant (CDC)
% Reference:
%   Mironov et al. TGRS 47(7):2059-2070, 2009
% 

C = clayfrac;
f = fGHz*1e9;
omega = 2*pi*f;

% fit of MBSDM spectroscopic parameters, correlated to C
nd = 1.634 - 0.539*C + 0.2748*C^2;
kd = 0.03952 - 0.04038*C;
mvt = 0.02863 + 0.30673*C;
e0b = 79.8 - 85.4*C + 32.7*C^2;
tb = 1.062e-11 + 3.450e-12*C;
sb = 0.3112 + 0.467*C;
su = 0.3631 + 1.217*C;
e0u = 100;
tu = 8.5e-12;

% Debye relaxation equations for bound and free water components
einf = 4.9;
eps0 = 8.854e-12;
eb = einf + (e0b - einf)/(1 - 1i*omega*tb) + 1i*sb/omega/eps0;
eu = einf + (e0u - einf)/(1 - 1i*omega*tu) + 1i*su/omega/eps0;

nbc = sqrt(eb); nb = real(nbc); kb = imag(nbc);
nuc = sqrt(eu); nu = real(nuc); ku = imag(nuc);

% refractive mixing dielectric model
if mv < mvt
    nm = nd + (nb - 1)*mv;
    km = kd + kb*mv;
else
    nm = nd + (nb - 1)*mvt + (nu - 1)*(mv - mvt);
    km = kd + kb*mvt + ku*(mv - mvt);
end

emr = nm^2 - km^2;
emi = 2*nm*km;

epsr = complex(emr,emi);
end