% Emission Model
% Divya Kumawat, 07/2023
%% Computes brightness temperature at Vertical and Horizontal Polarization for land-snow-vegetation system at L-band (1.4 GHz)

% Inputs
% Tcanopy: Temperature of vegetation (K)
% tau: Vegetation Optical Depth (VOD)
% omega: Vegetation Single scattering albedo
% h: Soil Roughness parameter (Q/H model)
% freq: Frequency (GHz)
% obs_angle: Observation angle (deg)
% Tsoil: Temperature of ground (K)
% rho_s: density of snowpack (Kg/m^3)
% espr_ground: dielectric constant of ground at L-band

% Outputs
% Tb: Vector of brightness temperature at Horizontal and Vertical Polarization

%==========================================================================

%Example
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% % % Inputs
% % % ----------------------------------------------------------------------
% % Tcanopy =   265;
% % tau = 0.1;
% % omega = 0.07;
% % h = 0.15;
% % freq= 1.4;
% % obs_angle = 40;
% % Tsoil = 265;
% % rho_s = 250;
% % epsr_ground = 5 -0.5i;
% % 
% % Tb= TO_snow(Tcanopy,tau,omega,h,freq,obs_angle,Tsoil,rho_s,epsr_ground);

function Tb = TO_snow(Tcanopy,tau,omega,h,freq,obs_angle,Tsoil,rho_s,espr_ground)

% % % Vegetation Microwave Transmissivity----------------------------------
angler = obs_angle.*pi./180 ;
exptau= exp(-real(tau)./cos(angler));

% % % Dielectric constant of snow----------------------------------------
Tsnow = Tcanopy;
T = Tsnow-273.15; Ps = rho_s/1000; f = 1.4;
epsr_ds = RelDielConst_DrySnow(T, Ps, f);
theta1 = deg2rad(obs_angle);

% % % Computing effective reflectivity and transmissivity of soil-snow system-------------------------------
eps1 = 1; eps2 = epsr_ds; eps3 = espr_ground; f = freq;
lam2 = 1i *20*pi/3*f*sqrt(eps2); % wave number in medium 2

theta2 = acos((1-(sqrt(eps1)/sqrt(eps2).*sin(theta1)).^2).^(1/2));
theta3 = acos((1-(sqrt(eps1)/sqrt(eps3).*sin(theta1)).^2).^(1/2));

rho12h = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos(theta2)) ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos(theta2));
rho12v = (sqrt(eps1).*cos(theta2)-sqrt(eps2).*cos(theta1)) ./ (sqrt(eps1).*cos(theta2) + sqrt(eps2)*cos(theta1));

rho23h = sqrt(exp(-h.*cos(theta2).*cos(theta2))).*((sqrt(eps2).*cos(theta2)-sqrt(eps3).*cos(theta3)) ./ (sqrt(eps2).*cos(theta2) + sqrt(eps3).*cos(theta3)));
rho23v = sqrt(exp(-h.*cos(theta2).*cos(theta2))).*((sqrt(eps2).*cos(theta3)-sqrt(eps3).*cos(theta2)) ./ (sqrt(eps2).*cos(theta3) + sqrt(eps3)*cos(theta2)));

a1 = rho12h; b1 = rho23h; c1 = -2*lam2.*cos(theta2);
gammah = 1+ (sign(a1*abs(b1)-1)*((-1+a1^2+abs(b1)^2-a1^2*abs(b1)^2)/(a1^2*abs(b1)^2-1)));

a = rho12v; b = rho23v; c = -2*lam2.*cos(theta2);
gammav = 1+ (sign(a*abs(b)-1)*((-1+a^2+abs(b)^2-a^2*abs(b)^2)/(a^2*abs(b)^2-1)));

% % % Computing Brightness temperature
TBv  =  Tsoil*(1-gammav).*exptau + Tcanopy.*(1-omega).*(1-exptau).*gammav.*exptau + Tcanopy.*(1-omega).*(1-exptau);
TBh = Tsoil*(1-gammah).*exptau + Tcanopy.*(1-omega).*(1-exptau).*gammah.*exptau + Tcanopy.*(1-omega).*(1-exptau);

Tb=[TBh;TBv];

end