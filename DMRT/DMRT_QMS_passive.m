function [TBv,TBh,deg0,ot,albedo,epsr_snow] = DMRT_QMS_passive(fGHz,dia,rho,tau,depth,Tsnow,Tg,epsr_ground,rough)
% passive remote sensing of layered snowpack
% solve passive DMRT equation using Discrete Ordinate-Eigenanalysis
% ( Guass-Legendre quadrature ) method, assuming QCA Mie Sticky
% spheres (QMS) scattering model
% bottom boundary could be rough, reflectivity described by Q/H, or
% Wegmuller and Matzler model (1999).

%
% Input parameters:
%   fGHz - frequncy in GHz, scalar
%   dia  - snow grain diameter, 1D array
%   rho  - snow density in gm/cc, 1D array
%   tau  - stickiness parameter in QCA model, 1D array
%   depth - snow depth of each layer in centimeter, 1D array, from top to bottom
%   Tsnow - snow temperature in Kelvin, 1D array
%   Tg    - ground temperature in Kelvin
%   epsr_ground -ground permittivity, could be complex
%   rough - specify the roughness
%           'QH' model: Q, H
%           'WM' (Wegmuller and Matzler 1999): s (rms height)
%
% Output:
%   TBv, TBh - array of Brightness Temperature in vertical and horizontal
%       polrization
%   deg0  - sampling angles of TB in air, in degree.
%   ot,albedo,epsr_snow: optical thickness, scattering albedo and snow
%       effective permittivity of each layer
%
% Ref:  Liang et al., TGRS, 46(11): 3663-3671, 2008
% Copyright: University of Washington, Electrical Engineering Department
% Revision Date: Aug. 1, 2014, Sep. 09, 2014.
%

% quadrature angles and integration weight
% global wi xi Nquad ndeg
Nquad = 32;
ndeg = Nquad/2;

[xi,wi] = GLNodeWt(Nquad);

mu0 = -xi(1:ndeg);
wi0 = wi(1:ndeg);

% k = 2*pi*fGHz/30; % 1/cm
eps0_ = 8.854e-12;
mu0_ = pi*4e-7;
c = 1/sqrt(eps0_*mu0_);
wave = c/fGHz*1e-7; % cm
k = 2*pi/wave;      % 1/cm

rho_ice = 0.917; % ice density in gm/cc
nlayer = length(depth);

% Calculate phase matrix from QCA
ke_a = zeros(1,nlayer);
Fm_a = zeros(Nquad,Nquad,nlayer);
Bm_a = zeros(Nquad,Nquad,nlayer);

epsr_snow = ones(nlayer,1);
ot = zeros(nlayer,1);
albedo = zeros(nlayer,1);

load('Phase_matrix_01.mat');
ka=kappaa; ks = kappas;
load('Phase_matrix_lookup_table.mat')
load('trans12f_01.mat')

for il = 1:nlayer
    %     epsr_ice = diel_ice(fGHz,Tsnow(il));
    %     fv = rho(il)/rho_ice;
    %     fv = max(0,min(1,fv));
    keff_f = interp1(den,keff,rho(il));

    %[keff,ka,ks,bmua,P11,P22,P33,P44,P34,P43] =...
    % PhaseMatrix12f_QMS(k,epsr_ice,dia(il),fv,tau(il),65); % Uncomment
    % this line to run Phase matrix module of DMRT

    P12 = zeros(size(P11));
    P21 = P12;

    %[Fm,Bm,ks2] = tran12f_FB(bmua,P11,P22,P33,P44,P34,P43,P12,P21,xi,wi); % Uncomment
    % this line to run full module of DMRT


    ke = ks + ka;
    epsr_eff = real(keff_f/k).^2;

    epsr_snow(il) = epsr_eff;
    ke_a(il) = ke;
    Fm_a(:,:,il) = Fm;
    Bm_a(:,:,il) = Bm;

    ot(il) = ke*depth(il);
    albedo(il) = ks/ke;

end

eps_mtx = [1;epsr_snow;epsr_ground].';


% calculate bottom interface reflectivity
rv = zeros(ndeg,1);
rh = zeros(ndeg,1);

for ii = 1:ndeg
    if strcmp(rough.model, 'QH')
        [rv(ii),rh(ii)] = rough_ref_QH(epsr_snow(end),epsr_ground,acos(mu0(ii)),rough.Q,rough.H);
    elseif strcmp(rough.model, 'WM')
        [rv(ii),rh(ii)] = rough_ref_WM(epsr_snow(end),epsr_ground,acos(mu0(ii)),k*sqrt(epsr_snow(end))*rough.s);
    else
        error('Unsupported surface reflectivity model');
    end
end

% calculate brightness temperature.
[TBv,TBh,deg0] = fun_DMRTpassive(depth,Tsnow,Tg,eps_mtx,k,ke_a,Fm_a,Bm_a,mu0,wi0,rv,rh);


end
