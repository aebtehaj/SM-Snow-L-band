function [Tb,gammav,gammah,espr_snow,espr_ground,Tb_DMRT,Tb_upwelling, Tb_downwelling_reflected] =tau_omega4(Tcanopy,tau,omega,h,freq,obs_angle,Tsnow,Tg,mv,clayfrac,rho_s)

% Soil roughness model
rough.model = 'QH';
rough.Q = 0;          
rough.H = h;

% Vegetation Microwave Transmissivity
angler = obs_angle.*pi./180 ;
exptau= exp(-tau./cos(angler));

x = [30,rho_s./1000,0.1];
n = 1;
tau_DMRT = 0.1;
[Tb_DMRT,espr_snow,espr_ground]= DMRT_QMS_PM(freq,obs_angle,Tsnow,Tg,mv,clayfrac,n,rough,tau_DMRT,x);
d0 = linspace(0.1,2,1000);
%gammav_m = zeros(1000); gammah_m = zeros(1000);
for i = 1:length(d0)
    [~, ~, gammav_m(i), gammah_m(i)] = Refl_TwoLayerComposite3(espr_snow, espr_ground, d0(i), 40, freq,h);
end
gammav = mean(gammav_m);
gammah = mean(gammah_m);

Tb_upwelling = Tcanopy.*(1-omega).*(1-exptau);
Tb_downwelling_reflected = [Tcanopy.*(1-omega).*(1-exptau).*gammav.*exptau;Tcanopy.*(1-omega).*(1-exptau).*gammah.*exptau];

TBv  = Tb_DMRT(1,1).*exptau + Tcanopy.*(1-omega).*(1-exptau).*gammav.*exptau + Tcanopy.*(1-omega).*(1-exptau);
TBh  = Tb_DMRT(2,1).*exptau + Tcanopy.*(1-omega).*(1-exptau).*gammah.*exptau + Tcanopy.*(1-omega).*(1-exptau);

Tb=[TBv;TBh];
end