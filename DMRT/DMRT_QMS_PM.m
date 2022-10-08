% Implementation of DMRT_QMS Module    
% Input parameters:
    %   freq - frequncy in GHz, vector
    %   rho_s  - snow density in gm/cc, 1D array
    %   Tsnow - snow temperature in Kelvin, 1D array
    %   Tg    - ground temperature in Kelvin, scalar
    
% Output:
    %   TBv, TBh - Brightness Temperature in vertical and horizontal
    %       polrization
    %
    % Ref:  Liang et al., TGRS, 46(11): 3663-3671, 2008
    % Copyright: University of Washington, Electrical Engineering Department
    % Revision Date: Aug. 1, 2014, Sep. 09, 2014.
    % DMRT_QMS_passive(fGHz,dia,rho,tau,depth,Tsnow,Tg,epsr_ground,rough)

function [Tb,epsr_snow,epsr_ground] = DMRT_QMS_PM(freq,obs_angle,Tsnow,Tg,mv,clayfrac,rho_s,h)

% Soil roughness model
rough.model = 'QH';
rough.Q = 0;          
rough.H = h;

n = 1;
tau_DMRT = 0.1;
depth = 30;
rho = rho_s./1000;
dia = 0.1;

tau = tau_DMRT;
tau = tau*ones(1,n);
for i = 1:length(freq)
    if Tg>273.15
     epsr_ground = soil_perm_MBSDM_Mironov(mv,clayfrac,freq(i));
    elseif Tg<273.15
      epsr_ground = 5 + 0.5i; % For frozen ground
    end
end

    [TBv,TBh,deg0,~,~,epsr_snow] = DMRT_QMS_passive(freq(i),dia,rho,tau,depth,Tsnow,Tg,epsr_ground,rough);
  
    Tb_v = spline(deg0,TBv,obs_angle(i));
    Tb_h = spline(deg0,TBh,obs_angle(i));
    Tb = [Tb_v;Tb_h];
    
end