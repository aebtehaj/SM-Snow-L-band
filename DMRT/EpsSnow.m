%   Snow Layer Radiative Transfer Properties  =================================

%  (The code is a direct adaptation of Mike Shwanke's Mathematica code
%  collected in 2S-EM Model.)
%  Divya Kumawat, 06/2022
%--------------------------------------------------

% Inputs =============================================
% RoS: snow mass-density [kg m^-3] 
% WS: volumetric snow liquid water content [-] 
% fGHz: frequency [GHz]
% Outputs ============================================
% epsS: permittivity of snow [-]. [Snow can be dry (WS=0, epsS=real) or wet (WS>0, epsS=complex)]

function epsS  =  EpsSnow(RoS,WS,fGHz)

roS=RoS/1000.0;
vfj=roS/0.917;
ehb=0.99913;
esb=1.4759;

if roS<=0.4
    epsSdry=1+1.5995*roS+1.861*roS^3;

elseif  roS>0.4
        epsSdry=((1-vfj)*ehb+vfj*esb)^3;
end

if  WS>0

    % permittivity 'ew' of water at TW = 273.15 K*)
    TW = 273.15;
    TETA=1-300/TW;
    e0=77.66-103.3*TETA;
    e1=0.0671*e0;
    f1=20.2+146.4*TETA+316*TETA*TETA;
    e2=3.52+7.52*TETA;
    f2=39.8*f1;
    ew=e2+(e1-e2)/(1-1i*fGHz/f2)+(e0-e1)/(1-1i*fGHz/f1);
    % permittivity 'epsSwet' of wet snow
    Aa=0.005; % depolarisation factors of prolate*)
    Ab=0.4975; % water inclusion (Matzler 1987)*)
    Ac=Ab;
    Ka=epsSdry/(epsSdry+Aa*(ew-epsSdry));
    Kb=epsSdry/(epsSdry+Ab*(ew-epsSdry));
    K=(Ka+2*Kb)/3;
    epsz=(1-WS)*epsSdry+WS*ew*K;
    epsn=1-WS*(1-K);
    epsSwet=epsz/epsn;

end


if  WS==0
    epsS=epsSdry;
elseif WS>0
        epsS=epsSwet;

end
end


