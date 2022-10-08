function epsr_ice = diel_ice(fGHz,T)
% dry ice dielectric emperical model by Matzler and Wegmuller
% real part valid for 243 < T < 273 K
% Reference: Matzler 2006.
% 

epsr_re = 3.1884 + 0.00091*(T-273.0); % 243 < T < 273 K

v = fGHz;
T0 = 300;
Ta = T0/T - 1;
alpha = (0.00504 + 0.0062*Ta)*exp(-22.1*Ta);
B1 = 0.0207; %K GHz-1
b = 335; % K
B2 = 1.16e-11; %GHz-3
betaM = B1/T*exp(b/T)/(exp(b/T) - 1)^2 + B2*v^2;
dbeta = exp(-9.963 + 0.0372*(T - 273.16));
beta = betaM + dbeta;

epsr_im = alpha/v + beta*v;
epsr_ice = epsr_re + 1i*epsr_im;
end