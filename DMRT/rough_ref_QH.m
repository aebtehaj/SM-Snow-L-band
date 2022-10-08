function [rv,rh] = rough_ref_QH(eps1,eps2,tai,Q,H)
% calculate reflectivity of rough interface using Q/H model
% Q depolarization factor, dimensionless
% H dimensionless parameter characterizing roughness height
% Q = 0, H = 0 degrades to flat interface
% an estimation of H, H = (2*k*s)^2, following KA, s is rms height
% 	k is wave number in medium 1
% Reference:
%   Wang and Choudhury, J. Geophysical Research, 86(C6): 5277-5282, 1981
% 

[rv0,rh0] = Fresnel(eps1,eps2,tai);
G = cos(tai)^2;
att = exp(- H*G);
rv = ((1 - Q)*rv0 + Q*rh0)*att;
rh = ((1 - Q)*rh0 + Q*rv0)*att;

end