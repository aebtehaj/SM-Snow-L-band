function [rv,rh,Rv,Rh] = Fresnel(eps1,eps2,tai)
% Calculate reflectivity of smooth interface between two medium eps1 and
% eps2, the incidence angle in medium 1 is tai
% In RT, assume eps1 to be real, and eps2 could be complex
% 

eb = eps2/eps1;
ub = 1;
Rv = (eb*cos(tai) - sqrt(ub*eb - sin(tai)^2))/(eb*cos(tai) + sqrt(ub*eb - sin(tai)^2));
Rh = (ub*cos(tai) - sqrt(ub*eb - sin(tai)^2))/(ub*cos(tai) + sqrt(ub*eb - sin(tai)^2));

rv = abs(Rv)^2;
rh = abs(Rh)^2;

end