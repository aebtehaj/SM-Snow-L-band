%Relative Dielectric Constant of Dry Snow
%Description: Code computes the real and imaginary parts of the relative
%dielectric constant of Dry Snow
%Input Variables:
    %T: Temperature in C
    %Ps: Dry Snow Density in g/cm^3
    %f: frequency in GHz
%Output Products:
    %epsr: real part of relative dielectric constant

%Example call: [A B] = RelDielConst_DrySnow(T,Ps,f)
%Computes the real and imaginary components of the permitivity of Dry Snow
    %based on the temperature value (T) in degrees C, dry snow density (Ps), and frequency
    %vector (f) and assigns them to vectors A and B respectively

%MATAB CODE

function [epsr_ds epsi_ds] = RelDielConst_DrySnow(T, Ps, f )

vi=  Ps /0.9167 ;


[epsr_ice, epsi_ice] = RelDielConst_PureIce(T,f);

if vi <= 0.45
    epsr_ds = 1 + 1.4667 .* vi + 1.435 .* vi.^3; % 0<vi<0.45
end
if vi > 0.45
    epsr_ds = (1+ 0.4759 .*vi).^3; % vi>0.45
end

epsi_ds = 0.34 * vi * epsi_ice ./(1- 0.42 * vi).^2;


end

