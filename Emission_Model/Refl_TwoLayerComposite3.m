%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. 
%%
%Code 2.4: Oblique Reflection @ Planar Boundries for a Two-Layer Composite
%Description: Code computes the reflection coefficients, transmission
    %coefficients, reflectivities and transmissivities for incidence in
    %lossless medium (medium 1) upon the planar boundary of a two layer
    %composite that may be lossless or lossy (middle layer is medium 2 and 
    %bottom layer is medium 3) at any incidence angle, for both h and v
    %polarizations
    
%Input Variables:
    %eps1: relative dielectric constant of medium 1, set here to 1
    %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
    %eps3 = eps3r-j*eps3i: relative dielectric constant of medium 3
    %d: vertical thickness of middle layer in meters
    %theta1: incidence angle in medium 1 in degrees
    %f: frequency in GHz
    
%Output Products:
    %rhoh: reflection coefficient for h pol
    %rhov: reflection coefficient for v pol
    %gammah:reflectivity for h pol
    %gammav: reflectivity for v pol
    
%Book Reference: Section 2-10

%Example call: [rhoh rhov gammah gammav] = Refl_TwoLayerComposite(eps2, eps3, d, theta1, f)

%MATLAB Code

function [rhoh, rhov, gammav, gammah] = Refl_TwoLayerComposite3(eps2, eps3, d, theta1, f,h)

theta1 = deg2rad(theta1);    
eps1 = 1;
    
    lam2 = 1i *20*pi/3*f*sqrt(eps2); % wave number in medium 2
    
    theta2 = acos((1-(sqrt(eps1)/sqrt(eps2).*sin(theta1)).^2).^(1/2));
    theta3 = acos((1-(sqrt(eps1)/sqrt(eps3).*sin(theta1)).^2).^(1/2));
    
    rho12h = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos(theta2)) ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos(theta2));
    rho12v = (sqrt(eps1).*cos(theta2)-sqrt(eps2).*cos(theta1)) ./ (sqrt(eps1).*cos(theta2) + sqrt(eps2)*cos(theta1));
    
    rho23h = (sqrt(eps2).*cos(theta2)-sqrt(eps3).*cos(theta3)) ./ (sqrt(eps2).*cos(theta2) + sqrt(eps3).*cos(theta3));
    rho23v = (sqrt(eps2).*cos(theta3)-sqrt(eps3).*cos(theta2)) ./ (sqrt(eps2).*cos(theta3) + sqrt(eps3)*cos(theta2));
    
    rhoh = (rho12h + rho23h.*exp(-2*lam2.*d.*cos(theta2)-0.5.*h.*cos(theta2).*cos(theta2)))./(1 + rho12h.*rho23h.*exp(-2*lam2.*d.*cos(theta2)-0.5.*h.*cos(theta2).*cos(theta2)));
    rhov = (rho12v + rho23v.*exp(-2*lam2.*d.*cos(theta2)-0.5.*h.*cos(theta2).*cos(theta2)))./(1 + rho12v.*rho23v.*exp(-2*lam2.*d.*cos(theta2)-0.5.*h.*cos(theta2).*cos(theta2)));

    gammah = abs(rhoh).^2;
    gammav = abs(rhov).^2;
    
       
  end