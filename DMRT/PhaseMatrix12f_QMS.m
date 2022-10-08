function [keff,kappaa,kappas,bmua,P11,P22,P33,P44,P34,P43] = ...
   PhaseMatrix12f_QMS(k,epsr,dia,fv,tau,Nmu)

% epsr = epsr_ice; dia = dia(il); tau = tau(il); Nmu = 65;
%fv = 0.1;
% Calculate effective permittivity, absorption, scattering coefficient and
% phase matrix in 1-2 frame using QCA Mie scattering of Sticky spheres
% (QMS) with Percus-Yevick pair distribution function.
% 
% Inputs:
%   k, wave number in 1/cm
%   epsr, complex dielectric constant of spheres, using exp(iwt) convention
%   dia, diameter of spheres in cm
%   fv, volume fraction of spheres
%   tau, stickiness parameter of the Percus-Yevick pair distribution
%       function.
%       tau > (2 - sqrt(2))/6 = 0.098
%       a larger tau means less stickiness 
%       typical value 0.1.
%   Nmu, number of discretization angle in 1-2 frame. 
%       default value 65.
% Outputs:
%   P11,P22,P33,P44,P34,P43: non-zero elements of the 4 by 4 phase matrix
%       in 1-2 frame, corresponding to the angles TA, use row vector. 
%       Unit is 1/cm.
%   bmua = cos(TA), TA is the sampling angles of phase matrix in 1-2 frames.
%   keff, kappaa, kappas: effective wave number, absorption coefficient,
%       and scattering coefficient in 1/cm.
% 
% Refereces:
%  Tsang, Kong and Ding, Scattering of EM Waves, vol 1.
%  Tsang, Kong, Ding and Ao, Scattering of EM Waves, vol 2.
%  Tsang and Kong, Scattering of EM Waves, vol 3.
% 

% %% Parameters
tol = 1e-4;
   
kp = sqrt(epsr)*k;
tauc=(2-sqrt(2))/6;   % vol. 2, pp.427
vol = 4*pi*(dia/2)^3/3;
rho = fv/vol;

a = dia/2;
ka = k*a;
kpa = kp*a;

Nmax = floor(k*dia) + 1;
Lmax = Nmax*(Nmax + 2);
IPMAX = 2*Nmax;

NS = 1;
nk = 4096;
 % 0.04s
%% Setup Calculation	    

[AA,BB] = COEFF(Nmax);
nzd = Nmax*2;


%ZKF=EFFECTIVE K UNDER FOLDY'S APPROXIMATION! 
[TnM,TnN] = TnMN(Nmax,ka,kpa); 
ztea = -real(TnM)-TnM.*conj(TnM);
ztha = -real(TnN)-TnN.*conj(TnN);

ZsumN = (2*(1:Nmax)+1)*(TnN + TnM);  
ZKF = k - 1i*rho*pi*ZsumN/(k*k);
 
y = (epsr-1)/(epsr+2);
facy = 1i*2*(ka)^3*y/3;
T1N = facy*(1+facy);

% %% It is impoortant to verify the accuracy of the pair function by intgertaion of r^2h(r) from 0 to infinity. 
% It should give H tiltd at p=0, divided by 4 pi.  Thus sumhb should be equal to comp
% This can be achived by increasing NK and RM  so that it goes to a large PMAX

% SHSPYnew1 extracts delta function but not stepfunction
% SHSPYnew extracts both delta function and  stepfunction and should be more accurate
%  
% [dis,grold,dr,HDEL,NR,NKMAX] = SHSPYnew1(tau,rho,fv,dia,vol,nk);
%         integb  = dis.*dis.*(grold - 1);
%         integhb = trapz(dis,integb);
%         sumhb   = integhb + HDEL*dia^2;
%  

[dis,gr,dr,HDEL,NR,NKMAX] = SHSPYnew(tau,rho,fv,dia,nk); % 0.08s

integ   = dis.*dis.*(gr-1);
integh  = trapz(dis,integ);
sumh    = integh + HDEL*dia^2;

pk = 0;
hkft    = hpft(tau,pk,rho,fv,dia,vol);
comp    = hkft/(4*pi);


% %% Muller’s method to find K
       
%page 259 of volume 3

ZY = (kp*kp-k*k)/(kp*kp+2*k*k);                                   
ZD = 1 - fv*ZY;                   
ZK3D = 2*k^3/(3*ZD);                                             
ZK2D = 3*k*k/ZD;                                                                                                       
ZX = fv*ZY*(1+1i*ZK3D*a^3*ZY*(1 + rho*hkft));                
ZKQL = sqrt(k^2+ZK2D*ZX);


LPP0 = -1i/(k*(ZKQL*ZKQL - k*k)) - dia^3/3;
LPP2 = -1i*(ZKQL/k)^2/(k*(ZKQL*ZKQL - k*k));
M0   = dia^3/3 + hkft/(4*pi);
last = 1+2*pi*rho*T1N*1i*(2+(ZKQL/k)^2)/(k*(ZKQL*ZKQL-k*k));
last = last+4*pi*rho*T1N*(-hkft/(4*pi));
ZSUM = MLSUM(k,ZKQL,gr,dis,dr,NS,NKMAX,IPMAX,dia,NR,HDEL);

ZD = SYSEQU(rho,AA,BB,ZSUM,TnM,TnN,NS,Nmax,IPMAX);

ZX(1) = ZKF;                                                 
ZX(2) = real(ZKQL);                                    
ZX(3) = real(ZKF) + 1i*imag(ZKF)/2;



for  II=1:3                                           
    ZKQL  = ZX(II);      
    ZSUM  = MLSUM(k,ZKQL,gr,dis,dr,NS,NKMAX,IPMAX,dia,NR,HDEL);
    ZD    = SYSEQU(rho,AA,BB,ZSUM,TnM,TnN,NS,Nmax,IPMAX);
    deter = det(ZD);
    ZF(II)= det(ZD);
end

RXR = 1;    RXI = 1;
ITR = 0;

while((RXR>tol)||(RXI)>tol)
    ITR = ITR+1;
    ZX = MULLER(ZX,ZF);                                    
    ZF(1) = ZF(2);                                           
    ZF(2) = ZF(3);                                           
    ZKQL  = ZX(3);                                              
    
    ZSUM = MLSUM(k,ZKQL,gr,dis,dr,NS,NKMAX,IPMAX,dia,NR,HDEL);
    ZD   = SYSEQU(rho,AA,BB,ZSUM,TnM,TnN,NS,Nmax,IPMAX);
    
    deter = det(ZD);
    ZF(3) = deter;                                            
    ZDX   = ZX(3)-ZX(2);
    
    RXR = abs(real(ZDX)/real(ZX(3)));                       
    RXI = abs(imag(ZDX)/imag(ZX(3)));
end

keff   = ZX(3);
QCAKR = real(keff);                                         
QCAKI = imag(keff);                                          
PHASE = k/QCAKR;                                            
TLOSS = 2*QCAKI/QCAKR;                 
fack  = k/real(keff);

% %% 16 Distainct Angles (not including negative ones)
phi = 0;

for  nn = 1:Nmax  %1001
    ny = nn + Nmax;
	GN = -pi*1i*rho*(2*nn+1)/k^2;
	   % replacing the last equation by a equation for nonzero right hand side
    ZD(nzd,nn) = GN*TnM(nn);
    ZD(nzd,ny) = GN*TnN(nn);
    RHS        = zeros(nzd,1);
	RHS(nzd,1) = keff - k;
end  %1001
TCOE = ZD\RHS;

reflc = 0;
absor = 0;
for nn = 1:Nmax
    ny = nn+Nmax;
    terma = (2*nn+1)*rho;
    term1 = (-1)^nn*(2*nn + 1);
    reflc = reflc + term1*(-TnM(nn)*TCOE(nn)+TnN(nn)*TCOE(ny))*pi*i*rho/(k*k*(keff + k));
	absor=absor+terma*(abs(TCOE(nn))^2*ztea(nn));
 	absor=absor+terma*(abs(TCOE(ny))^2*ztha(nn));
end
tranm = 1 - reflc;
absor = absor*2*pi/(k*tranm)^2;
absor = absor*k /real(keff);
kappaa = abs(absor);

% %% In 1-2 Frame

rzkq=real(keff);
first = zeros(1,Nmax);
second = zeros(1,Nmax);
for nn = 1:Nmax %1003
    ny = nn + Nmax;
	QN = (2*nn+1)/(nn*(nn+1));
    first(nn)=QN*TnM(nn)*TCOE(nn);
    second(nn)=QN*TnN(nn)*TCOE(ny);
end  %1003

facr = -1i*sqrt(1/(k*rzkq))/(1-reflc);
% Nmu = 65; % 64;
bmua = linspace(-1,1,Nmu);
P11 = zeros(1,Nmu);
P22 = P11;
P33 = P11;
P44 = P11;
P34 = P11;
P43 = P11;

for imu = 1:Nmu  %1004
    bmu  = bmua(imu);
    pk   = sqrt(rzkq^2 + k*k - 2*k*rzkq*bmu);
    hkft = hpft(tau,pk,rho,fv,dia,vol);
	UU   = rho*(1 + rho*hkft);
    
    [pin,taun]=ptn(Nmax,bmu);
    f11=0;  f22=0;
    for n=1:Nmax  %1005
        f11 = f11+first(n)*taun(n)+second(n)*pin(n);
        f22 = f22+first(n)*pin(n)+second(n)*taun(n);
    end  %1005
    f11 = f11*facr;
    f22 = f22*facr;
    f12 = f11*conj(f22);
    
    P11(imu) = abs(f11)^2*UU;
    P22(imu) = abs(f22)^2*UU;
    P33(imu) = real(f12)*UU;
    P44(imu) = P33(imu);
    P34(imu) = -imag(f12)*UU;
    P43(imu) = -P34(imu);
end  %1004

% calculate kappas
Psum   = P11 + P22;
kappas = pi*trapz(bmua,Psum);

% update imaginary part of keff, vol.3 pp. 302.
kappae = kappaa + kappas;
keff = complex(real(keff),kappae/2);

end
% keff
% kappaa
% kappas
% figure; histogram(bmua);
% figure; histogram(P11);
% figure; histogram(P22);
% figure; histogram(P33);
% figure; histogram(P44);
% figure; histogram(P34);
% figure; histogram(P43);
% ===============================================
%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%

function [AA,BB]=COEFF(NMAX)                                
%     THIS SUBROUTINE CALCULATES AA = A(1,N|-1,NU|P)A(N,NU,P) AND       
%     BB = A(1,N|-1,NU|P,P-1)B(N,NU,P). THE RESULTS ARE STORED IN       
%     THE MATRICES AA(N,NU,P) AND BB(N,NU,P)
% However in Fortran P counts from zero, here we add it to 1
% DIMENSION AA(NMAX,NMAX,0:IPMAX),BB(NMAX,NMAX,0:IPMAX)      
IPMAX1=2*NMAX+1;

AA=zeros(NMAX,NMAX,IPMAX1);
BB=zeros(NMAX,NMAX,IPMAX1);
     
for  N=1:NMAX                                                    
  for  NU=1:NMAX                                                   
      RN=N ;                                                  
      RNU=NU ;                                                
      CN=RN*(RN+1);                                                
      CNU=RNU*(RNU+1)  ;                                            
      CDEN=sqrt(CN*CNU);                                             
      IPBEG=abs(N-NU)  ;                                            
      IPEND=N+NU      ;                                              
      for IP=IPBEG:IPEND                                          
          RIP=IP;                                          
          CIP=RIP*(RIP+1);                                        
          W1=Wigner3j(RN,RNU,RIP,1,-1,0);                          
          NT=N+NU+IP;                                               
          ITEST=NT-floor((NT/2))*2 ;                                        
          if (ITEST ==1)                                    
              C1=NT+1;                                        
              C2=RNU+RIP-RN ;                                        
              C3=RN+RIP-RNU ;                                        
              C4=RN+RNU-RIP+1 ;                                    
              CNUM=(2*IP+1)*sqrt(C1*C2*C3*C4) ;                    
              CC=CNUM/CDEN  ;                                        
              W0=Wigner3j(RN,RNU,RIP-1,0,0,0)  ;                    
              BB(N,NU,IP+1)=CC*W1*W0   ;                               
              AA(N,NU,IP+1)=0   ;                                   
          else                                                     
              CNUM=(2*IP+1)*(CN+CNU-CIP)  ;                        
              CC=CNUM/CDEN  ;
              W0=Wigner3j(RN,RNU,RIP,0,0,0) ;  
              AA(N,NU,IP+1)=CC*W1*W0   ;                              
              BB(N,NU,IP+1)=0   ;                                    
          end                                                   
      end                                                    
  end  
end
end

function wigner = Wigner3j(j1,j2,j3,m1,m2,m3)

t1 = j2 - m1 - j3;
t2 = j1 + m2 - j3;
t3 = j1 + j2 - j3;
t4 = j1 - m1;
t5 = j2 + m2;

tmin = max( 0, max( t1, t2 ) );
tmax = min( t3, min( t4, t5 ) );

wigner = 0;

for t = tmin:tmax
    wigner = wigner + (-1)^t / ( factorial(t) * factorial(t-t1) * factorial(t-t2) ...
    * factorial(t3-t) * factorial(t4-t) * factorial(t5-t) );
end

wigner = wigner * (-1)^(j1-j2-m3) ...
    * sqrt( factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) / factorial(j1+j2+j3+1)...
    * factorial(j1+m1) * factorial(j1-m1) * factorial(j2+m2) * factorial(j2-m2) * factorial(j3+m3) * factorial(j3-m3) );
end

% ==============================================
function [TnM,TnN]=TnMN(NMAX,ka,kpa)
% page 98 of volume 1
% CALCULATE T-MATRIX ELEMENTS FROM N=1 TO NMAX

TnM = zeros(NMAX,1);
TnN = zeros(NMAX,1);
Jnp = zeros(NMAX,1);
Jsnp = zeros(NMAX,1);
Hnp = zeros(NMAX,1);


%use matlab built-in functions of bessel and change to nsphericla bessel function
order = 1.5:1:(NMAX+0.5);
Jn=sqrt(pi/(2*ka))*besselj(order,ka);
Jsn=sqrt(pi/(2*kpa))*besselj(order,kpa);
Hn=sqrt(pi/(2*ka))*besselh(order,1,ka);
Jnp(1)=sqrt(pi/(2*ka))*besselj(0.5,ka)-2*Jn(1)/ka;
Jsnp(1)=sqrt(pi/(2*kpa))*besselj(0.5,kpa)-2*Jsn(1)/kpa;
Hnp(1)=sqrt(pi/(2*ka))*besselh(0.5,1,ka)-2*Hn(1)/ka;
if  NMAX>=2
    for II=2:NMAX
        Jnp(II)=Jn(II-1)-(II+1)*Jn(II)/ka;
        Jsnp(II)=Jsn(II-1)-(II+1)*Jsn(II)/kpa;
        Hnp(II)=Hn(II-1)-(II+1)*Hn(II)/ka;
    end
end

for II=1:NMAX
    NUM=Jsn(II)*(Jn(II)+ka*Jnp(II))-Jn(II)*(Jsn(II)+kpa*Jsnp(II));
    DEN=Jsn(II)*(Hn(II)+ka*Hnp(II))-Hn(II)*(Jsn(II)+kpa*Jsnp(II));
    TnM(II)=-NUM/DEN;
    NUM=(kpa)^2*Jsn(II)*(Jn(II)+ka*Jnp(II))-(ka)^2*Jn(II)*(Jsn(II)+kpa*Jsnp(II));
    DEN=(kpa)^2*Jsn(II)*(Hn(II)+ka*Hnp(II))-(ka)^2*Hn(II)*(Jsn(II)+kpa*Jsnp(II));
    TnN(II)=-NUM/DEN;
end

end

% ================================================================
function [DIS,GR,DR,HDEL,NR,NKmax]=SHSPYnew(TAU,RHO,FRAC,DIA,NK)
% new routines all merged in one
%	********************************************************************
%	* This program calculates the pair distribution function of sticky *
%	* hard spheres under the Percus-Yevick approximation               *
%	* DIA     = diameters of NS species                                *
%	* RHO     = number density of spheres                              *
%	* FRAC    = fractional volume                                      *
%	* GR(NK)  = pair distribution function                             *
%	* DIS(NK) = array contains the distancs from the origin            *
%  RM is the maximum distance for the pair function
RM = 64;

%note that DK*R<<1  for fourier transform or sine transform  intgeral to be valid
% also NK*DK should be big enough so that the regular intgerand has decayed
XI = FRAC;
XI1 = 1-XI;
RNU = TAU+XI/XI1;
RGAMA = XI*(1+XI/2)/3/XI1/XI1;
RLAM = 6*(RNU-sqrt(RNU*RNU-RGAMA))/XI;
RMU = RLAM*XI*XI1;
if (RMU>(1+2*XI)) return; end
F = zeros(NK,1);
FREG = zeros(NK,1);
CST1 = XI/(1-XI);
CST2 = 1-RLAM*XI+3*CST1;
CST3 = 3-RLAM*(1-XI);
AS = CST1*CST2;
BS = CST1*CST3;
DS = CST1;

T1 = -BS+3*DS;
T2 = 6*AS-9*DS*DS;
T3 = -(BS*BS+6*DS);
T4 = (-2*BS+6*DS)^2;
DK = pi/(RM*DIA);
for IK=2:NK
    PK=(IK-1)*DK;
    X=PK*DIA/2;
    SNX=sin(X);
    CSX=cos(X);
    PSIX=SNX/X;
    PHIX=3*(SNX-X*CSX)/X/X/X;
    PP1=CST1*(CST2*PHIX+CST3*PSIX)+CSX;
    PP2=CST1*X*PHIX+SNX;
    CFT=PP1*PP1+PP2*PP2;
    CFTINV=1/CFT;
    F(IK)=X*(CFTINV-1);
    FASYI=T2*CSX*CSX/X;
    FASYI=FASYI+T3*SNX*SNX/X;
    FASYI=FASYI+SNX*SNX*CSX*CSX*T4/X;
    FASY=2*SNX*CSX*T1+FASYI;
    FREG(IK)=F(IK)-FASY;
    XA(IK)=X;
end
  %  at  IR=1+ NK/(2*RM), R is eqaul to DIA
DR = 2*RM*DIA/NK;
NR = 1+ NK/(2*RM);
%note also DK=pi/(RM*DIA),  results not accurate for large R  DK*maxR=0.4

maxR = 16*DIA;
NKmax = floor(maxR/DR)+1;
DIS=zeros(NKmax,1);
    GR=zeros(NKmax,1);
    termr = (2/DIA)^2/(2*pi*pi*RHO);
for IR=1:NKmax
    R=(IR-1)*DR;
    DIS(IR)=R;
    if (IR>=NR)
        Y=2*R/DIA;
        term=termr/R;
        %intgeration by Simpson's rule in fortrab, here trappezoidal 's rule
        for IK=1:NK
            if  (IK==1)
                XA(IK)=0;
                integ(IK)=-T2*Y;  
            else
                integ(IK)=sin(XA(IK)*Y)*FREG(IK);
            end
        end
        res=trapz(XA,integ);
        HRreg=res*term;
        %integrals
        [integ1,integ2,integ3]=integr(Y);
        HRasy=(T2*integ1+T3*integ2+T4*integ3);
        HR=HRasy*term+HRreg;
        GR(IR)=1+HR;
    else
        GR(IR)=0;
    end
end
HDEL=pi*termr*T1/4;
end
    
function [integ2,integ3,integ4]=integr(Y)
integ2=istandard(Y)/2+istandard(Y+2)/4+istandard(Y-2)/4;
integ3=istandard(Y)/2-istandard(Y+2)/4-istandard(Y-2)/4;
integ4=istandard(Y)/8-istandard(Y+4)/16-istandard(Y-4)/16;
end

function in=istandard(Y)
if  Y>0
    in=pi/2;
elseif Y==0
    in=0;
else % if  (Y<0)
    in=-pi/2;
end
end

% ==================================================
function HKFT=hpft(TAU,PK,RHO,FRAC,DIA,VOL)
%page 427 of volume 2
%first solve for t, RLAM is the t, RMU is (8.4.18) on page 426
% pk is p
XI=FRAC;
XI1=1-XI;
RNU=TAU+XI/XI1;
RGAMA=XI*(1+XI/2)/(3*XI1*XI1);
RLAM=6*(RNU-sqrt(RNU*RNU-RGAMA))/XI;
RMU=RLAM*XI*XI1;
if  (RMU>(1+2*XI))
    return
end
% next do (8.4.19) on page 426
CST1=XI/(1-XI);
CST2=1-RLAM*XI+3*CST1;
CST3=3-RLAM*(1-XI);
X=PK*DIA/2;
if  (PK==0) 
    PP1=CST1*(4-RLAM+3*CST1)+1;
    PP2=0;
else
    SNX=sin(X);
    CSX=cos(X);
    PSIX=SNX/X;
    PHIX=3*(SNX-X*CSX)/X/X/X;
    PP1=CST1*(CST2*PHIX+CST3*PSIX)+CSX;
    PP2=CST1*X*PHIX+SNX;
end
    CFT=PP1*PP1+PP2*PP2;
    CFTINV=1/CFT;
    HKFT=(CFTINV-1)/RHO;
    %HKFT is H tilde at p
end

% ========================================
function [pin,taun]=ptn(Nmax,x)
% pi tau functions
% page 543 of volume 2

pin = zeros(Nmax,1);
taun = zeros(Nmax,1);
if  (Nmax==1)
    pin(1)=1;
    taun(1)=x;
end

if  (Nmax==2)
    pin(1) = 1;
    taun(1) = x;
    n = 2;
    pin(n) = 3*x;
    taun(n) = n*x*pin(n)-(n+1)*pin(n-1);
end

if (Nmax>=3)
   pin(1)=1;
   taun(1)=x;
   n = 2;
   pin(n)=3*x;
   taun(n)=n*x*pin(n)-(n+1)*pin(n-1);
   for n=3:Nmax
        pin(n)=((2*n-1)*x*pin(n-1)-n*pin(n-2))/(n-1);
        taun(n)=n*x*pin(n)-(n+1)*pin(n-1);
   end
end
end

% =========================================
function  ZSUM=MLSUM(BK,ZK,GR,DIS,DR,NS,NK,IPMAX,DIA,NR,HDEL)
%C	ZSUM= volume 3  Lp+Mp on pgae 253 and page 254
for IR=NR:NK 
    ind=(IR-NR)+1;                                         
    R=DIS(IR);                                                    
    ZBKR=BK*R   ;                                             
    ZKR=ZK*R  ;
    %use matlab built-in functions of bessel and change to nsphericla bessel function
    id=0;
    [ZSBJ,bp]=sbessel(IPMAX,ZKR,id);
    [ZSHH,hp]=shankel(IPMAX,ZBKR,id);
    for KP=1:(IPMAX +1)                                             
        ZPHJ(KP,ind)=R*R*ZSBJ(KP)*ZSHH(KP)   ;                      
    end
end                                                      
  %intgertaion of (6.1.39) of page 254 of volume 3 using simpson'srule 
  %here we use trapezoidal rule
for IP=0:IPMAX                                                                                                                                   
    for IK=NR:NK 
        ind=(IK-NR)+1;
        RA(ind)=DIS(IK);
        integ(ind)=(GR(IK)-1)*ZPHJ(IP+1,ind) ;
     end
     ZRES(IP+1)=trapz(RA,integ);
 end
% note from page 257 of volume 3 , real part of ZSUMF is causing attnuation while imaginary part 
% is contributing to chnaging the realpart of the effective propgation constant
% note that the imaginary part of ZSUMF is dominated by lpp
    R=DIA;
    ZBKR=BK*R;                                              
	ZKR=ZK*R;
    id=1;
    [zsbj,jpp]=sbessel(IPMAX,ZKR,id);
    [zshh,hpp]=shankel(IPMAX,ZBKR,id);
for ip=0:IPMAX
    ipa=ip+1;
    lpp=-R*R/(ZK*ZK-BK*BK);                                        
	lpp=lpp*(BK*hpp(ipa)*zsbj(ipa)-ZK*zshh(ipa)*jpp(ipa));
    delpart=HDEL*R*R*zshh(ipa)*zsbj(ipa);
    Zint=ZRES(ipa);
    Mpart=delpart+Zint ;
    ZSUMold=Zint+lpp;
    ZSUMF=ZSUMold+delpart;
    ZSUM(ipa)=ZSUMF;
end
end

function [b,bp]=sbessel(np,arg,id)
%calcuate spheeical bessel functions from order 0 to order np
%spherical bessel functions for matlab by using matlab built in functions
%if id=0 only want functions
% if id=1, want both functions and derivatives
number=np+1;
order=0:np;
orderB=order+0.5;
bp=zeros(1,number);
% b=zeros(1,number);
b=sqrt(pi/(2*arg))*besselj(orderB,arg);
if  (id==1)
   for in=1:number
       ip=order(in);
       if  (ip==0)
          bp(in)=-b(in+1);
       else
            bp(in)=b(in-1)-(ip+1)*b(in)/arg;
       end
    end
end
end

function [b,bp]=shankel(np,arg,id)
%calcuate spherical hankel functions from order 0 to order np
%spherical bessel functions for matlab by using matlab built in functions
%if id=0 only want functions
% if id=1, want both functions and derivatives
number=np+1;
order=0:np;
orderB=order+0.5;
bp=zeros(1,number);
% b=zeros(1,number);
b=sqrt(pi/(2*arg))*besselh(orderB,1,arg);
if  (id==1)
   for in=1:number
       ip=order(in);
       if  (ip==0)
         bp(in)=-b(in+1);
       else
        bp(in)=b(in-1)-(ip+1)*b(in)/arg;
       end
   end
end
end

%  ========================================
function ZD=SYSEQU(RHO,AA,BB,ZSUM,ZTE,ZTH,NS,NMAX,IPMAX)
nzd=NMAX*2;
ZD=zeros(nzd);
for NU=1:NMAX
  NUX=NU    ;                                                
  NUY=NUX+NMAX ;                                                   
  for  N=1:NMAX
      NX=N ;                                               
      NY=NX+NMAX  ;                                                
      IPBEG=abs(N-NU) ;                                        
      IPEND=N+NU    ;                                            
      for IP=IPBEG:IPEND   
          ipa=IP+1;
          ZD(NUX,NX)=ZD(NUX,NX)+AA(N,NU,ipa)*ZSUM(ipa)  ;    
          ZD(NUY,NX)=ZD(NUY,NX)+BB(N,NU,ipa)*ZSUM(ipa) ;     
          ZD(NUX,NY)=ZD(NUX,NY)+BB(N,NU,ipa)*ZSUM(ipa) ;     
          ZD(NUY,NY)=ZD(NUY,NY)+AA(N,NU,ipa)*ZSUM(ipa) ;     
      end     
      FN=2*pi*RHO*(2*N+1)   ;
      ZD(NUX,NX)=ZD(NUX,NX)*FN*ZTE(N);                        
      ZD(NUY,NX)=ZD(NUY,NX)*FN*ZTE(N);                        
      ZD(NUX,NY)=ZD(NUX,NY)*FN*ZTH(N);                        
      ZD(NUY,NY)=ZD(NUY,NY)*FN*ZTH(N) ;                       
  end
end
for  NU=1:nzd                                                  
  ZD(NU,NU)=1+ZD(NU,NU) ;                                       
end
end

% ===================================
 function ZY=MULLER(ZX,ZF)                                          
%(A)	INPUT
%	ZX=THREE ROOTS
%	ZF=ARE THE THREE CORRESPOMDING DETERMINANTAL VALUES
%(B) OUTPUT
%	ZY(1)=ZX(2)
%	ZY(2)=ZX(3)
%	ZY(3)=NEW ROOT FROM MULLER'S METHOD
%    MAKE SURE SOLUTION 3 IS BETTER THAN SOLUTION 2                    
  ZRTO=ZF(3)/ZF(2);                                                  
  if (abs(ZRTO)> 1)                                       
      ZFTMP=ZF(2);                                                   
      ZF(2)=ZF(3);                                                   
      ZF(3)=ZFTMP;                                                   
      ZXTMP=ZX(2);                                                   
      ZX(2)=ZX(3);                                                   
      ZX(3)=ZXTMP;                                                   
  end                                                        
  ZQ=(ZX(3)-ZX(2))/(ZX(2)-ZX(1));                                    
  ZP=1+ZQ;                                                          
  ZA=ZQ*ZF(3)-ZQ*ZP*ZF(2)+ZQ*ZQ*ZF(1);                               
  ZB=(2*ZQ+1)*ZF(3)-ZP*ZP*ZF(2)+ZQ*ZQ*ZF(1);                       
  ZC=ZP*ZF(3);                                                       
  ZD=sqrt(ZB*ZB-4.*ZA*ZC);                                          
  U1=abs(ZB+ZD);                                                    
  U2=abs(ZB-ZD);   
  if (U1 > U2)                                             
      ZU=ZB+ZD;                                                      
  else                                                           
      ZU=ZB-ZD;                                                      
  end   
%      Z=ZX(3)-(ZX(3)-ZX(2))*ZC/ZU ,2/24/99                             
  Z=ZX(3)-2*(ZX(3)-ZX(2))*ZC/ZU;                                   
  ZY(1)=ZX(2);                                                       
  ZY(2)=ZX(3);                                                       
  ZY(3)=Z;  
 end
 
%  =========================================================
