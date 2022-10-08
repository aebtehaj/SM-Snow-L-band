function [TBv,TBh,deg0] = fun_DMRTpassive(t,Ts,Tg,eps_mtx,k,ke,Fm,Bm,mu0,wi0,rv,rh)
% solve passive DMRT equation using Discrete Ordinate-Eigenanalysis 
% ( Guass-Legendre quadrature ) method, support multilayer snowpack
% configuration.
% Inputs:
%   t, thickness of each layer in cm, length N
%   Ts, snow temperature of each layer in Kelvin, length N
%   Tg, ground temperature in Kelvin, length N
%   eps_mtx, effective permittivity array, length N + 2
%   k, free space wave number in 1/cm
%   ke, extinction coefficient in each layer, 1/cm
%   Fm, Bm, Forward and Backward phase matrix (2x2) in body frame.
%   mu0, cosines of quadrature angles, 0 -> pi/2, half, column vector
%   wi0, corresponding quadrature weights
%   rv, rh, reflectivity of bottom snow-soil boundary at the quadrature
%       angles.
% Outputs:
%   deg0, quadrature angles in air
%   TBv, TBh, brightness temperatures on deg0
% 
% Ref: Liang et al., TGRS, 46(11): 3663-3671, 2008 
% 


depth = cumsum([0;t(:)]); % position of boundaries

% global wi xi Nquad ndeg 
% mu0 = -xi(1:ndeg);  wi0 = wi(1:ndeg);

ndeg = length(mu0);
Nquad = 2*ndeg;

thetap = acos(mu0);  % thetam = acos(-mu0);
% thetapd = thetap*180/pi;
%=====================================================
% Calculate Reflectivity & Emissivity
%=====================================================
kmtx = k*sqrt(eps_mtx);
nlayer = size(kmtx,2)-2;
%     Ts = Tsnow;
%     depth = depth_guess;

% optdep = ke*depth_guess;
r = zeros(2*ndeg,2*nlayer);
for ii = 1:nlayer
    %kr = kmtx(ii+1)*sin(thetap);
    r(:,2*ii-1) = rflct(k,thetap,eps_mtx(ii+1),eps_mtx(ii));
%    mu0 = -xi(1:ndeg); 

%     r(:,2*ii)   = rflct(k,kr,eps_mtx(ii+1),eps_mtx(ii+2));
    r(:,2*ii)   = rflct(k,thetap,eps_mtx(ii+1),eps_mtx(ii+2));
%     rQH = rflctQH(Q,H,k,thetap,eps_mtx(ii+1),eps_mtx(ii+2));
%     r(:,2*ii)  = rQH;

    % ke_d(ii) = ke(ii)*(depth(ii+1) - depth(ii));
end
% modify the bottom reflectivity with Q/H model.
% r(:,2*ii) = rflctQH(Q,H,k,thetap,eps_mtx(ii+1),eps_mtx(ii+2));
r(:,2*ii) = [rv(:);rh(:)];

t = ones(size(r)) - r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=====================================================
% RTE using Discrete Ordinate Method
%=====================================================
awi = diag([wi0;wi0]);
%miu = diag([mu0;mu0]);  miu_1 = inv(miu);
miu_1 = diag(1./[mu0;mu0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculte eigenvectors and eigenvalue
%=====================================================
% Load Phase Matrix
%=====================================================
% Matrix = zeros(nlayer*2*Nquad);
Matrix = sparse(nlayer*2*Nquad,nlayer*2*Nquad);
T = zeros(nlayer*2*Nquad,1);
jj = 0;
for ii = 1:2*nlayer  
  if jj ~= round(ii/2)
    jj = round(ii/2);
    ked = ke(jj)*eye(32);
    W = -ked + Fm(:,:,jj)*awi - Bm(:,:,jj)*awi;
    A = -ked + Fm(:,:,jj)*awi + Bm(:,:,jj)*awi;
    X = miu_1*W*miu_1*A;
    [V,D] = eig(X);
    alpha = sqrt(D*ones(Nquad,1));
  end
  
  if ii==1
     A1 = A; V1 = V; alpha1 = alpha;
  end
%--------------------------------------------------------------------------
% Match Boundry Condition
%----------diag----------------
if mod(ii,2)
[Nu Nd] = Nmatrix(0,depth(jj+1)-depth(jj),miu_1,Nquad,alpha,A,V); 
Nz1 = Nd-diag(r(:,ii))*Nu;
% depth(jj+1)-depth(jj);
else
[Nu Nd] = Nmatrix(depth(jj)-depth(jj+1),0,miu_1,Nquad,alpha,A,V);   
Nz2 = Nu-diag(r(:,ii))*Nd;
% depth(jj)-depth(jj+1);
index = (jj-1)*2*Nquad+1:jj*2*Nquad;
Matrix(index,index) = [Nz1;Nz2];
end

%----------others---------------
if ((ii+1)/jj == 2) && (ii ~= 1)
    Tr1 = interpolation(Nquad,eps_mtx(jj),eps_mtx(jj+1),mu0);
    index11 = (ii-2)*Nquad+1:(ii-1)*Nquad;
    index12 = (ii-1)*Nquad+1:(ii+1)*Nquad;
    Matrix(index11,index12) = -diag(t(:,ii))*Tr1*Nu;
    Tr2 = interpolation(Nquad,eps_mtx(jj+1),eps_mtx(jj),mu0);
    index21 = (ii-1)*Nquad+1:ii*Nquad;
    index22 = (ii-3)*Nquad+1:(ii-1)*Nquad;
    Matrix(index21,index22) = -diag(t(:,ii-1))*Tr2*Ndpast;
    T((ii-2)*Nquad+1:(ii-1)*Nquad) = -2*Ts(jj)*(t(:,ii-1) - Tr1*t(:,ii)); 
    T((ii-1)*Nquad+1:ii*Nquad) =  -2*Ts(jj)*(t(:,ii) - Tr2*t(:,ii-1));  
end
Ndpast = Nd;
end
T(1:Nquad)=-2*Ts(1)*t(:,1);
T(end-Nquad+1:end) = 2*t(:,end)*(Tg-Ts(nlayer));

 P = Matrix\T;
 P = P(1:2*Nquad);

IuH = 0;
for ii = 1:Nquad
   IuH = IuH + 1/2*((eye(Nquad,Nquad) + 1/alpha1(ii)*miu_1*A1)*V1(:,ii)*P(ii)...
       + (eye(Nquad,Nquad) - 1/alpha1(ii)*miu_1*A1)*V1(:,ii).*exp(-alpha1(ii)*depth(2))*P(Nquad+ii));
end

Iu = IuH + Ts(1);
TB = diag(t(:,1))*Iu;
TBv = TB(1:Nquad/2,1);
TBh = TB(Nquad/2+1:Nquad,1);

% Observation
ta_crt = asin(k/kmtx(2)); % Brewster angle in snow
ta = acos(mu0);           % quadrature angle in snow  
ta_id = ta < ta_crt;      % find thoes angles smaller than Brewster angle  
deg0 = asin(kmtx(2)*sin(ta(ta_id))/k)/pi*180; % quadrature angles in air
TBv = TBv(ta_id);
TBh = TBh(ta_id);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Observation
% ob_anglei = linspace(10,90,90);
% for id = 1:90
% ob_angle = ob_anglei(id);
%    theta_obs = real(asin(kmtx(2)/k.*sin(acos(mu0)))*180/pi);
%    lb = max(find(theta_obs<ob_angle));		% Define lowerbound and upperbound
%    ub = lb+1;
%    % linear interpolation
%    Tb_v(id) = TBv(lb)+(ob_angle-theta_obs(lb))/(theta_obs(ub)-theta_obs(lb))*(TBv(ub)-TBv(lb));
%    Tb_h(id) = TBh(lb)+(ob_angle-theta_obs(lb))/(theta_obs(ub)-theta_obs(lb))*(TBh(ub)-TBh(lb));
% end
% ob_angle = 54; 
%    theta_obs = real(asin(kmtx(2)/k.*sin(acos(mu0)))*180/pi);
%    lb = find(theta_obs<ob_angle, 1, 'last' );		% Define lowerbound and upperbound
%    ub = lb+1;
%    % linear interpolatio
%    TBv = TBv(lb)+(ob_angle-theta_obs(lb))/(theta_obs(ub)-theta_obs(lb))*(TBv(ub)-TBv(lb));
%    TBh = TBh(lb)+(ob_angle-theta_obs(lb))/(theta_obs(ub)-theta_obs(lb))*(TBh(ub)-TBh(lb));
end

function r12 = rflct(k,thetap,e1,e2)
k1 = sqrt(e1)*k;    
k2 = sqrt(e2)*k;

k1z = sqrt(k1.^2 - (k1.*sin(thetap)).^2);
k2z = sqrt(k2.^2 - (k1.*sin(thetap)).^2);

if (imag(k1z)<0)
    k1z=-k1z;
end
if (imag(k2z)<0)
    k2z=-k2z;
end

%-----------------------------------------------------------
r12h = abs((k2z-k1z)./(k2z+k1z)).^2;
r12v = abs((e1*k2z-e2*k1z)./(e1*k2z+e2*k1z)).^2;
r12 = [r12v;r12h];
end

% function r12 = rflctQH(Q,H,k,thetap,e1,e2)
% % coherent reflectivity with roughness using Q/H model
% % Ref: Wang et al., TGRS, GE-21(1):44-51, 1983
% 
% k1 = sqrt(e1)*k;    
% k2 = sqrt(e2)*k;
% 
% k1z = sqrt(k1.^2 - (k1.*sin(thetap)).^2);
% k2z = sqrt(k2.^2 - (k1.*sin(thetap)).^2);
% 
% if (imag(k1z)<0)
%     k1z=-k1z;
% end
% if (imag(k2z)<0)
%     k2z=-k2z;
% end
% 
% %-----------------------------------------------------------
% r12h = abs((k2z-k1z)./(k2z+k1z)).^2;
% r12v = abs((e1*k2z-e2*k1z)./(e1*k2z+e2*k1z)).^2;
% r12v = ((1-Q)*r12v + Q*r12h).*exp(-H*cos(thetap).^2);
% r12h = ((1-Q)*r12h + Q*r12v).*exp(-H*cos(thetap).^2);
% r12 = [r12v;r12h];
% end

function [Nu,Nd] = Nmatrix(d1,d2,miu_1,n,alpha,A,V)

Nuplus = zeros(n,n);Ndplus = zeros(n,n);Numinus = zeros(n,n);Ndminus = zeros(n,n);
for i = 1:n
    Nuplus(:,i)  = (eye(n,n) + 1/alpha(i)*miu_1*A)*V(:,i).*exp(alpha(i)*d1); 
    Numinus(:,i) = (eye(n,n) - 1/alpha(i)*miu_1*A)*V(:,i).*exp(-alpha(i)*d2); 
    Ndplus(:,i)  = (eye(n,n) - 1/alpha(i)*miu_1*A)*V(:,i).*exp(alpha(i)*d1); 
    Ndminus(:,i) = (eye(n,n) + 1/alpha(i)*miu_1*A)*V(:,i).*exp(-alpha(i)*d2); 
end

Nu = [Nuplus,Numinus];   Nd = [Ndplus,Ndminus];   
end

% =======================================================
function Trans = interpolation(n,euni,eps,miu0)
% all angle in degree

thetai = real(asin(real(sqrt(euni/eps)).*sin(acos(miu0)))*180/pi);
theta_uni = acos(miu0)*180/pi;

if euni == eps
    Trans = eye(n);
else
    Trans = zeros(n/2,n/2); 
    for i = 1:n/2
       ub = find(theta_uni>thetai(i), 1 );
       lb = find(theta_uni<thetai(i), 1, 'last' );

       if ub == 1 
           Trans(i,1) = 1;
       elseif lb == 16 
           Trans(i,16) = 1;
       else
       % linear interpolation
       coefa = (thetai(i)-theta_uni(lb))/(theta_uni(ub)-theta_uni(lb));
       Trans(i,lb)= 1 - coefa;  Trans(i,ub) = coefa;
       end
    end
    Trans = [Trans,zeros(16,16);zeros(16,16),Trans];
end
end

% ==========================================================
