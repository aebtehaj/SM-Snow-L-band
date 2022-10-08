function [F,B,kappas] = tran12f_FB(bmua1,P11,P22,P33,P44,P34,P43,P12,P21,xi,wi)
%transform from 1-2 to principal for passive for phase matrix of big P

% global wi xi Nquad ndeg
Nquad = length(xi);
ndeg = Nquad/2;
mu  = [-xi(1:ndeg);xi(1:ndeg)];
aac = [wi(1:ndeg);wi(1:ndeg)];

% mu = [mu0; -mu0];
% aac = [wi0; wi0];
% ndeg = length(mu0);
% Nquad = 2*ndeg;

% remember as follows
% first half is upward going  and from small angle to large angle
% second half downward but in ordrer of angle as the first half

phsa = xi*pi;
nphs = Nquad;
nphsh= nphs/2+1;
            
ipt=0;
for idegi=1:Nquad %7107
    tai=acos(mu(idegi));
    cti=cos(tai);
    sti=sin(tai);
    degi(idegi)=tai*180/pi;
    for idegs=1:Nquad  %3101
        tas=acos(mu(idegs));
        cts=cos(tas);
        sts=sin(tas);
        degs(idegs)=tas*180/pi;
        for  iphs=1:nphs  %3102
            ipt=ipt+1;
            phs=phsa(iphs);
            cps=cos(phs);
            sps=sin(phs);
            cbt(ipt)=cti*cts+sti*sts*cps;
        end %3102
    end %3101
end %7107
P11i = interp1(bmua1,P11,cbt,'spline') ;
P22i = interp1(bmua1,P22,cbt,'spline') ;
P12i = interp1(bmua1,P12,cbt,'spline') ;
P21i = interp1(bmua1,P21,cbt,'spline') ;
            
P33i = interp1(bmua1,P33,cbt,'spline') ;
P44i = interp1(bmua1,P44,cbt,'spline') ;
P34i = interp1(bmua1,P34,cbt,'spline') ;
P43i = interp1(bmua1,P43,cbt,'spline') ;

P11A = zeros(Nquad, ndeg);
P12A = zeros(Nquad, ndeg);
P21A = zeros(Nquad, ndeg);
P22A = zeros(Nquad, ndeg);

ipt=0;
for idegi=1:ndeg %7107
    tai=acos(mu(idegi));
    cti=cos(tai);
    sti=sin(tai);
    
    ki=[sti 0 cti];
    vi=[cti 0 -sti];
    hi=[0 1 0];
    
    for idegs=1:Nquad  %3101
        tas=acos(mu(idegs));
        cts=cos(tas);
        sts=sin(tas);
        
        P = zeros(4);
        for  iphs=1:nphs  %3102
            icode=3;
            if  ((iphs==1)&&(idegs==idegi));   icode=1;    end;    % forward scattering
            if  ((iphs==nphsh)&&(abs(idegs-idegi)==ndeg));icode=2;end; % backscattering
            
            phs=phsa(iphs);
            cps=cos(phs);
            sps=sin(phs);
            ipt=ipt+1;
            
            P12f=[P11i(ipt) P12i(ipt) 0 0;P21i(ipt) P22i(ipt) 0 0;...
                0 0 P33i(ipt) P34i(ipt);0 0 P43i(ipt) P44i(ipt)];        
            
            ks=[sts*cps sts*sps cts];
            vs=[cts*cps cts*sps -sts];
            hs=[-sps cps 0];
            
            if  ((icode==1)||(icode==2))
                onei=-hi;
            else
                q=cross(ks,ki);
                absq=sqrt(dot(q,q));
                onei=q/absq;
            end
            oness=onei;
            alphai=atan2(dot(hi,onei),dot(vi,onei));
            alphas=atan2(dot(hs,oness),dot(vs,oness ));
            s2i=sin(alphai)^2;
            c2i=cos(alphai)^2;
            s2ai=sin(2*alphai);
            c2ai=cos(2*alphai);
            s2s=sin(alphas)^2;
            c2s=cos(alphas)^2;
            s2as=sin(2*alphas);
            c2as=cos(2*alphas);
            Pti=[c2i s2i s2ai/2 0;s2i c2i -s2ai/2 0;-s2ai s2ai c2ai 0;0 0 0 1];
            Pts=[c2s s2s -s2as/2 0;s2s c2s s2as/2 0;s2as -s2as c2as 0;0 0 0 1];
            
            P =  P + Pts*P12f*Pti.*wi(iphs).*pi;           
        end %3102
        P11A(idegs,idegi)=P(1,1);
        P12A(idegs,idegi)=P(1,2);
        P21A(idegs,idegi)=P(2,1);
        P22A(idegs,idegi)=P(2,2);

    end %3101
end %7107

% getting the kappas.  The two pol should be the same
kappasv = aac'*(P11A + P21A);
kappash = aac'*(P12A + P22A);

kappas = [kappasv; kappash];

fvv = P11A(1:ndeg,:);
fhv = P21A(1:ndeg,:);
fvh = P12A(1:ndeg,:);
fhh = P22A(1:ndeg,:);

bvv = P11A(1+ndeg:end,:);
bhv = P21A(1+ndeg:end,:);
bvh = P12A(1+ndeg:end,:);
bhh = P22A(1+ndeg:end,:);


F = [fvv fvh;fhv fhh];
B = [bvv bvh;bhv bhh];   

end