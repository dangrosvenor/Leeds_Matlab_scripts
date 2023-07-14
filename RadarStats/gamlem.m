function [lam,nx0,nD,mD,D]=gamlem(n,q,rho,type)
%returns lambda and nx0 coefficients in LEM gamma distribution
%usage: [lam,nx0]=gamlem(n,q,rho,type). n is the number concentration in #/kg, 
%q is the mixing ratio in kg/kg, rho is the air density (kg/m^3)
%type is either 'r','s','g', or 'i' for rain, snow, graupel or ice


cx(1)=523.6; %rain
cx(2)=52.36; %snow
cx(3)=366.5; %graupel - changed to this from 261.8 due to tropical graupel modification
cx(4)=104;  %ice

dx=3;

alx(1)=2.5;
alx(2)=2.5;
alx(3)=2.5;
alx(4)=0;

switch type
	case 'r'
        it=1;
	case 's'
        it=2;
	case 'g'
        it=3;
	case 'i'
        it=4;        
end


lam=( n.*cx(it)*gamma(1+alx(it)+dx)./(q.*gamma(1+alx(it))) ).^(1/dx);
a=find(lam==Inf);
lam(a)=0;
nx0=rho.*n.*lam.^(1+alx(it))/gamma(1+alx(it));

minD=0.1e-6;
maxD=20e-1;
nbins=5000;

dD=(log10(maxD)-log10(minD))/nbins;

D=10.^([log10(minD):dD:log10(maxD)]);
D=repmat(D,[length(n) 1]);
D=D';

mD=cx(it).*D.^dx; %M(D)
nD=repmat(nx0',[1 size(D,1)])'.*D.^alx(it).*exp(-repmat(lam',[1 size(D,1)])'.*D); %dN/dD(D) (#/m/m3(air))

D=D(:,1);
    