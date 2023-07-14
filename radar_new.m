% 1=rain, 2=snow, 3=graupel, 4=ice (constants - see p.19 of les doc part two.

function ztot=Radar_new(Grid,TwoD,ih1,ih2)


cx(1)=523.6; %rain
cx(2)=52.36; %snow
cx(3)=366.5;   %261.8; - changed to this due to tropical graupel modification
cx(4)=104; %ice

dx=3;

nax(1)=1.1e15;
nax(2)=2.0e27;
nax(3)=5e25;

nbx(1)=0;
nbx(2)=-3.5;
nbx(3)=-4;

alx(1)=2.5;
alx(2)=2.5;
alx(3)=2.5;
alx(4)=0;



if ~exist('ih1')
    ih1=1;
    ih2=length(Grid.Z);
end


[r,c,p]=size(TwoD.P);
PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
Tempera=TwoD.TH2.*((PRESS)./100000).^0.286;

% Calculate density
RHO=(PRESS./287)./Tempera;
RHO=RHO(ih1:ih2,:);




ztot=0;

for it=1:4
    
	switch it
	case 4 %ice
        lab='Ice';
        im=6;
        in=7;
	case 2 %snow
        lab='Snow';
        im=4;
        in=9;
	case 3 %graupel
        lab='Graupel';
        im=5;
        in=8;
	case 1 %rain
        lab='Rain';
        im=3;
        in=3;        
	end
    
    n=TwoD.Q(ih1:ih2,:,in);
    q=TwoD.Q(ih1:ih2,:,im);
    
	if strmatch(lab,'Rain');
       ii=find(q<1e-8);
	   lam=( nax(it)*cx(it)*gamma(1+alx(it)+dx)./(RHO.*q) ).^(1/(1+alx(it)+dx-nbx(it))); %single moment        
       lam(ii)=1e10;
       nxo=nax(it).*lam.^nbx(it);
       Z(it).z=paul_zr(cx(it),nxo,lam,alx(it),dx);  %no 0.19 for liquid (rain)
   else          
       lam=( n.*cx(it)*gamma(1+alx(it)+dx)./(q.*gamma(1+alx(it))) ).^(1/dx);
       ii=find(q < 1e-8 | n < 1 );
       lam(ii)=1e10;
       nxo=RHO.*n.*lam.^(1+alx(it))/gamma(1+alx(it));
       Z(it).z=0.19*paul_zr(cx(it),nxo,lam,alx(it),dx);    %0.19 factor for ice (= K^2/0.93)    
	end
    
  
   
    ztot=ztot+Z(it).z;   
end
   
        
ztot=ztot;
   
   
   
   