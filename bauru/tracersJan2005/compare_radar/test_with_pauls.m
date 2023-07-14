%test radar calcs against Paul's

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


ihm=1;
alphax=alx(ihm);

ih1=1;
ih2=250;



[r,c,p]=size(TwoD.P);
PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
Tempera=TwoD.TH2.*((PRESS)./100000).^0.286;

% Calculate density
RHO=(PRESS./287)./Tempera;




hm='graupel';
hm='snow';
%hm='ice';

switch hm
case 'graupel'
    ihm=3;
	Q=TwoD.Q(:,:,5);
	NQ=TwoD.Q(:,:,8);
    
	lamgra=( TwoD.Q(ih1:ih2,:,8).*cx(3)*gamma(1+alx(3)+dx)./(TwoD.Q(ih1:ih2,:,5).*gamma(1+alx(3))) ).^(1/dx);
	%ii=find(TwoD.Q(ih1:ih2,:,5)<1e-8 | TwoD.Q(ih1:ih2,:,8)<1);
	%lamgra(ii)=1e10;
	nxogra=RHO.*TwoD.Q(ih1:ih2,:,8).*lamgra.^(1+alx(3))/gamma(1+alx(3));
	z=1e18*0.224*nxogra.*lamgra.^(-7-alx(3))*gamma(7+alx(3));
	%iii=find(zgra<=0);
	%zgra(iii)=1e-10;
    


case 'snow'
    ihm=2;
	Q=TwoD.Q(:,:,4);
	NQ=TwoD.Q(:,:,9);
    
   lamsnow=( TwoD.Q(ih1:ih2,:,9).*cx(2)*gamma(1+alx(2)+dx)./(TwoD.Q(ih1:ih2,:,4).*gamma(1+alx(2))) ).^(1/dx);
 %  ii=find(TwoD.Q(ih1:ih2,:,4)<1e-8 | TwoD.Q(ih1:ih2,:,9)<1 );
 %  lamsnow(ii)=1e10;
   nxosnow=RHO.*TwoD.Q(ih1:ih2,:,9).*lamsnow.^(1+alx(2))/gamma(1+alx(2));
   z=1e18*0.224*nxosnow.*lamsnow.^(-7-alx(2))*gamma(7+alx(2));
 %  iii=find(zsnow<=0);
 %  zsnow(iii)=1e-10;
   

    
    
case 'ice'
    ihm=4;
	Q=TwoD.Q(:,:,6);
	NQ=TwoD.Q(:,:,7);    
    
	lamice=( TwoD.Q(ih1:ih2,:,7).*cx(4)*gamma(1+alx(4)+dx)./(TwoD.Q(ih1:ih2,:,6).*gamma(1+alx(4))) ).^(1/dx);
%   ii=find(TwoD.Q(ih1:ih2,:,6)<1e-8 | TwoD.Q(ih1:ih2,:,7)<1 );
%   lamice(ii)=1e10;
   nxoice=RHO.*TwoD.Q(ih1:ih2,:,7).*lamice.^(1+alx(4))/gamma(1+alx(4));
   z=1e18*0.224*nxoice.*lamice.^(-7-alx(4))*gamma(7+alx(4));
 %  iii=find(zice<=0);
 %  zice(iii)=1e-10;
      
end


alphax=alx(ihm)
cx=cx(ihm)


[LAMBDA,NX0]=Pauls_lam_nxo(Q,NQ,alphax,cx,dx,RHO); %these are ok


Zx=0.19*(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./...
    (RHO.*LAMBDA.^(1+alphax+2.*dx)))./((1E-3).^6));







   
   