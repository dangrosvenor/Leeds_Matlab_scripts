% 1=rain, 2=snow, 3=graupel, 4=ice (constants - see p.19 of les doc part two.

function ztot=Radar(Grid,TwoD)

cx(1)=523.6;
cx(2)=52.36;
cx(3)=261.8;
cx(4)=104;

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



[PRE,RHO,TEMP,QSAT]=Press2_3(Grid,TwoD);
   
   ii=find(TwoD.Q(:,:,3)<1e-8);
   lamrain=(nax(1)*cx(1)*gamma(1+alx(1)+dx)./(RHO.*TwoD.Q(:,:,3))).^(1/(1+alx(1)+dx-nbx(1))); %single moment 
   lamrain(ii)=1e10;
   nxorain=nax(1).*lamrain.^nbx(1);
   zrain=1e18*nxorain.*lamrain.^(-7-alx(1)).*gamma(7+alx(1));
   iii=find(zrain<=0);
   zrain(iii)=1e-10;
   
   lamice=( TwoD.Q(:,:,7).*cx(4)*gamma(1+alx(4)+dx)./(TwoD.Q(:,:,6).*gamma(1+alx(4))) ).^(1/dx);
   ii=find(TwoD.Q(:,:,6)<1e-8);
   lamice(ii)=1e10;
   nxoice=RHO.*TwoD.Q(:,:,7).*lamice.^(1+alx(4))/gamma(1+alx(4));
   zice=1e18*0.224*nxoice.*lamice.^(-7-alx(4))*gamma(7+alx(4));
   iii=find(zice<=0);
   zice(iii)=1e-10;
   
   lamsnow=( TwoD.Q(:,:,9).*cx(2)*gamma(1+alx(2)+dx)./(TwoD.Q(:,:,4).*gamma(1+alx(2))) ).^(1/dx);
   ii=find(TwoD.Q(:,:,4)<1e-8);
   lamsnow(ii)=1e10;
   nxosnow=RHO.*TwoD.Q(:,:,9).*lamsnow.^(1+alx(2))/gamma(1+alx(2));
   zsnow=1e18*0.224*nxosnow.*lamsnow.^(-7-alx(2))*gamma(7+alx(2));
   iii=find(zsnow<=0);
   zsnow(iii)=1e-10;
   
   lamgra=( TwoD.Q(:,:,8).*cx(3)*gamma(1+alx(3)+dx)./(TwoD.Q(:,:,5).*gamma(1+alx(3))) ).^(1/dx);
   ii=find(TwoD.Q(:,:,5)<1e-8);
   lamgra(ii)=1e10;
   nxogra=RHO.*TwoD.Q(:,:,8).*lamgra.^(1+alx(3))/gamma(1+alx(3));
   zgra=1e18*0.224*nxogra.*lamgra.^(-7-alx(3))*gamma(7+alx(3));
   iii=find(zgra<=0);
   zgra(iii)=1e-10;
   
   ztot=(zrain+zsnow+zgra+zice);
   
  
   %break
   %ii=find(isnan(zrain));
   %zrain(ii)=1e-10;
   %ii=find(isnan(zsnow));
   %zsnow(ii)=1e-10;
   %ii=find(isnan(zgra));
   %zgra(ii)=1e-10;
   %ii=find(isnan(zice));
   %zice(ii)=1e-10;
   
   
        
