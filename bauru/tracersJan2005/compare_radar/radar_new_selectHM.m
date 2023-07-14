% 1=rain, 2=snow, 3=graupel, 4=ice (constants - see p.19 of les doc part two.

function z=Radar_new_selectHM(Grid,n,q,ih1,ih2,it,vol)


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


% [r,c,p]=size(TwoD.P);
% PRESS=permute(repmat(Grid.PREFN,[1 1 c]),[1 3 2])+TwoD.P;
% Tempera=TwoD.TH2.*((PRESS)./100000).^0.286;
% 
% % Calculate density
% RHO=(PRESS./287)./Tempera;
%RHO=RHO(ih1:ih2,:,:);




ztot=0;


    
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
        in=99;   %so will be ignored      
	end

kstep=50;    
kchunks=[1:kstep:250+1];
for k=1:length(kchunks)-1
    k
    kselect=[kchunks(k):kchunks(k+1)-1];
        
    RHO=repmat(Grid.RHON(kselect),[1 size(q,1) size(q,2)]);  %using this approximation to rho (instead of 2d pressure and temp
    RHO=permute(RHO,[2 3 1]); 
 
 
	if strmatch(lab,'Rain');
       ii=find(q(:,:,kselect)<1e-8);
%	   lam=( nax(it)*cx(it)*gamma(1+alx(it)+dx)./(Grid.RHON(k).*q(:,:,k)) ).^(1/(1+alx(it)+dx-nbx(it))); %single moment        
	   lam=( nax(it)*cx(it)*gamma(1+alx(it)+dx)./(RHO.*q(:,:,kselect)) ).^(1/(1+alx(it)+dx-nbx(it))); %single moment        
       
       lam(ii)=1e10;
%       nxo=nax(it).*lam.^nbx(it);
       z(:,:,kselect)=paul_zr(cx(it),nax(it).*lam.^nbx(it),lam,alx(it),dx);  %note no 0.19 factor for rain
       
   else  
       if length(vol)==0
           c=cx(it);
       else
           c=pi/6.*q(:,:,kselect)./vol(:,:,kselect); %mass constant from grauple density
       end %replacing cx(it) with pi*rho_graupel/6
       
      lam=( n(:,:,kselect).*c*gamma(1+alx(it)+dx)./(q(:,:,kselect).*gamma(1+alx(it))) ).^(1/dx);           


       ii=find(q(:,:,kselect) < 1e-8 | n(:,:,kselect) <= 1 );
       lam(ii)=1e10;
       
       if length(vol)~=0
           c(ii)=0;
       end

       %nxo=Grid.RHON(k).*n.*lam.^(1+alx(it))/gamma(1+alx(it));
%       z(:,:,k)=0.19*paul_zr(cx(it),Grid.RHON(k).*n(:,:,k).*lam.^(1+alx(it))/gamma(1+alx(it)),lam,alx(it),dx);       
       z(:,:,kselect)=0.19*paul_zr(c,RHO.*n(:,:,kselect).*lam.^(1+alx(it))/gamma(1+alx(it)),lam,alx(it),dx);       

	end

end %k
   




   
        

   
   
   
   