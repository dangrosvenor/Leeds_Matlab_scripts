function [Zx]=ZfactorTripleMom(NX0,LAMDA,RHO,str,POINT,TwoD,queryData,run2)
% compute reflectivity factor for species x

if(strcmp('graupel',str) ==1)
    Q=find([queryData{:,7}]==run2);
    % assignment of variables
    cx=queryData{Q,14}; %366.5;
    alphax=queryData{Q,15}; %2.5;
    dx=queryData{Q,16}; %3.;   
  
    if queryData{Q(1),5}==0 %no prognostic graupel
        Zx=(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./(RHO.*LAMDA.^(1+alphax+2.*dx)))./((1E-3).^6));
        numero=0
    else %prognostic graupel       
       rho(:,:,:)= TwoD.Q(:,:,POINT.IQG(1),:)./TwoD.Q(:,:,POINT.IQGV(1),:);
       Zx=(rho./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./(LAMDA.^(1+alphax+2.*dx)))./((1E-3).^6));
       numero=1;
    end


elseif(strcmp('ice',str) ==1)
        % assignment of variables
        cx=104.;
        alphax=0.;
        dx=3.;   
        %Zx=real((NX0.*gamma(1+alphax+6)./(RHO.*LAMDA.^(1+alphax+6)))./((1E-3).^6));
        Zx=(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./(RHO.*LAMDA.^(1+alphax+2.*dx)))./((1E-3).^6));
        
elseif(strcmp('snow',str) ==1)
        % assignment of variables
        cx=52.36;
        alphax=2.5;
        dx=3.;    
        %Zx=real((NX0.*gamma(1+alphax+6)./(RHO.*LAMDA.^(1+alphax+6)))./((1E-3).^6));
        Zx=(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./(RHO.*LAMDA.^(1+alphax+2.*dx)))./((1E-3).^6));
    
end

Zx(find(isinf(Zx(:))))=NaN;