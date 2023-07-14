function [ztot]=calc_radar_UM_LEM_etc(model,q_in,n_in,RHO)

switch model
    case 'LEM'

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
        
        kfac(1) = 1; %No k-value modification for rain
        kfac(2) = 0.19; %0.19 factor for ice (= K^2/0.93)
        kfac(3) = 0.19; %0.19 factor for ice (= K^2/0.93)
        kfac(4) = 0.19; %0.19 factor for ice (= K^2/0.93)  
        
    case 'UM'

    % Mass diameter relation  m(D) = cx * D^(dx)
        cx(1)=pi*997/6; %rain
        cx(2)=pi*100/6; %snow
        cx(3)=pi*500/6; %graupel
        cx(4)=pi*200/6; %ice
        cx(5)=pi*997/6; %liquid
        
        dx=3;

   % n(D) = nx0 * D^(alx) * exp(-lambda_x * D)
   % For double moment scheme :-
   %    nx0 = rho * N * lambda_x^(1+alx) / gamma(1+alx)
   %    nax and nbx not needed for double moment schemes
%         nax(1)=NaN;
%         nax(2)=NaN;
%         nax(3)=NaN;
% 
%         nbx(1)=NaN;
%         nbx(2)=NaN;
%         nbx(3)=NaN;

        alx(1)=2.5;
        alx(2)=2.5;
        alx(3)=2.5;
        alx(4) = 0;
        alx(5) = 0; %for liquid
        
   % Constant factor for k value
        kfac(1) = 1; %No k-value modification for rain
        kfac(2) = 0.19; %0.19 factor for ice (= K^2/0.93)
        kfac(3) = 0.19; %0.19 factor for ice (= K^2/0.93)
        kfac(4) = 0.19; %0.19 factor for ice (= K^2/0.93)        
        kfac(5) = 1; %No k-value modification for liquid water         
end



ztot=0;

for it=1:length(q_in) %

    q = q_in{it};
    n = n_in{it};





    %                     if strmatch(lab,'Rain single moment');
    %                         ii=find(q<1e-8);
    %                         lam=( nax(it)*cx(it)*gamma(1+alx(it)+dx)./(RHO.*q) ).^(1/(1+alx(it)+dx-nbx(it))); %single moment
    %                         lam(ii)=1e10;
    %                         nxo=nax(it).*lam.^nbx(it);
    %                         Z(it).z=paul_zr(cx(it),nxo,lam,alx(it),dx);  %no 0.19 for liquid (rain)
    %                     else
    lam=( n.*cx(it)*gamma(1+alx(it)+dx)./(q.*gamma(1+alx(it))) ).^(1/dx);
    ii=find(q < 1e-8 | n < 1 );
    lam(ii)=1e10;
    nxo=RHO.*n.*lam.^(1+alx(it))/gamma(1+alx(it));
    Z(it).z = kfac(it)*paul_zr(cx(it),nxo,lam,alx(it),dx);    %0.19 factor for ice (= K^2/0.93), =1 for liquid and rain
    %                    end

    Z(it).z(ii) = 1e-100; %set "zeroes" to a very low value


    ztot=ztot+Z(it).z;
end





