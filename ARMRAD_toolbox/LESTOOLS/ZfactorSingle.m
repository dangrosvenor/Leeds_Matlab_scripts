function [Zx]=ZfactorSingle(NX0,LAMDA,RHO,str,queryData)
% compute reflectivity factor for species x

if(strcmp('rain',str) ==1)
    % assignment of variables
    cx=queryData{8}; %523.6;
    alphax=queryData{9}; %2.5;
    dx=queryData{10}; %3.;   
    nbx=0;
    nax=queryData{20}; %1.1E15;
end

%Zx=real(nax.*gamma(1+alphax+6)./(RHO)./(LAMDA.^(1+alphax+6))./((1E-3).^6));
Zx=(6.*cx./pi./1000).^2.*real(nax.*gamma(1+alphax+2.*dx)./(RHO)./(LAMDA.^(1+alphax+2.*dx))./((1E-3).^6));

