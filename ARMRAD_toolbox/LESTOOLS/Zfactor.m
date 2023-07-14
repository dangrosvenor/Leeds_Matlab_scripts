function [Zx]=Zfactor(NX0,LAMDA,RHO,str,queryData)
% compute reflectivity factor for species x

if(strcmp('ice',str) ==1)
    % assignment of variables
    cx=queryData{17}; %104.;
    alphax=queryData{18}; %0.;
    dx=queryData{19}; %3.;   
end
if(strcmp('graupel',str) ==1)
    % assignment of variables
    cx=queryData{14}; %366.5;
    alphax=queryData{15}; %2.5;
    dx=queryData{16}; %3.;   
end
if(strcmp('snow',str) ==1)
    % assignment of variables
    cx=queryData{11}; %52.36;
    alphax=queryData{12}; %2.5;
    dx=queryData{13}; %3.;    
end
%Zx=real((NX0.*gamma(1+alphax+6)./(RHO.*LAMDA.^(1+alphax+6)))./((1E-3).^6));
Zx=(6.*cx./pi./1000).^2.*real((NX0.*gamma(1+alphax+2.*dx)./...
    (RHO.*LAMDA.^(1+alphax+2.*dx)))./((1E-3).^6));

% On the z factor...
% int N(D_m)D_m^6 dD, but for the current configuration, dx=3 so particles are spherical
% 
% pi./6.*Dm^3.*1000=D^3.*cx
% Dm^6 = (D^3.*cx.*6./pi./1000)^2