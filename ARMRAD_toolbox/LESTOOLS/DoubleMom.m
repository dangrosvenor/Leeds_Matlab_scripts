function [NX0,LAMDA]=DoubleMom(NQ0I2,Q0I2,RHO,str,queryData)

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

%LAMDA=real((NQ0I2.*cx.*gamma(1+alphax+dx)./(RHO.*Q0I2)).^(1./dx));
LAMDA=real((NQ0I2.*cx.*gamma(1+alphax+dx)./(Q0I2.*gamma(1+alphax))).^(1./dx));
NX0=real(RHO.*NQ0I2.*(LAMDA.^(1+alphax))./gamma(1+alphax));
%find(NX0 == Inf);
%NX0(ans)=0.;