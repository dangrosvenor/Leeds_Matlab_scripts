function [NX0,LAMDA]=SingleMom(Q0I2,RHO,str,queryData)

if(strcmp('rain',str) ==1)
    % assignment of variables
    cx=queryData{8}; %523.6;
    alphax=queryData{9}; %2.5;
    nbx=0;
    nax=queryData{20}; %1.1E15;
    dx=queryData{10}; %3.;   
end
%Rna_G=5.0E25,Rnb_G=-4.0,alph_G=2.5 c_G=261.8,d_G=3.0
%Rna_R=1.1E15,Rnb_R=0.0,alph_R=2.5 c_R=523.6,d_R=3.0
LAMDA=real((nax*cx.*gamma(1+alphax+dx)./(RHO.*Q0I2)).^(1./(1+alphax+dx-nbx)));
NX0=real(nax.*LAMDA.^nbx);
%find(NX0 == Inf);
%NX0(ans)=0.;