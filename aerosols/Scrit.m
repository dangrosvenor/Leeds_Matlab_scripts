function [S,R,Mav]=Scrit(T,rd)
%Calculates critical supersaturation in % and critical radius in m as
%function of dry radius
%T in K and rd=dry radius in m.

%calculations are based in cgs units

rd=rd.*100; %convert rd from m to cm

%Ammonium sulphate
mu=3; %no. ions of salt that dissociate
Ms=132.1; %molecular weight of salt g/mol % ammonium sulphate (NH4)*2 SO4 = 132 g/mol
rhoS=1.769; %density of salt g/cm^3 %ammonium sulphate=1.77 
Mw=18.02; %molecular weight of water

ms=4/3*pi.*rd.^3*rhoS; %mass of salt in g 


%ms=1e-14;

A=3.3e-5./T; %units = cm
B=4.3*mu.*ms/Ms;

S=100*2/3*sqrt(A.^3/3./B); %critical supersaturation
R=sqrt(3.*B./A); %critical radius

%S=0.89;
%R=A.*( 4*Ms*1./(27*Mw*rhoS*mu.*(S./100).^2) ).^(1/3);


R=R./100; %convert to m

%S=1 + A./rd - B./rd.^3;




