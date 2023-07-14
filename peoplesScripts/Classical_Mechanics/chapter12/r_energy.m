%r_energy.m  - finds the angular momentum and energy of a rigid body about an
%axis of rotation given the angular speed. Can also find the L and E
%for rotations about the principal axes, given their inertia and directions.
clear; format compact;
h=pi/180;         %factor to convert degrees to radians
mw=sqrt(2);w=[cos(45*h);cos(45*h);0]*mw;  %rot speed given
fprintf('w= %5.3f %5.3f %5.3f, magnitude=%5.3f\n',w,mw)
I=[[1/3,-1/4,0];[-1/4,1/3,0];[0,0,2/3]]   %inertia tensor about origin
L=I*w;T=w'*I*w/2; %angular momentum, energy about the given axis
fprintf('L= %5.3f %5.3f %5.3f, T=%5.3f\n',L,T)
Ip=[[1/12,0,0];[0,7/12,0];[0,0,2/3]]      %principal axes moments
w1=[cos(135*h);cos(135*h);cos(90*h)]*mw ; %1st principal axis rotation
fprintf('1st p axis rotation: w1= %5.3f %5.3f %5.3f\n',w1)
%sqrt(w1'*w1)  %can check that the magnitude of rotation is unchanged
L1=Ip*w1;      %angular momentum about 1st princ. axis
T1=w1'*Ip*w1/2;%kinetic energy about 1st princ. axis
fprintf('L1= %5.3f %5.3f %5.3f, T1=%5.3f\n',L1,T1)
w2=[cos(135*h);cos(45*h);cos(90*h)]*mw;   %2nd principal axis rotation
fprintf('2nd p axis rotation: w2= %5.3f %5.3f %5.3f\n',w2)
%sqrt(w2'*w2)  %can check that the magnitude of rotation is unchanged
L2=Ip*w2;      %angular momentum about 2nd princ. axis
T2=w2'*Ip*w2/2;%kinetic energy about 1st princ. axis
fprintf('L2= %5.3f %5.3f %5.3f, T2=%5.3f\n',L2,T2)
w3=[cos(90*h);cos(90*h);cos(0*h)]*mw;     %3rd principal axis rotation
fprintf('3rd p axis rotation: w3= %5.3f %5.3f %5.3f\n',w3)
%sqrt(w3'*w3)  %can check that the magnitude of rotation is unchanged
L3=Ip*w3;      %angular momentum about 3rd princ. axis
T3=w3'*Ip*w3/2;%kinetic energy about 1st princ. axis
fprintf('L3= %5.3f %5.3f %5.3f, T3=%5.3f\n',L3,T3)
