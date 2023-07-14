%ruthercm_cross.m - program to do the scattering crossection versus scattering
%angle in Rutherford Scattering - modified version of Ch10's ruther_cross1.m
%here we include CM motion effect and compare to the regular case. The projectile
%particle here is gold, and the target particle is also gold.
clear; warning off;
m=1;                     %projectile mass in units of gold particle mass
vb=0.0044;               %velocity units (a_b/tau_b) (in inits of c=light speed)
ma=197*939e6;            %Gold particle mass energy in eV
Ene=5e6;                 %initial projectile energy in eV
v0=sqrt(2*Ene/m/ma)/vb;  %initial speed in units of vb
za=79;zt=79;             %za=projectile, zt=target charges (79->Gold)
K=za*zt;                 %dimensionless force constant
thsc=.4:.05:pi/2;
%Scatt. Crossection - stationary target
sigma=K^2*2*pi*sin(thsc)./sin(thsc/2).^4/4/m^2/v0^4;
%Scatt. Crossection - recoiling target
sigma_cm=K^2*8*pi*sin(thsc).*cos(thsc)./sin(thsc).^4/m^2/v0^4;
%plot(thsc,sigma,'ko',thsc,sigma_cm,'b.')
semilogy(thsc,sigma,'ko',thsc,sigma_cm,'b.')
xlabel('\Theta (Radians)','FontSize',14)
ylabel('\sigma(\Theta) (a_b^2)','FontSize',14)
str=cat(2,'Gold on Gold Scattering Crossection versus Scattering Angle');
title(str,'FontSize',12)
h=legend('Stationary Target','Recoiling Target'); set(h,'FontSize',14)