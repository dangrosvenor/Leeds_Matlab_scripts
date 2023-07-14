%ruther_cross1.m - program to do the scattering cross-section versus scattering
%angle in Rutherford Scattering
clear; warning off;
m=1;                     %projectile mass in units of alpha particle mass
vb=0.01965;              %velocity units (a_b/tau_b) (in inits of c=light speed)
ma=3730e6;               %alpha particle mass energy in eV
Ene=5e6;                 %initial projectile energy in eV
v0=sqrt(2*Ene/m/ma)/vb;  %initial speed in units of vb
za=2;zt=79;              %za=projectile, zt=target charges (2->alpha, 79->gold)
K=za*zt;                 %dimensionless force constant
thsc=.4:.1:pi;
%Scatt. Cross Section
sigma=K^2*2*pi*sin(thsc)./sin(thsc/2).^4/4/m^2/v0^4;
plot(thsc,sigma,'k')
xlabel('\Theta (Radians)','FontSize',14)
ylabel('\sigma(\Theta) (a_b^2)','FontSize',14)
title('Rutherford Scattering Cross-Section versus Scattering Angle')
