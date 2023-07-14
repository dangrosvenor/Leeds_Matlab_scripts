%ruther_cross2.m - program to do the scattering cross-section versus atomic
%number in Rutherford Scattering
clear; warning off;
m=1;                     %projectile mass in units of alpha particle mass
vb=0.01965;              %velocity units (a_b/tau_b) (in inits of c=light speed)
ma=3730e6;               %alpha particle mass energy in eV
Ene=5e6;                 %initial projectile energy in eV
v0=sqrt(2*Ene/m/ma)/vb;  %initial speed in units of vb
b=20;                    %impact parameter
za=2;                    %za=projectile(2->alpha)
ztmin=13; ztmax=79;      %from Aluminum to Gold
j=0;                     %initial counter
for iz=ztmin:2:ztmax
    j=j+1;
    zt(j)=iz;            %target array
    K=za*zt(j);          %dimensionless force constant
    ye=-sign(K)*sqrt(1+m^2*v0^4*b^2/K^2);% eccentricity
% can find min.max angle limits for the assymptotes of r for this ye
% these occur at the zeros of the denominator of the r(theta) equation
    %thmin=fzero(inline('ye*cos(x)+1'),[-pi/2,0],[],ye)+1.e-5;
    %thmax=fzero(inline('ye*cos(x)+1'),[0,pi/2],[],ye)-1.e-5;
    %th(j)=thmax-thmin;   %angle from +x axis
    %thsc(j)=pi-th(j);    %total scattering angle produced by target species
%best is to find alpha angle assymptote - faster way, done only once
    alpha=fzero(inline('ye*cos(x)+1'),[0,pi/2],[],ye)-1.e-5;
    thsc(j)=pi-2*alpha;   %total scattering angle produced by target species
%Scatt. Cross-Section   
    sigma(j)=K^2*2*pi*sin(thsc(j))/sin(thsc(j)/2)^4/4/m^2/v0^4;
end
str=cat(2,'E_k= ',num2str(Ene,3),' eV',', b= ',num2str(b,3),' a_b');
subplot(2,2,1)
plot(zt(:),thsc(:),'k.','MarkerSize',5)
xlabel('Z_t','FontSize',14),ylabel('\Theta (Radians)','FontSize',14)
title('\alpha Particle Scattering for various Targets','FontSize',12)
subplot(2,2,2)
plot(zt(:),sigma(:),'b.','MarkerSize',5)
xlabel('Z_t','FontSize',14),ylabel('\sigma(\Theta) (a_b^2)','FontSize',14)
title(str,'FontSize',12)
subplot(2,2,3)
plot(thsc(:),sigma(:),'r.','MarkerSize',5)
xlabel('\Theta (Radians)','FontSize',14)
ylabel('\sigma(\Theta) (a_b^2)','FontSize',14)
