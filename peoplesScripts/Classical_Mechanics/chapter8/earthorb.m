%earthorb.m - script to draw Earth's ellipse with minimum radius rmim and
%eccentricity e given
clear;
m=5.98e24;       %Earth's mass in kg
G=6.67e-11;      %Universal gravitational constant
M=1.991e30;      %Sun's mass in kg
a=1.496e8;       %average Earth-Sun distance in km, or semimajor axis
e=0.017;         %Earth's eccentricity
r0=a*(1-e);      %r0 in km, same as rmin
K=-G*M*m/1000^3; %in km^3.kg/s^2
vc=sqrt(-K/r0/m);%speed in km/s (circular orbit speed)
v0=vc*sqrt(1+e); %speed at r0
th=[0:0.05:2*pi];%range variable
r=r0*(1+e)./(1+e*cos(th));%orbit formula, with ro=rmin
rmax=r0*(v0/vc)^2/(2-(v0/vc)^2);%same as ro*(1+e)/(1-e)
v_rmax=r0*v0/rmax;%using angular momentum conservation
polar(th,r,'k.');hold on
polar(0,r0,'bo')
polar(pi,rmax,'ro')
str1=cat(2,'Earth Orbit: e = ',num2str(e,3),', (r_0, r_{max}) = (',...
    num2str(r0,4),',',num2str(rmax,4),')km');
title(str1,'FontSize',12)
str2=cat(2,' (v_0,v_{rmax})=(',num2str(v0,4),',',num2str(v_rmax,4),') km/s');
xlabel(str2,'FontSize',12)
