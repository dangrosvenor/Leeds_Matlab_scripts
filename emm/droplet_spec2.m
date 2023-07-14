 Dg=0.090870e-6; %mode diameter
 sig2=6.223699;
 sig=20*Dg/2;
 
 Nbins=10;
 N=1*569.64368e6; %total no. conc m^-3
 
 
 M=0.026e-3; %actual from EMM kg/m3
 dbar= 2*( M/N/(4/3*pi*1000) )^(1/3); %mean diameter of drop of mass M/N
 
Dg=dbar;
 
sig=0.05*Dg;
 
Dlow=Dg/10/sig2;
Dhigh=Dg*10*sig2;
Dspacing=( log(Dhigh)-log(Dlow) )/Nbins; %log spacing

logD=[log(Dlow):Dspacing:log(Dhigh)]; %equally spaced in log space

Dmid=(logD(2:end) + logD(1:end-1) )/2;
Dmid=exp(Dmid);



%	nd=lognormal(Dmid,N,sig,Dg)*Dspacing; %total no. in each bin
nd = N/sqrt(2*pi) * exp( - ((Dmid/2-Dg/2)/2/sig).^2 /2 );

%figure;
%plot(Dmid*1e6,nd,'bx-');
m=1000*nd.*4./3.*pi.*(Dmid/2).^3;

%set(gca,'xlim',[0 1]);




%dbar=20e-6;
N=5000e6;
sig=6;

M=0.026e-3; %actual from EMM kg/m3
Dg = 2*( M/N/(4/3*pi*1000) )^(1/3); 

%Dg=25e-6;


Nbins=10;
Dlow=Dg/10;
Dhigh=Dg*10;
Dspacing=( log10(Dhigh)-log10(Dlow) )/Nbins; %log spacing

Dmid=( log10(3*Dg)-10*Dspacing:Dspacing:log10(3*Dg) );
Dmid=10.^(Dmid);

nd=lognormal(Dmid,N,sig,Dg)*Dspacing; %total no. in each bin

sig=0.3*Dg;
nd = N/sqrt(2*pi) * exp( - ((Dmid/2-Dg/2)/2/sig).^2 /2 );

m=1000*nd.*4./3.*pi.*(Dmid/2).^3;

figure;
plot(Dmid*1e6,nd,'bx-');




lwc=0.0026;
N=500e6;


