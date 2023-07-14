clear N Sc logD Dg sig Ntot Dlow Dhigh M

%ammonium sulphate properties
mu=3; %no. ions of salt that dissociate
Ms=132.1; %molecular weight of salt g/mol
rhoS=1.769; %density of salt g/cm^3 %ammonium sulphate=1.77 
Mw=18.02; %molecular weight of water


%define various modes

%average background nulcei mode
Dg(1)=0.016e-6; %mode diameter (m)
sig(1)=1.7;
Ntot(1)=6400e6; %total no. conc m^-3

%average background accumulation mode
Dg(2)=0.076e-6; %mode diameter (m)
sig(2)=2.0;
Ntot(2)=2300e6; %total no. conc m^-3


%create overall number distribution in doubling mass bins

for i=1:length(Dg) %run over all distributions to be added
    Dlow(i)=Dg(i)/10/sig(i); %max and min diameter based on population properties as in Nenes&Seinfeld paper
    Dhigh(i)=Dg(i)*10*sig(i);
    %convert from diameter to mass
    Mlow(i)=4/3*pi.* (Dlow(i)*100/2).^3 * rhoS/1000; %mass of salt in kg
    Mhigh(i)=4/3*pi.* (Dhigh(i)*100/2).^3 * rhoS/1000; %mass of salt in kg
end

T=298; %define T quite high so that max critical supersat should cover all environments (at lower temps supersat will increase above 5.5%)
A=3.3e-5./T; %units = cm
smax=5.5; %max critical supersat required
B=A^3/3 * (200/3/smax)^2;
msmin=B*Ms/4.3/mu/1000; %corresponding max mass

msmax=max(Mhigh);
%msmin=min(Mlow); %find extreme values over all distributions

n=40;
%n=ceil(1+log10(msmax/msmin)/log10(2)); %no. bins for mass doubling from min mass to max mass
X=10^(1/(n-1) * log10(msmax/msmin)); %factor to increase mass by to get n bins

M(1)=msmin;
for i=2:n
    M(i)=M(i-1)*X; %create mass bins with doubling
    D2=log( (M(i)*3/4/pi/rhoS*1000)^(1/3)*2/100 ); %large diameter in m
    D1=log( (M(i-1)*3/4/pi/rhoS*1000)^(1/3)*2/100 ); %small diameter in m
    N(i-1)=0; %running total of number in mass bin
    Mtot(i-1)=0;
    for j=1:length(Dg) %run over all distributions to be added
        N(i-1) = N(i-1) + quadl(@lognormal2,D1,D2,[],[],Ntot(j),sig(j),Dg(j)); %integrate lognormal distribution over mass interval (converted to log diameter)
        Mtot(i-1) = Mtot(i-1) + quadl(@DmtotDm,M(i-1),M(i),[],[],Ntot(j),sig(j),Dg(j),rhoS); %add total mass of aerosol in mass bin
    end %quad integrates to a certain tolerance with multiple calls to lognormal2
    Mav(i-1)=Mtot(i-1)/N(i-1); %average aerosol mass in bin
end
    