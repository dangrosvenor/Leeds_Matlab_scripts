clear N Sc logD ms

T=276; %283
Nbins=1e4;


%average background nulcei mode
Dg=0.016e-6; %mode diameter
sig=1.7;
Ntot=6400e6; %total no. conc m^-3

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);



%average background accumulation mode
Dg=0.076e-6; %mode diameter
sig=2.0;
Ntot=2300e6; %total no. conc m^-3

[Sc(2).s,Dcrit,N(2).n,Dspacing,logD(2).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);




mu=3; %no. ions of salt that dissociate
Ms=132.1; %molecular weight of salt g/mol
rhoS=1.769; %density of salt g/cm^3 %ammonium sulphate=1.77 
Mw=18.02; %molecular weight of water


rd=exp(logD(1).d).*100; %convert rd from m to cm
ms(1).m=4/3*pi.*rd.^3*rhoS/1000; %mass of salt in kg 

rd=exp(logD(2).d).*100; %convert rd from m to cm
ms(2).m=4/3*pi.*rd.^3*rhoS/1000; %mass of salt in kg

%ms in order of decreasing mass

msmin=ms(1).m(end);
msmax=ms(1).m(1);
for i=2:length(ms)
    msmin=min([msmin,ms(i).m(end)]); %find min and max masses
    msmax=max([msmax,ms(i).m(1)]);
end

n=ceil(1+log10(msmax/msmin)/log10(2)); %no. bins for mass doubling from min mass to max mass

M(1)=msmin;
for i=2:n
    M(i)=M(i-1)*2; %create mass bins with doubling
end





