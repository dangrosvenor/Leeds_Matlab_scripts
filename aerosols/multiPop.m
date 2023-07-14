%calculates aerosol distributions. bargraph_aerosol plots bar graph of distributions
clear N Sc logD

T=276; %283
Nbins=1e4;


aerosol_case='urban';
aerosol_case='background';
%aerosol_case='marine';
%aerosol_case='clean';
%aerosol_case='SM2';
%aerosol_case='urban';

switch aerosol_case
    
case 'background'
%average background nulcei mode
Dg=0.016e-6; %mode diameter
sig=1.7;
Ntot=6400e6; %total no. conc m^-3

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(1).n=N(1).n/Dspacing;


%average background accumulation mode
Dg=0.076e-6; %mode diameter
sig=2.0;
Ntot=2300e6; %total no. conc m^-3

[Sc(2).s,Dcrit,N(2).n,Dspacing,logD(2).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(2).n=N(2).n/Dspacing;

case 'urban'
%urban
%nulei mode
Dg=0.014e-6;
Ntot=106000e6;
sig=1.8;

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(1).n=N(1).n/Dspacing;

%accumulation mode
Dg=0.054e-6;
sig=2.16
Ntot=32000e6;

[Sc(2).s,Dcrit,N(2).n,Dspacing,logD(2).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(2).n=N(2).n/Dspacing;

case 'marine'
%urban

%nuclei mode
Dg=0.01e-6;
Ntot=340e6;
sig=1.6;

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(1).n=N(1).n/Dspacing;

%accumulation mode
Dg=0.07e-6;
sig=2.;
Ntot=60e6;

[Sc(2).s,Dcrit,N(2).n,Dspacing,logD(2).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(2).n=N(2).n/Dspacing;

case 'clean'
%average background nulcei mode
Dg=0.016e-6; %mode diameter
sig=1.6;
Ntot=1000e6; %total no. conc m^-3

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(1).n=N(1).n/Dspacing;


%average background accumulation mode
Dg=0.0686e-6; %mode diameter
sig=2.1;
Ntot=800e6; %total no. conc m^-3

[Sc(2).s,Dcrit,N(2).n,Dspacing,logD(2).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(2).n=N(2).n/Dspacing;

case 'SM2'
%average background nulcei mode
Dg=0.02e-6; %mode diameter
sig=2.5;
Ntot=1000e6; %total no. conc m^-3

[Sc(1).s,Dcrit,N(1).n,Dspacing,logD(1).d]=makebinsNenes(Dg,sig,Ntot,Nbins,T);
dN(1).n=N(1).n/Dspacing;

Sc(2).s=Sc(1).s;
N(2).n=zeros(size(N(1).n));
end

%Combined_Sdist3;  %not necessary for the later parts of this script

snew30=snew;


nsbins=333;
Sstart=min([Sc(1).s Sc(2).s]);
Send=max([Sc(1).s(end) Sc(2).s(end)]);
%Send=5.5/100;

Sspacing=( log(Send) - log(Sstart) ) / nsbins;
snew2=[log(Sstart):Sspacing:log(Send)];
snew2=exp(snew2);

Ntot=zeros(length(snew2)-1);
for j=1:length(Sc)
    N(j).n(end+1)=0;
	for i=1:length(snew2)-1
        is=find( Sc(j).s>snew2(i) & Sc(j).s<=snew2(i+1) );
        Ntot(i)=Ntot(i)+sum(N(j).n(is));
    end
end

