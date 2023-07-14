function rate_bigg_final=microphysics_WRF_biggs(tc,QC3D)
%final_contact_rate = microphysics_WRF_biggs(tc,QC3D)
%gives the rate of contact IN production in #/L/s - multiply by a time to get
%a total number
%tc  - temperature celsius
%QC3D - cloud MR (g/m^3)

%P - pressure in mb
P=1000; %not actually a function of P either - must cancel out

T3D = tc + 273.15;
P = P*100;

RHO = P./(287.15*T3D); %this is the equation of state for RHO
%RHO=density(P,T3D); %but is a function of density if pressure is held constant...
%must take into account the ideal gas equation somewhere

QC3D = QC3D./RHO/1000;   %convert to kg/kg
%RHO=0.2;


NDCNST=250; %currently set constant (250 per cc) in version 3.0.1.1 of WRF
NC3D=NDCNST.*1.E6./RHO; %convert to #/kg - but note, Biggs not dependent on this - see below

RHOW=997;
CONS26=pi/6.*RHOW;


%DUM3 = P./(287.15*T3D); %this is just calcualting RHO 287.15 = R/mu = 8.314472/28.97e-3
PGAM=0.0005714*(NC3D./1.E6./RHO)+0.2714;
PGAM=1./(PGAM.^2)-1;
PGAM=max(PGAM,2);
PGAM=min(PGAM,10);


CDIST1 = NC3D./gamma(PGAM+1);
LAMC=(CONS26.*NC3D.*GAMMA(PGAM+4)./(QC3D.*GAMMA(PGAM+1))).^(1/3);

LAMMIN = (PGAM+1.)/60.e-6;
LAMMAX = (PGAM+1.)/1.e-6;

iless=find(LAMC<LAMMIN);
LAMC(iless) = LAMMIN(iless);

imore=find(LAMC>LAMMAX);
LAMC(imore) = LAMMAX(imore);


AIMM=0.66;
BIMM=100;
CONS40=pi/6.*BIMM;

%%% immersion nucleation (of cloud droplets to form ice) - Bigg's

rate_bigg=CONS40*exp(log(CDIST1)+log(gamma(PGAM+4))-3*log(LAMC)).*exp(AIMM*(273.15-T3D));
rate_bigg_final = rate_bigg.*RHO/1000;


% N.B. doesn't depend on droplet number as cancels out in rate_bigg formula
% (from formulae for PGAM, CDIST1 and LAMC)

% aside - mass of ice frozen from Biggs does depend on droplet size
%CONS39=pi*pi/36.*RHOW*BIMM;
%MNUCCC = CONS39*exp(log(CDIST1)+log(gamma(7.+PGAM))-6.*log(LAMC)).*exp(AIMM*(273.15-T3D));