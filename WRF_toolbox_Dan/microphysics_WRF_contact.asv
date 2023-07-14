function final_contact_rate = microphysics_WRF_contact(tc,P,QC3D)
%final_contact_rate = microphysics_WRF_contact(tc,P,QC3D)
%gives the rate of contact IN production in #/L/s - multiply by a time to get
%a total number
%tc  - temperature celsius
%P - pressure in mb
%QC3D - cloud MR (g/m^3)

T3D = tc + 273.15;
P = P*100;
RHO=density(P,T3D);
QC3D = QC3D./RHO/1000;  %convert to kg/kg


NDCNST=250; %currently set constant (250 per cc) in version 3.0.1.1 of WRF
NC3D=NDCNST.*1.E6./RHO; %convert to #/kg
tc=T3D-273.15;



RIN=0.1e-6;
CONS37=4.*pi*1.38e-23/(6.*pi*RIN);
RHOW=997;
CONS26=pi/6.*RHOW;


DUM2 = 1.496e-6*T3D.^1.5./(T3D+120);
MU = DUM2./RHO;
NACNT=exp(-2.80+0.262*(-tc))*1000; %per m3
DUM = 7.37.*T3D./(288.*10.*P)/100;
DAP = CONS37.*T3D.*(1.+DUM/RIN)./MU;

DUM3 = P./(287.15*T3D);
PGAM=0.0005714*(NC3D./1.E6./DUM3)+0.2714;
PGAM=1./(PGAM.^2)-1;
PGAM=max(PGAM,2);
PGAM=min(PGAM,10);

CDIST1 = NC3D./gamma(PGAM+1);
LAMC=(CONS26.*NC3D.*GAMMA(PGAM+4)./(QC3D.*GAMMA(PGAM+1))).^(1/3);

% LAMMIN, 60 MICRON DIAMETER
% LAMMAX, 1 MICRON

LAMMIN = (PGAM+1.)/60.e-6;
LAMMAX = (PGAM+1.)/1.e-6;

iless=find(LAMC<LAMMIN);
LAMC(iless) = LAMMIN(iless);
NC3D(iless)= exp(3.*log(LAMC(iless))+log(QC3D(iless))+    ...
    log(GAMMA(PGAM(iless)+1.))-log(GAMMA(PGAM(iless)+4.)))/CONS26;

imore=find(LAMC>LAMMAX);
LAMC(imore) = LAMMAX(imore);
NC3D(imore) = exp(3.*log(LAMC(imore))+log(QC3D(imore))+  ...
    log(GAMMA(PGAM(imore)+1.))-log(GAMMA(PGAM(imore)+4.)))/CONS26;


%NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
%                    GAMMA(PGAM(K)+2.)/                         &
%                    LAMC(K)

%this is a rate in #/kg/s
rate_contact=2.*pi*DAP.*NACNT.*CDIST1.*gamma(PGAM+2)./LAMC;
%DT=0.5*60;
%DT=500;
%final_contact_rate = rate_contact*DT.*RHO/1000; %work out for arbitrary time
final_contact_rate = rate_contact.*RHO/1000; %rate in #/L/s



AIMM=0.66;                            
BIMM=100;
CONS40=pi/6.*BIMM;                           

%%% immersion nucleation (of cloud droplets to form ice) - Bigg's #/kg/s
rate_bigg=CONS40*exp(log(CDIST1)+log(gamma(PGAM+4))-3*log(LAMC)).*exp(AIMM*(273.15-T3D));



