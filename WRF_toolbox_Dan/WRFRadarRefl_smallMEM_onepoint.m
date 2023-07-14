function Z=WRFRadarRefl_smallMEM_onepoint(nc,time,scheme,iz,ilat,ilon)

%nc=netcdf('/data/pjc/WRF/WRFV2/test/em_real/wrfout_d01_2006-02-05_18:00:00')
%time=31

p =WRFUserARW(nc,'p',time,999,999,iz); p=p(ilat,ilon);
tc=WRFUserARW(nc,'tc',time,999,999,iz); tc=tc(ilat,ilon);
rho_a=p.*100./287./(tc+273.15);

clear p;

%clear p QC QI;

rho_w=1000;
rho_s=100;
rho_g=400;
rho_i=890;

switch scheme
case 'kessler'
% This uses equations from LEM documentation
QR=WRFUserARW(nc,'QRAIN',time); QR=QR.var;

na=8.e6;
nb=0;
lambda=(na.*cx.*gamma(1+3)./rho_a./QR).^(1./(1+3));
n0=na.*lambda.^nb;

Zr=n0.*(6.*523.6./pi./rho_w).^2./lambda.^(1+2.*3) ...
    .*gamma(1+2.*3)./((1E-3).^6);

Zr(find(isnan(Zr(:))))=1e-30;
Z=abs(Zr);

case 'thompson'
Nt_c=100e6;

mu_c=min(15,1000e6./Nt_c+2);
mu_r=0.;
mu_s=0.6357;
mu_g=0.;
mu_i=0.;

% Snow is sum of two gammas N(D)=M2**4/M3**3 * (Kap0*exp(-M2*Lam0*D/M3)
% + Kap1*(M2/M3)**mu_s*D**mu_s*exp(-M2*Lam1*D/M3))
Kap0=490.6;
Kap1=17.46;
Lam0=20.78;
Lam1=3.29;

% Y-intercept for rain and graupel. vary depending on mixing ratio.
gonv_min=1e4;
gonv_max=1e6;
ronv_min=2e6;
ronv_max=2e9;
ronv_sl =1./4;
ronv_r0 =0.1e-3;
ronv_c0 =ronv_sl./ronv_r0;
ronv_c1 =(ronv_max-ronv_min)*0.5;
ronv_c2 =(ronv_max+ronv_min)*0.5;

% Mass power law relations: mass = am*D^bm
am_r=pi*rho_w/6;
bm_r=3;
am_s=0.069;
bm_s=2;
am_g=pi*rho_g/6;
bm_g=3;
am_i=0.069;
bm_i=2;

% Snow
QS=nc{'QSNOW'}(time,iz,ilat,ilon);
%QI=WRFUserARW(nc,'QICE',time);  QI=QI.var;
%QC=WRFUserARW(nc,'QCLOUD',time);QC=QC.var;

oams=1./am_s;
smo2=(QS.*rho_a).*oams;   %big

% loga_=sa(1)+sa(2).*tc0+sa(3)+sa(4).*tc0+sa(5).*tc0.*tc0...
%     +sa(6)+sa(7).*tc0.*tc0+sa(8).*tc0+sa(9).*tc0.*tc0.*tc0+sa(10);
% a_=10.^loga_;
% b_=sb(1)+sb(2).*tc0+sb(3)+sb(4).*tc0+sb(5).*tc0.*tc0+sb(6)+sb(7).*tc0.*tc0 ...
%     +sb(8).*tc0+sb(9).*tc0.*tc0.*tc0;+sb(10);
% smo1=a_.*smo2.^b_;

% Snow is sum of two gammas N(D)=M2**4/M3**3 * (Kap0*exp(-M2*Lam0*D/M3)
% + Kap1*(M2/M3)**mu_s*D**mu_s*exp(-M2*Lam1*D/M3))

% Calculate the third moment

n=3;
sa=[5.065339 -0.062659 -3.032362 0.029469 -0.000285 0.312550 0.000204 ...
     0.003199 0.000000 -0.015952];
sb=[0.476221 -0.015896 0.165977 0.007468 -0.000141 0.060366 0.000079 ...
     0.000594 0.000000 -0.003577];

%tc=(tc.var);

log_a=sa(1)+sa(2).*tc+sa(3).*n+sa(4).*tc.*n+sa(5).*tc.^2+sa(6).*n.^2+sa(7).*tc.^2.*n+ ...
      sa(8).*tc.*n.^2+sa(9).*tc.^3+sa(10).*n.^3;
a    =10.^(log_a);
b    =sb(1)+sb(2).*tc+sb(3).*n+sb(4).*tc.*n+sb(5).*tc.^2+sb(6).*n.^2+sb(7).*tc.^2.*n+...
      sb(8).*tc.*n.^2+sb(9).*tc.^3+sb(10).*n.^3;

clear tc;
M3=a.*smo2.^b;  %big


% Z_s is also the sum of the two moments...
N0_exp=smo2.^4./(M3).^3.*Kap0;
lam_exp=smo2.*Lam0./(M3);
Zs=N0_exp.*(6.*am_s./pi./rho_w).^2./lam_exp.^(1+0+2.*bm_s) ...
    .*gamma(1+0+2.*bm_s);

N0_exp=smo2.^4./(M3).^3.*Kap1.*(smo2./(M3)).^mu_s;
lam_exp=smo2.*Lam1./(M3);
Zs=Zs+N0_exp.*(6.*am_s./pi./rho_w).^2./lam_exp.^(1+mu_s+2.*bm_s) ...
    .*gamma(1+mu_s+2.*bm_s);
Zs=0.19.*Zs./((1E-3).^6);


clear QS smo2;

% graupel
%QG=WRFUserARW(nc,'QGRAUP',time);QG=QG.var;
QG=nc{'QGRAUP'}(time,iz,ilat,ilon);
cge=[bm_g+1, mu_g+1, bm_g+mu_g+1];
cgg=gamma(cge);
obmg=1./bm_g;
oge1=1./cge(1);
ogg1=1./cgg(1);
ogg2=1./cgg(2);
% Y-intercept and slope for graupel
N0_min=gonv_max;
N0_exp=100*rho_a./(QG.*rho_a);
N0_exp=max(gonv_min,min(N0_exp,gonv_max));
N0_min=min(N0_exp,N0_min);
N0_exp=N0_min;
lam_exp=(N0_exp.*am_g*cgg(1)./(QG.*rho_a)).^oge1;
lamg=lam_exp.*(cgg(3).*ogg2.*ogg1).^obmg;
ilamg=1./lamg;
N0_g=N0_exp./(cgg(2).*lam_exp).*lamg.^cge(2);
% Z_g
Zg=0.19.*(6.*am_g./pi./rho_w).^2.*N0_g./lamg.^(1+mu_g+2.*bm_g) ...
    .*gamma(1+mu_g+2.*bm_g)./((1E-3).^6);
clear QG;

% Rain
%QR=WRFUserARW(nc,'QRAIN',time); QR=QR.var;
QR=nc{'QRAIN'}(time,iz,ilat,ilon);

cre=[bm_r+1,mu_r+1,bm_r+mu_r+1,bm_r+mu_r+1];
crg=gamma(cre);
ore1=1./cre(1);
org1=1./crg(1);
org2=1./crg(2);
obmr=1./bm_r;

N0_min=ronv_max;
N0_exp=ronv_c1.*tanh(ronv_r0.*(ronv_r0-(QR.*rho_a))) + ronv_c2;
N0_min=min(N0_exp,N0_min);
N0_exp=N0_min;
lam_exp=(N0_exp.*am_r.*crg(1)./(QR.*rho_a)).^ore1;
lamr=lam_exp.*(crg(3).*org2.*org1).^obmr;
mvd_r=(3+mu_r+0.672)./lamr;
ind=find(mvd_r>3e-3);
if(length(ind))
    mvd_r(ind)=3e-3;
	lamr(ind)=(3+mu_r+0.672)./3e-3;
	lam_exp(ind)=lamr(ind).*(crg(3).*org2.*org1).^bm_r;
	N0_exp(ind)=org1.*(QR(ind).*rho_a(ind))./am_r.*lam_exp(ind).^cre(1);
end
N0_r=N0_exp./(crg(2).*lam_exp).*lamr.^cre(2);
ilamr=1./lamr;

% Z_r
Zr=N0_r.*(6.*am_r./pi./rho_w).^2./lamr.^(1+mu_r+2.*bm_r) ...
    .*gamma(1+mu_r+2.*bm_r)./((1E-3).^6);
clear QR;

% Now calculate radar reflectivity...
% Remember moments are per m3
Zr(find(isnan(Zr(:))))=1e-30;
Zs(find(isnan(Zs(:))))=1e-30;
Zg(find(isnan(Zg(:))))=1e-30;
Z=abs(Zr)+abs(Zs)+abs(Zg);
otherwise
error('error scheme not coded');
end
