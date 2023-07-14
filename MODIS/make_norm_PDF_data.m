%create data from a normal PDF
Tau2=[0:0.1:100];
Re2=[0:0.1:30];

bins=Re2;
%   mu=12.5;
   mu=16.5;
%bins=Tau2; mu=20;

sigma=15;
sigma=3;

Ntot=1e4;

fpdf=normpdf(bins,mu,sigma);

Npoints = round(Ntot*fpdf);

iloc=1;
for i=1:length(Npoints)
    dat_norm(iloc:iloc+Npoints(i)-1)=bins(i);
    iloc=iloc+Npoints(i);
end

[N,H,W,k,Q,cw]=MODIS_N_H_func(15,dat_norm*1e-6,'calc',0);
%[N,H,W,k,Q,cw]=MODIS_N_H_func(dat_norm,12.5e-6,'calc',0);

N(N>Nd_bins(end))=Nd_bins(end-1);
N(isnan(N))=Nd_bins(end-1);
Npdf=ndHistc(squeeze(N)', Nd_bins);