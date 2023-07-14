load('/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat','Cloud_Fraction_Liquid');
load('/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat','MLAT','MLON');
load('/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat','N_time3_37')
load('/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat','Cloud_Optical_Thickness_Liquid_Mean')
load('/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat','Cloud_Effective_Radius_37_Liquid_Mean');

ilats=find(MLAT>=20 & MLAT<=30);
ilons=find(MLON>=-140 & MLON<=-130);

tau=Cloud_Optical_Thickness_Liquid_Mean.timeseries3(ilats,ilons,:);
Nd=N_time3_37(ilats,ilons,:);
cf=Cloud_Fraction_Liquid.timeseries3(ilats,ilons,:);
re=Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(ilats,ilons,:);
A=tau./(tau+7.7);
lwp = 5/9 .* re .* tau; %or 5/9 * rhow * re*1e-6 * tau *1e3 (for g/m2).

%% AMSRE data
%see load_amsre_saved_data.m - using case as for amsre_case ='Global 2008 Aqua only';
amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2007_to_2009_20130430T150542.mat';
load(amsre_loadfile,'lwp_amsre');
lwp_amsre=1e3*flipdim(lwp_amsre,1); %convert to g/m2
load(amsre_loadfile,'month_amsre','year_amsre');
itime=find(year_amsre==2008 & month_amsre>=6 & month_amsre<=8 );

%For last index, think that 1 is daytime and 2 is nighttime
lwp_amsre2 = lwp_amsre(ilats,ilons,itime,2);


nd_bins=[0:5:350];
albedo_bins=[0.2:0.02:0.8];
albedo_bins=[0.2:0.01:0.8];

icf=find(cf>0.99);


np=ndhistc_run([Nd(icf),A(icf)],nd_bins,albedo_bins);

%Mean values for each Nd bin
qh=np;
qh(:,end+1)=0;
qh(end+1,:)=0;

pdf_norm=qh;
Ybins=nd_bins;
Xbins=albedo_bins;
mid_Ybins = 0.5 * ( Ybins(1:end-1) + Ybins(2:end) );
mid_Xbins = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );

X2d=repmat(mid_Xbins,[size(pdf_norm,1)-1 1]);
X_mean = sum(pdf_norm(1:end-1,1:end-1).*X2d,2)./sum(pdf_norm(1:end-1,1:end-1),2);

Y2d=repmat(mid_Ybins,[size(pdf_norm,2)-1 1])';
Y_mean = sum(pdf_norm(1:end-1,1:end-1).*Y2d,1)./sum(pdf_norm(1:end-1,1:end-1),1);


np_Nd=ndhistc_run([Nd(icf)],nd_bins);
np_Nd_norm = 0.2*max(X_mean) * np_Nd / max(np_Nd); %normalise and multiply by fraction of max albedo to fit on figure

figure
plot(mid_Ybins,X_mean,'ro');
xlabel('N_d (cm^{-3})');
ylabel('Bin mean cloud albedo');
hold on
plot(mid_Ybins,np_Nd_norm);
grid on
legend({'Mean albedo','Normalised frequency*0.2'});

%plot overall 2d histogram
figure
dpcolor(nd_bins,albedo_bins,np');
shading flat
xlabel('N_d (cm^{-3})');
ylabel('Cloud albedo');
colorbar;


%Restricting LWP based on MODIS - N.B. - this is NOT the way to do it since
%restricting LWP and Nd fixes tau and therefore albedo...
ilwp=find(lwp>=45 & lwp<50 & cf>0.99);
np2=ndhistc_run([Nd(ilwp),A(ilwp)],nd_bins,albedo_bins);


figure
dpcolor(nd_bins,albedo_bins,np2');
shading flat
xlabel('N_d (cm^{-3})');
ylabel('Cloud albedo');
title('CF>99% and 45<=LWP<50');
colorbar

y=1.7*nd_bins.^(1/3); yA=y./(y+7.7);
plot(nd_bins,yA,'w-');

% figure
% dpcolor(nd_bins,albedo_bins,np2');
% shading flat
% xlabel('N_d (cm^{-3})');
% ylabel('Cloud albedo');
% set(gca,'yscale','log');
% set(gca,'xscale','log');
% title('CF>99% and 45<=LWP<50');
% colorbar


%% Using AMSRE data (independent LWP)
%Find the modal LWP for the most data
dlwp=20;
lwp_bins=[0:dlwp:maxALL(lwp_amsre2*1.01)];
np_lwp=ndhistc_run([lwp_amsre2(icf)],lwp_bins);
[maxval,imax]=max(np_lwp);
lwp_mode = lwp_bins(imax);

ilwp=find(lwp_amsre2>=lwp_mode & lwp_amsre2<lwp_mode+dlwp & cf>0.99);
np2=ndhistc_run([Nd(ilwp),A(ilwp)],nd_bins,albedo_bins);

figure
dpcolor(nd_bins,albedo_bins,np2');
shading flat
xlabel('N_d (cm^{-3})');
ylabel('Cloud albedo');
title(['CF>99%, ' num2str(lwp_mode) ' <= LWP < ' num2str(lwp_mode+dlwp)]);
colorbar

%Calculate theoretical tau for known Nd and LWP values (using adiabatic
%relationships - i.e. tau propto Nd^1/3 at constant LWP)
[tau_theory]=MODIS_tau_func_N_LWP( (lwp_mode+dlwp/2)*1e-3,nd_bins*1e6); %,CTT,k,P)
%y=1.7*nd_bins.^(1/3); yA=y./(y+7.7);
yA=tau_theory./(tau_theory+7.7);
hold on
plot(nd_bins,yA,'w-');


% Do albedo means for this
qh=np2;
qh(:,end+1)=0;
qh(end+1,:)=0;

pdf_norm=qh;
Ybins=nd_bins;
Xbins=albedo_bins;
mid_Ybins = 0.5 * ( Ybins(1:end-1) + Ybins(2:end) );
mid_Xbins = 0.5 * ( Xbins(1:end-1) + Xbins(2:end) );

X2d=repmat(mid_Xbins,[size(pdf_norm,1)-1 1]);
X_mean = sum(pdf_norm(1:end-1,1:end-1).*X2d,2)./sum(pdf_norm(1:end-1,1:end-1),2);

Y2d=repmat(mid_Ybins,[size(pdf_norm,2)-1 1])';
Y_mean = sum(pdf_norm(1:end-1,1:end-1).*Y2d,1)./sum(pdf_norm(1:end-1,1:end-1),1);

figure
plot(mid_Ybins,X_mean,'ro');
xlabel('N_d (cm^{-3})');
ylabel('Mean cloud albedo in bin');
set(gca,'ylim',[0.2 0.8]);
hold on
plot(nd_bins,yA,'r-');
title(['CF>99%, ' num2str(lwp_mode) ' <= LWP < ' num2str(lwp_mode+dlwp)]);




