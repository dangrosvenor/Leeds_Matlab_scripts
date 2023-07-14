


idomain = find( ~( Plat_global>=thresh_LAT(1) & Plat_global<thresh_LAT(2) & Plon_global>=thresh_LON(1) & Plon_global<thresh_LON(2) ) );

idays = find( ~(Ndays2>=thresh_ndays) );



dat_modis=N_time3;
dat_modis(ihtot)=NaN;
dat_modis=permute(dat_modis,[3 1 2]);
%dat_modis=dat_modis(:,:);
dat_modis2=dat_modis(:,:);

dat_modis2(:,idomain)=NaN;
dat_modis2(:,idays)=NaN;

Pall = dat_modis2(time_inds_average,:);
Pall = Pall(:);
Pall2 = Pall;
Pall2(isnan(Pall))=[];
Nd_all3 = reshape(dat_modis2,size(dat_modis));
Nd_all3 = Nd_all3(time_inds_average,:,:);
[Nd_mean,Nnums,Nd_std] = meanNoNan(Nd_all3,1);
me3D = repmat(Nd_mean,[1 1 size(Nd_all3,1)]); me3D = permute(me3D,[3 1 2]);
std3D = repmat(Nd_std,[1 1 size(Nd_all3,1)]); std3D = permute(std3D,[3 1 2]);



dat_modis=Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
dat_modis(ihtot)=NaN;
dat_modis=permute(dat_modis,[3 1 2]);
dat_modis2=dat_modis;

dat_modis2(:,idomain)=NaN;
dat_modis2(:,idays)=NaN;

Tau_all = dat_modis2(time_inds_average,:);
Tau_all = Tau_all(:);
Tau_all2 = Tau_all;
Tau_all2(isnan(Tau_all))=[];
Tau_all3 = reshape(dat_modis2,size(dat_modis));
Tau_all3 = Tau_all3(time_inds_average,:,:);
[Tau_mean,Nnums,Tau_std] = meanNoNan(Tau_all3,1);
meTau3D = repmat(Tau_mean,[1 1 size(Nd_all3,1)]); meTau3D = permute(meTau3D,[3 1 2]);


dat_modis=Cloud_Effective_Radius_Liquid_Mean.timeseries3;
dat_modis(ihtot)=NaN;
dat_modis=permute(dat_modis,[3 1 2]);
dat_modis2=dat_modis;

dat_modis2(:,idomain)=NaN;
dat_modis2(:,idays)=NaN;

Re_all = dat_modis2(time_inds_average,:);
Re_all = Re_all(:);
Re_all2 = Re_all;
Re_all2(isnan(Re_all))=[];
Re_all3 = reshape(dat_modis2,size(dat_modis));
Re_all3 = Re_all3(time_inds_average,:,:);
[Re_mean,Nnums,Re_std] = meanNoNan(Re_all3,1);
meRe3D = repmat(Re_mean,[1 1 size(Nd_all3,1)]); meRe3D = permute(meRe3D,[3 1 2]);


% [P,Npoints] = meanNoNan(dat_modis2(time_inds_average,:),1); %tim
% 
% [Pmean,NPmean] = meanNoNan(P(:),1); %tim




mean_P = mean(Pall2);
std_P = std(Pall2);

mean_Tau = mean(Tau_all2);
mean_Re = mean(Re_all2);

istd=find(Pall>mean_P+2*std_P);


iNd = find( Nd_all3 - me3D > 2*std3D );
dTau = Tau_all3(iNd) ./ meTau3D(iNd);
dN_Tau = sqrt(dTau);
dRe = Re_all3(iNd) ./ meRe3D(iNd);
dN_Re = dRe.^(-5/2);

dN_Tau2 = Nd_all3(iNd) .* 0.5 .* (Tau_all3(iNd)-meTau3D(iNd))./meTau3D(iNd);
dN_Re2 = Nd_all3(iNd) .* 2.5 .* (Re_all3(iNd)-meRe3D(iNd))./meRe3D(iNd);




