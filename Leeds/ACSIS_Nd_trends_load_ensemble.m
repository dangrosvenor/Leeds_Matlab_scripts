function []=ACSIS_Nd_trends_load_ensemble()

output_period = 'all';
%output_period = 'recent';

savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm_' output_period];

ens_runs={'r1i1p1f2_u-bc179','r3i1p1f2_u-bc370','r5i1p1f3_u-az513','r7i1p1f3_u-az524','r9i1p1f2_u-bc470'...
'r2i1p1f2_u-bc292','r4i1p1f2_u-bb075','r6i1p1f3_u-az515','r8i1p1f2_u-bb277'};


load_type = 'mat';
load_type = 'merged netCDF';

var_UM = 'Nd_cf_weighted_UKESM_ztop';
var_UM = 'Nd_cf_weighted_UKESM';
var_UM = 'calipso_total_cloud_amount';


savefile = [savefile_pre_str '_' var_UM '.mat'];

for iens=1:length(ens_runs)
    switch output_period
        case 'all'
            um_case=['UKESM/' ens_runs{iens} '/output/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
        case 'recent'
            um_case=['UKESM/' ens_runs{iens} '/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
    end
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type); %data is ordered [time lat lon]. 180 times (monthly over 15 years)
    if iens==1
       Nd_ens = NaN*ones([length(ens_runs) size( dat_global.dat)]);
    end
    %store all the Nd data in a big array
    Nd_ens(iens,:,:,:) = dat_global.dat;
end

%Calculate ensemble mean and the inter-ensemble std dev
[Nd_ens_mean,Nd_ens_Ndatap Nd_ens_std] = meanNoNan(Nd_ens,1);


%% Calculate monthly means, etc.
gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;
gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);

%manually set the time since matlab time doesn't work for the 360 day
%calendar (I think this is why the times are weird).
years_ukesm = repmat([2000:2014],[12 1]);
years_ukesm = years_ukesm';
years_ukesm_1d = years_ukesm(:,1);
months_ukesm = repmat([1:12],[length(years_ukesm) 1]);
days_ukesm = ones(size(months_ukesm));
time_ukesm = datenum(years_ukesm,months_ukesm,days_ukesm);       


Nd_std_ukesm = NaN*ones([length(years_ukesm_1d) size(Nd_ens_mean,2) size(Nd_ens_mean,3)]);
clear Nd_annual Nd_ukesm
for iy=1:length(years_ukesm_1d)
    %for im=1:12
    tind_01 = (iy-1)*12 + 1;
    tind_02 = tind_01+11;
    Nd_annual(iy,:,:) = meanNoNan(Nd_ens_mean(tind_01:tind_02,:,:),1);
    Nd_ukesm(iy,:,:,:) = Nd_ens_mean(tind_01:tind_02,:,:);
    %calculate the annual means keeping each ensemble separate
    Nd_annual_ens(:,iy,:,:) = meanNoNan(Nd_ens(:,tind_01:tind_02,:,:),2);
    %end
    
%     for ilat=1:size(Nd_ens_mean,2)
%         for ilon=1:size(Nd_ens_mean,3)
%             dat_std = Nd_ens_std(tind_01:tind_02,ilat,ilon);
%             dat_mean = Nd_ens_mean(tind_01:tind_02,ilat,ilon);
%             dat_Ndatap = Nd_ens_Ndatap(tind_01:tind_02,ilat,ilon);
%             [Nd_annual_std(iy,ilat,ilon), me_out]=std_combine2(dat_std',dat_mean',dat_Ndatap');
%         end
%     end
    
end


%save(savefile,'/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat','-V7.3'...
%    ,'Nd_annual','Nd_ukesm','Nd_annual_ens','Nd_ens','Nd_ens_mean','Nd_ens_Ndatap','Nd_ens_std','gcm_Plat2D_UM','gcm_Plon2D_UM'...
%    ,'gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','years_ukesm_1d',);

save(savefile,'-V7.3'); %save all variables



