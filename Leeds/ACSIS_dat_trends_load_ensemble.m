function []=ACSIS_dat_trends_load_ensemble(var_UM)

output_period = 'all';
%output_period = 'recent';
%output_period = 'SW_down_TOA_partial';

ens_runs={'r1i1p1f2_u-bc179','r3i1p1f2_u-bc370','r5i1p1f3_u-az513','r7i1p1f3_u-az524','r9i1p1f2_u-bc470'...
'r2i1p1f2_u-bc292','r4i1p1f2_u-bb075','r6i1p1f3_u-az515','r8i1p1f2_u-bb277'}; ens_str='';

%PI control - just need to calculate the time mean and std. dev (only one
%member). Runs from 1960 to 2709 giving 750 years, or 9000 months.
ens_runs={'picontrol_r1i1p1f2_u-aw310'}; ens_str='_PI_control'; output_period = 'all_PI';

savefile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_' output_period ens_str];

load_type = 'mat';
load_type = 'merged netCDF';

ind_start_offset = 0; %offset for first time index of data (since SW TOA starts with a partial year)

savefile = [savefile_pre_str '_' var_UM '.mat'];
for iens=1:length(ens_runs)
    switch output_period
        case 'all'
            um_case=['UKESM/' ens_runs{iens} '/output/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
            start_year = 1850;
            end_year = 2014;
        case 'all_PI'
            um_case=['UKESM/' ens_runs{iens} '/output/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
            start_year = 1960;
            end_year = 2709;
        case 'recent'
            um_case=['UKESM/' ens_runs{iens} '/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
            start_year = 2000;
            end_year = 2014;
        case 'SW_down_TOA_partial'
            um_case=['UKESM/' ens_runs{iens} '/output/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
            start_year = 1984;
            end_year = 2014;
            ind_start_offset = 8;
            
    end
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    opts.isort_time=1;
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type,[],[],var_UM,opts); %data is ordered [time lat lon]. 180 times (monthly over 15 years)
    if iens==1
       dat_ens = NaN*ones([length(ens_runs) size( dat_global.dat)]);
    end
    
    %deal with the ones that need masks (cloud fraction ,etc)
    switch var_UM
        case {'calipso_low_cloud_amount','calipso_mid_cloud_amount','calipso_high_cloud_amount','calipso_total_cloud_amount'}
            var_UM_mask = [var_UM '_mask'];
            dat_mask = UM_load_merged_netCDF(dirUM,var_UM_mask,run_type,load_type,[],[],var_UM_mask,opts);
            %store all the dat data in a big array
            dat_ens(iens,:,:,:) = dat_global.dat ./ dat_mask.dat;
        otherwise
            %store all the dat data in a big array
            dat_ens(iens,:,:,:) = dat_global.dat;            
    end
            
end

%Calculate ensemble mean and the inter-ensemble std dev
[dat_ens_mean,dat_ens_Ndatap dat_ens_std] = meanNoNan(dat_ens,1);


%% Calculate monthly means, etc.
gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;
%Had this the wrong way around :-
%[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);
[gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);

%manually set the time since matlab time doesn't work for the 360 day
%calendar (I think this is why the times are weird).
years_ukesm = repmat([start_year:end_year],[12 1]);
years_ukesm = years_ukesm';
years_ukesm_1d = years_ukesm(:,1);
months_ukesm = repmat([1:12],[length(years_ukesm) 1]);
days_ukesm = ones(size(months_ukesm));
time_ukesm = datenum(years_ukesm,months_ukesm,days_ukesm);       


dat_std_ukesm = NaN*ones([length(years_ukesm_1d) size(dat_ens_mean,2) size(dat_ens_mean,3)]);
clear dat_annual dat_ukesm 
for iy=1:length(years_ukesm_1d)
    %for im=1:12
    tind_01 = ind_start_offset + (iy-1)*12 + 1;
    tind_02 = tind_01+11;
    dat_annual(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
    dat_ukesm(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
    %calculate the annual means keeping each ensemble separate
    dat_annual_ens(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
    %end
    
    %DJF means
    if iy==1 %skip the first year since don't have December
        dat_annual_DJF(iy,:,:) = NaN*ones([size(dat_ens_mean,2) size(dat_ens_mean,3)]);
        %dat_ukesm_DJF(iy,:,:,:) = NaN;
        dat_annual_ens_DJF(:,iy,:,:) = NaN*ones([size(dat_ens,1) size(dat_ens,3) size(dat_ens,4)]);
    else
        tind_01 = ind_start_offset + (iy-1)*12 + 0; %+1 is Jan, so +0 is Dec
        tind_02 = tind_01+2; %Feb
        dat_annual_DJF(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
        %dat_ukesm_DJF(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
        %calculate the annual means keeping each ensemble separate
        dat_annual_ens_DJF(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
    end
    
    % MAM
    tind_01 = ind_start_offset + (iy-1)*12 + 3; %+1 is Jan
    tind_02 = tind_01+2; %May
    dat_annual_MAM(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
    %dat_ukesm_MAM(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
    %calculate the annual means keeping each ensemble separate
    dat_annual_ens_MAM(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
    
    % JJA
    tind_01 = ind_start_offset + (iy-1)*12 + 6; %+1 is Jan
    tind_02 = tind_01+2; %Aug
    dat_annual_JJA(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
    %dat_ukesm_JJA(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
    %calculate the annual means keeping each ensemble separate
    dat_annual_ens_JJA(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);
    
     % SON
    tind_01 = ind_start_offset + (iy-1)*12 + 9; %+1 is Jan
    tind_02 = tind_01+2; %Nov
    dat_annual_SON(iy,:,:) = meanNoNan(dat_ens_mean(tind_01:tind_02,:,:),1);
    %dat_ukesm_SON(iy,:,:,:) = dat_ens_mean(tind_01:tind_02,:,:);
    %calculate the annual means keeping each ensemble separate
    dat_annual_ens_SON(:,iy,:,:) = meanNoNan(dat_ens(:,tind_01:tind_02,:,:),2);   
    
end


%save(savefile,'/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat','-V7.3'...
%    ,'Nd_annual','Nd_ukesm','Nd_annual_ens','Nd_ens','Nd_ens_mean','Nd_ens_Ndatap','Nd_ens_std','gcm_Plat2D_UM','gcm_Plon2D_UM'...
%    ,'gcm_Plat2D_edges_UM','gcm_Plon2D_edges_UM','years_ukesm_1d',);

save(savefile,'-V7.3'); %save all variables


