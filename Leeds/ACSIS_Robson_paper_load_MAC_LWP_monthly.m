function ACSIS_Robson_paper_load_MAC_LWP_monthly()

%Load in monthly cloud fraction anomaly data from Norris Nature 2016 paper

file_dir = '/home/disk/eos5/d.grosvenor/MAC_LWP/monthly_MAC/'

lwp_type = 'lwp'; %LWP only
lwp_type = 'tlwp'; %LWP+RWP only

switch lwp_type
    case 'lwp'
        save_file = [file_dir 'MAC_monthly_LWP.mat'];
        files_nc = dir([file_dir 'maclwp_cloudlwpave*.nc']); %diurnal avearge LWP (no rain)
        var_name = 'cloudlwp';
    case 'tlwp'
        save_file = [file_dir 'MAC_monthly_TLWP.mat'];
        files_nc = dir([file_dir 'maclwp_totallwpave*.nc']); %diurnal avearge LWP (no rain)
        var_name = 'totallwp';
end



%% 

for ifile=1:length(files_nc)
    
nc  = netcdf([file_dir files_nc(ifile).name]);

years(ifile) = str2num(files_nc(ifile).name(20:23));
%time(ifile) = nc{'time'}(:); %"Month index (0-11); 0=January, ..., 11=December"
%time_matlab_in(ifile) = datenum('01-Jan-1970') + time;

if ifile==1
    lat= nc{'lat'}(:); dlat=abs(mean(diff(lat)));
    lon = nc{'lon'}(:); dlon=abs(mean(diff(lon(1:2)))); %0 to 360 convention whereas we want -180 to +180.
    [gcm_Plon2D_360,gcm_Plat2D_360] = meshgrid(lon,lat);
    i180 = find(lon>180);
    lon(i180) = lon(i180)-360;
    [gcm_Plon2D,gcm_Plat2D] = meshgrid(lon,lat);
end

lwp_load = nc{var_name}(:); % g/m2 size = [12 180 360]
lwp_err_load = nc{[var_name '_error']}(:); % g/m2 (1-sigma error)

nc=close(nc);

lwp(ifile,:,:,:) = screen_fill(lwp_load);
lwp_err(ifile,:,:,:) = screen_fill(lwp_err_load);

end


nyears = length(years);
for iyear=1:nyears
    %istart = (iyear-1)*12+1;
    [lwp_annual(iyear,:,:) N stdddev] = meanNoNan(lwp(iyear,:,:,:),2);
    sum_sq = meanNoNan(lwp_err(iyear,:,:,:).^2,2,'sum');
    lwp_err_annual(iyear,:,:) = 1./N .* sqrt(sum_sq);
        %combinign the errors 1/n * sqrt(sum(sig_i^2))
end

%% Re-grid to UM grid
UM_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
load(UM_file,'gcm_Plon2D_UM','gcm_Plat2D_UM');
i180 = find(gcm_Plon2D_UM<0);
gcm_Plon2D_UM_360 = gcm_Plon2D_UM;
gcm_Plon2D_UM_360(i180) = gcm_Plon2D_UM_360(i180) + 360;


nlat = size(gcm_Plon2D_UM,1);
nlon = size(gcm_Plon2D_UM,2);
lwp_UM_monthly = NaN * ones([nyears 12 nlat nlon]);
lwp_UM_annual = NaN * ones([nyears nlat nlon]);
lwp_err_UM_monthly = NaN * ones([nyears 12 nlat nlon]);
lwp_err_UM_annual = NaN * ones([nyears nlat nlon]);

for iyear=1:nyears   
    dat = squeeze(lwp_annual(iyear,:,:));            
    %lwp_UM_annual(itime,:,:) = griddata(gcm_Plat2D_360,gcm_Plon2D_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);    
    lwp_UM_annual(iyear,:,:) = griddata(gcm_Plat2D,gcm_Plon2D,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);  
    
    dat = squeeze(lwp_err_annual(iyear,:,:));
    lwp_err_UM_annual(iyear,:,:) = griddata(gcm_Plat2D,gcm_Plon2D,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);        
end


for iyear=1:nyears
    for im=1:12
        dat = squeeze(lwp(iyear,im,:,:));
        %lwp_UM(itime,:,:) = griddata(gcm_Plat2D_360,gcm_Plon2D_360,dat,gcm_Plat2D_UM,gcm_Plon2D_UM_360);
        lwp_UM_monthly(iyear,im,:,:) = griddata(gcm_Plat2D,gcm_Plon2D,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        dat = squeeze(lwp_err(iyear,im,:,:));        
        lwp_err_UM_monthly(iyear,im,:,:) = griddata(gcm_Plat2D,gcm_Plon2D,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);
    end
end



save(save_file,'-V7.3');

function [dat_out] = screen_fill(dat)

inan = find(dat<-998);
dat(inan) = NaN;
dat_out = dat;