function SW_DEEPC_data_read_FUNC()
%Takes the net downward SW at TOA from the Allan dataset and converts it to
%the upwelling SW TOA flux using incoming SW at TOA from the UM (r1
%ensemble member)
icoarse_grain=0;

dat_dir = '/home/disk/eos15/d.grosvenor/eos8/DEEPC_CERES/';
save_file = [dat_dir 'sw_up_TOA_calculated.mat'];
%save_file = [dat_dir 'sw_up_TOA_calculated_non_coarse_grained_non_GEOD.mat'];
%save_file = [dat_dir 'sw_up_TOA_calculated_non_coarse_grained.mat'];

%Net downward SW rad flux at TOA.
sw_down_net_file = 'DEEPC_TOA_ASR_v03.0_198501-201601_GEOD.nc';
%sw_down_net_file = 'DEEPC_TOA_ASR_v03.0_198501-201601.nc';
nc=netcdf([dat_dir sw_down_net_file]);
lat_deepc = nc{'lat'}(:); dlat=abs(mean(diff(lat_deepc)));
lon_deepc = nc{'lon'}(:); dlon=abs(mean(diff(lon_deepc(1:2)))); %0-360 convection
lon_deepc(lon_deepc>180)=lon_deepc(lon_deepc>180)-360;
sw_deepc = nc{'ASR'}(1:360,:,:); %use data to end of 2014.
time = nc{'time'}(1:360); %days since 1985-1-1; starts Jan 1985, ends Jan 2016
time_matlab = datenum('01-Jan-1985') + time;

[lon2d_deepc,lat2d_deepc] = meshgrid(lon_deepc,lat_deepc);

%Load in the UKESM downwelling SW TOA data.
load_type = 'merged netCDF';
var_UM = 'SW_down_TOA';
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/UKESM/r1i1p1f2_u-bc179/output/']; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa';
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

 %size dat_global.dat = [380         144         192];
 %380 is 8 months of 1983 and all 12 months for 1984-2014
 %So ignore the 1st 8 months + 1 year
dlat_UM = abs(mean(diff(dat_global.gcm_Plat2D_UM(:,1))));
dlon_UM = abs(mean(diff(dat_global.gcm_Plon2D_UM(1,1:2))));
yr_start_UM = 1985;
yr_end_UM = 2014;
sw_in_TOA_UM = dat_global.dat(21:end,:,:);  %1985-2014


%coarse-grain the DEEP-C data to approximate the resolution of the UM data
if icoarse_grain==1
    M_coarse_grain=round(dlat_UM/dlat);
    N_coarse_grain=round(dlon_UM/dlon);
else
    M_coarse_grain=1;
    N_coarse_grain=1;
end

Plat2D_coarse = reduce_matrix_subsample_mean(lat2d_deepc,M_coarse_grain,N_coarse_grain);
Plon2D_coarse = reduce_matrix_subsample_mean(lon2d_deepc,M_coarse_grain,N_coarse_grain);
%gcm_Plon2D_GENERIC = reduce_matrix_subsample_mean(gcm_Plon2D_UM,M_coarse_grain,N_coarse_grain);
%[gcm_Plon2D_edges_GENERIC,gcm_Plat2D_edges_GENERIC] = get_edges_lat_lon(gcm_Plon2D_GENERIC,gcm_Plat2D_GENERIC);
clear dat_coarse sw_net_coarse
for it=1:length(time)
    dat_coarse(it,:,:) = reduce_matrix_subsample_mean(squeeze(sw_deepc(it,:,:)),M_coarse_grain,N_coarse_grain);
    sw_net_coarse(it,:,:) = griddata(Plat2D_coarse,Plon2D_coarse,squeeze(dat_coarse(it,:,:)),dat_global.gcm_Plat2D_UM,dat_global.gcm_Plon2D_UM);
end

%sw_down_net = sw_down_TOA - SW_up_TOA
sw_up_toa = sw_in_TOA_UM - sw_net_coarse;

save(save_file,'-V7.3','sw_up_toa','sw_net_coarse','M_coarse_grain','N_coarse_grain','yr_start_UM','yr_end_UM',...
    'sw_in_TOA_UM','sw_deepc','time_matlab','lat2d_deepc','lon2d_deepc');




