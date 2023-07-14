%Read Ryan_day's CTH files - can't seem to read them on the usual (older) versions of
%Matlab, so run this on e.g. r2013a


%variables (also ncdump works on Pynchon, but not on Challenger...)

% ncdump -h /home/disk/eos10/rmeast/efold/climatologies/Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_DG_hif_zb.nc
% netcdf Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_DG_hif_zb {
% dimensions:
% 	latitude = 100 ;
% 	longitude = 360 ;
% 	time = 1461 ;
% 	cseasons = 4 ;
% variables:
% 	double cth(time, longitude, latitude) ;
% 	double days(time) ;
% 	double latitude(latitude) ;
% 	double longitude(longitude) ;
% 	double months(time) ;
% 	double seasons(time) ;
% 	double years(time) ;
% 	double cth_seas(cseasons, longitude, latitude) ;
% }

file_nc = 'Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_DG_hif_zb.nc';
path_nc = ['/home/disk/eos10/rmeast/efold/climatologies/' file_nc];
file_save = ['/home/disk/eos8/d.grosvenor/CTH_Ryan/' file_nc '.mat'];

Ryan_day.cth=ncread(path_nc,'cth');
Ryan_day.days=ncread(path_nc,'days');
Ryan_day.lat=ncread(path_nc,'latitude');
Ryan_day.lon=ncread(path_nc,'longitude');
Ryan_day.months=ncread(path_nc,'months');
Ryan_day.seasons=ncread(path_nc,'seasons');
Ryan_day.years=ncread(path_nc,'years');
Ryan_day.cth_seas=ncread(path_nc,'cth_seas');

save(file_save,'Ryan_day');


%Night file
file_nc = 'Daily_1x1_JHISTO_CTH_c6_night_v2_calboxes_DG_hif_zb.ncDaily_1x1_JHISTO_CTH_c6_night_v2_calboxes_DG_hif_zb.nc';
path_nc = ['/home/disk/eos10/rmeast/efold/climatologies/' file_nc];
file_save = ['/home/disk/eos8/d.grosvenor/CTH_Ryan/' file_nc '.mat'];

Ryan_night.cth=ncread(path_nc,'cth');
Ryan_night.days=ncread(path_nc,'days');
Ryan_night.lat=ncread(path_nc,'latitude');
Ryan_night.lon=ncread(path_nc,'longitude');
Ryan_night.months=ncread(path_nc,'months');
Ryan_night.seasons=ncread(path_nc,'seasons');
Ryan_night.years=ncread(path_nc,'years');
Ryan_night.cth_seas=ncread(path_nc,'cth_seas');

save(file_save,'Ryan_night');



