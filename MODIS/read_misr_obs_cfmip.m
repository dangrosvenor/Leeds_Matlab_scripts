savemem=1;  %flag so that certain fields are not read in to save memory (mainly 3D fields)

%nc_dir = '/home/disk/margaret/hillmanb/nobackup/climate_observations/calipso/cmor/monthly/';
nc_dir = '/home/disk/margaret/hillmanb/nobackup/climate_observations/misr/cmor/climo/';
%nc_dir2 = '/home/disk/eos8/d.grosvenor/CPT/CAM5/';
savedir = '/home/disk/eos1/d.grosvenor/modis_work/plots/';


comp='UWchallenger';


nc_inst_file = ['clMISR_ANN.nc']; gcm_str='MISR_CMOR_ANN';



%nc_grid=netcdf([nc_dir nc_grid_file],'nowrite');
nc_inst=netcdf([nc_dir nc_inst_file],'nowrite');

%gcm_phalf_ref=nc_grid{'phalf'}(:);
gcm_lat_full = nc_inst{'lat'}(:);
gcm_lon_full = nc_inst{'lon'}(:);

gcm_slat_full = -90:1:90;
gcm_slon_full = 0:1:360;

%gcm_lat_full([1 end]) ---> ans = [-90 90]  &  size(gcm_lat_full) = [96 1], i.e. cell edges
%gcm_slat_full([1 end]) ---> ans= [-89.0526 89.0526]  & size(gcm_lat_full)= [95 1], i.e. cell centers

%gcm_lon_full([1 end]) --> ans = [0 357.500]  &  size(gcm_lon_full) = [144 1], likely cell edges
%gcm_slon_full([1 end]) --> ans = [-1.2500 356.2500]  &  size(gcm_lon_full) = [144 1], likely cell centers

lat_range = [-40 10];
lon_range = [-140 -50]+360;

ilat=find(gcm_slat_full>=lat_range(1) & gcm_slat_full<lat_range(2));
%go one either side
%ilat = [max([1 ilat(1)-1]); ilat; min([length(gcm_lat_full) ilat(end)+1])];
%ilat = [max([1 ilat(1)-1]):min([length(gcm_lat_full) ilat(end)+1])];

%ilat = [ilat(1):ilat(end)];

ilon=find(gcm_slon_full>=lon_range(1) & gcm_slon_full<lon_range(2));
%go one either side
%ilon = [max([1 ilon(1)-1]):min([length(gcm_lon_full) ilon(end)+1])];


%ilon = [ilon(1):ilon(end)];

gcm_lat = gcm_lat_full(ilat);
gcm_lon = gcm_lon_full(ilon);
%gcm_slat = gcm_slat_full(ilat);
%gcm_slon = gcm_slon_full([ilon ilon(end)+1]);

gcm_slat = gcm_slat_full([ilat(1):ilat(end)+1]);
gcm_slon = gcm_slon_full([ilon(1):ilon(end)+1]);

%ilat=[ilat(2:end); ilat(end)+1];
%ilon=[ilon; ilon(end)+1];

%have monthly output for COSP - each file is one month
gcm_idays=[1];

% gcm_time_read = nc_inst{'time'}(gcm_idays); %"days since 2000-01-01 00:00:00"
% gcm_time_matlab = gcm_time_read + datenum('01-Jan-0000'); %convert to Matlab time
% [Y,MO,D,H,MI,S] = datevec(gcm_time_matlab); %this outputs numbers for the date components
% 
% gcm_time_days = D;
% 
% %will be consistent with the actual calendar day rather than nearest day
% 
% %the hour of the day
% gcm_time_UTC = H;
% %the month
% gcm_month = MO;
% 
% gcm_decimal_days = gcm_time_read - gcm_time_read(1);
% %days of year since beginning of the year

% daynum_timeseries3 = floor(gcm_time_matlab - datenum(Y(1),1,1) + 1);
daynum_timeseries3 = 1; 
 eval(['daynum_timeseries3_' gcm_str ' = daynum_timeseries3;']);
%zero difference means day 1




%clMISR(time, tau, cth, lat, lon) ;
%clMISR:_FillValue = 1.e+20f ;
%clMISR:long_name = "Cloud Fraction as Caculated by the MISR Simulator" (%)

CF_vs_tau_Z_misr = read_calipso_cmor(nc_dir,nc_inst_file,'clMISR',ilat,ilon,'(:,:,:,ilat,ilon)'); 
%Matlab needs edge arrays to be equal sized
%to data array for some reason! It ignores the last indices of the array
%for plotting. So an option is to add a NaN to the end of each index - this
%is done in read_calipso_cmor

cth_misr = nc_inst{'cth'}(:);
tau_misr = nc_inst{'tau'}(:);
%NOTE = "First entry value of -1 is for NO RETRIEVAL cases where cloud is detected but
%cth and/or tau retrieval fails" ;

%think the above are bin-centers
%tau bin edges are as for ISCCP I think
tau_misr_edges = [0 0.3 1.3 3.6 9.4 23 60 379];
cth_misr_edges = 1000*[0:0.5:3 4 5 7 9 11 13 15 17 23];

Npix_misr = read_calipso_cmor(nc_dir,nc_inst_file,'MISR_total_samples',ilat,ilon,'(ilat,ilon)'); 



% -----------------------------------
%  lat lon and time indices
% -----------------------------------
%cell edges (contained in gcm_slat)
% Plat=gcm_slat;
% Plon=gcm_slon;
% 
% Plat=gcm_lat_full;
% Plon=gcm_lon_full;
% 
% Plat=gcm_slat_full;
% Plon=gcm_slon_full;

Plat=gcm_slat;
Plon=gcm_slon;


i180=find(Plon>180);
Plon2=Plon;
Plon2(i180)=Plon2(i180)-360;

%Plon=[Plon2(1:end); Plon2(1)];
Plon=Plon2;



[gcm_Plon2D_edges,gcm_Plat2D_edges]=meshgrid(Plon,Plat);


%cell centres - contained in gcm_lat - seems that the fields are stored as point values at gcm_lat, which runs to lat=+/-90
% - so gcm_lat are effectively the cell face positions. And we ignore the
% first and last lat value (since cells can't extend beyond 90 degree lat)
Plat=gcm_lat;
Plon=gcm_lon;

i180=find(Plon>180);
Plon(i180)=Plon(i180)-360;

[gcm_Plon2D,gcm_Plat2D]=meshgrid(Plon,Plat);

%dlon = 1.87, dlat = 2.5. Both constant throughout the grid for CAM5

% dlat = diff(gcm_lat);
% dlon = diff(gcm_lon);
% Plat=gcm_lat+[dlat; dlat(end)]/2;
% Plon=gcm_lon+[dlon; dlon(end)]/2;
% 
% i180=find(Plon>180);
% Plon(i180)=Plon(i180)-360;
% 
% [gcm_Plon2D_edges,gcm_Plat2D_edges]=meshgrid(Plon,Plat);










