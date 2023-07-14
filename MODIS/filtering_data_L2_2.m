Plat2 = split_matrix(lat,5);
Plon2 = split_matrix(lon,5);
cfL2 = split_matrix(cf,5);

%lat grid samples the centre of every 5x5 pixel block
%so know that pixels 1-5 across in e.g. re array are correspond to the
%first lat value
%make arrays of lat and lon for each 1 km pixel

row_L2_inds = [sampling_across(1)-2:sampling_across(2)+2];
%here will be ignoring the pixels at the end of each swath as they don't
%have a lat,lon geolocation value
col_L2_inds = [sampling_along(1)-2:sampling_along(2)+2];


Plat = Plat2(row_L2_inds,col_L2_inds);
Plon = Plon2(row_L2_inds,col_L2_inds);


%not sure about this indexing.... assuming lat=lat_full(4:5:1353,4:5:2029)
%size(tau)=[1354 2030], so ignoring the last point and starting from index
%4. But should just be for plotting, so should be ok
%row_L2_inds = [4:1349];
%col_L2_inds = [4:2029];

%row_L2_inds = [4:size(Plon,1)+3];
%col_L2_inds = [4:size(Plon,2)+3];

row_L2_inds = [sampling_across(1):sampling_across(2)];
%here will be ignoring the pixels at the beginning and end of each swath as they don't
%they are outside of the values for the 5km grids - would need to
%extrapolat
col_L2_inds = [sampling_along(1):sampling_along(2)];


t_top2 = NaN*ones(size(re_un));
t_top2(row_L2_inds,col_L2_inds) = split_matrix_int(t_top,5);

p_top2 = NaN*ones(size(re_un));
p_top2(row_L2_inds,col_L2_inds) = split_matrix_int(p_top,5);

solarZA_L2 = NaN*ones(size(re_un));
solarZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(solar_zenith,5);

sensorZA_L2 = NaN*ones(size(re_un));
sensorZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(sensor_zenith,5);


Wflag='calc';
[N,H,W,k,Q,cw]=MODIS_N_H_func(tau,re*1e-6,Wflag,NaN,t_top2);

%re_un and tau_un are given in % error (checked in the hdf attributes)
 percent_error_Nd = sqrt( (0.5*tau_un).^2 + (-5/2*re_un).^2 );
 re_un_abs = re_un.*re/100;
 tau_un_abs = tau_un.*tau/100;


%phase_flag = squeeze(qapq_1km(8,:,:));   %1621 phase flag
phase_flag = squeeze(qapp_1km(1,:,:));    %"Primary" phase flag
surface_flag = squeeze(mask_1km(6,:,:));
tau_bounds = squeeze(qapq_1km(3,:,:)); %tau out of bounds flag


% cf = Cloud_Fraction_Liquid.data;
% WMOD=(5/6)*Cloud_Water_Path_Liquid_Mean.data/1000; %convert to kg/m2
% tau = Cloud_Optical_Thickness_Liquid_Mean.data;
% %                reff = Cloud_Effective_Radius_Liquid_Mean.data; %convert to metres
% 
% try
% 
% cf_time=Cloud_Fraction_Liquid.timeseries(:,:);
% NP_time=Cloud_Fraction_Liquid_Pixel_Counts.timeseries(:,:)./cf_time;
% 
% tau_time = Cloud_Optical_Thickness_Liquid_Mean.timeseries(:,:);
% reff_time = Cloud_Effective_Radius_Liquid_Mean.timeseries(:,:)*1e-6; %convert to metres
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% 
% [N_time,H,W,k,Q,cw]=MODIS_N_H_func(tau_time,reff_time,Wflag,0);
% 
% 
% %timeseries3
% cf_time3=Cloud_Fraction_Liquid.timeseries3;
% NP_time3=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3./cf_time3;
% SZA_time3=Solar_Zenith_Mean.timeseries3;
% sensZA_time3 = Sensor_Zenith_Mean.timeseries3;
% 
% tau_time3 = Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
% reff_time3 = Cloud_Effective_Radius_Liquid_Mean.timeseries3*1e-6; %convert to metres
% Wflag='calc'; %calculate LWP using the Eq. 6 in Bennartz (2007)
% 
% [N_time3,H,W,k,Q,cw]=MODIS_N_H_func(tau_time3,reff_time3,Wflag,0);
% 
% catch
% end
