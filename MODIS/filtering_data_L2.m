lat = inpaint_nans(lat); %extrapolates to fill in NaN values in case of missing rows, etc.
lon = inpaint_nans(lon); %extrapolates to fill in NaN values (inpaint_nans is located in peoples scripts folder)

Plat_L2 = split_matrix_int(lat,5);
Plon_L2 = split_matrix_int(lon,5);
Plat_L2_5km = lat;
Plon_L2_5km = lon;
cfL2 = split_matrix(cf,5);

%not sure about this indexing.... assuming lat=lat_full(4:5:1353,4:5:2029)
%size(tau)=[1354 2030], so ignoring the last point and starting from index
%4. But should just be for plotting, so should be ok
%row_L2_inds = [4:1349];
%col_L2_inds = [4:2029];
%row_L2_inds = [4:size(Plon_L2,1)+3];
%col_L2_inds = [4:size(Plon_L2,2)+3];

%think that this should be the correct sampling - as suggested from the
%variable attributes
row_L2_inds = [sampling_across(1):sampling_across(2)];
%here will be ignoring the pixels at the beginning and end of each swath as they don't
%have a lat,lon, etc (5km variables) values and so we would need to
%extrapolate
col_L2_inds = [sampling_along(1):sampling_along(2)];

%these are the indices of the points for all of the pixels in each 5x5km
%square of Plat_L2
row_L2_inds2 = [sampling_across(1)-2:sampling_across(2)+2];
col_L2_inds2 = [sampling_along(1)-2:sampling_along(2)+2];



t_top2 = NaN*ones(size(re_un));
t_top2(row_L2_inds,col_L2_inds) = split_matrix_int(t_top,5);

p_top2 = NaN*ones(size(re_un));
p_top2(row_L2_inds,col_L2_inds) = split_matrix_int(p_top,5);

solarZA_L2 = NaN*ones(size(re_un));
solarZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(solar_zenith,5);

solarAz_L2 = NaN*ones(size(re_un));
solarAz_L2(row_L2_inds,col_L2_inds) = split_matrix_int(solar_azimuth,5);

sensorZA_L2 = NaN*ones(size(re_un));
sensorZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(sensor_zenith,5);

sensorAz_L2 = NaN*ones(size(re_un));
sensorAz_L2(row_L2_inds,col_L2_inds) = split_matrix_int(sensor_azimuth,5);

%Relative azimuth is defined as zero when the sensor is looking
            %into the Sun - as if there was no scattering (just
%             %transmission). Backscatter is then at 180 degrees. So do the
%             difference of sensor and solar and subtract 180 so that if
%             the difference is 180 (forward scatter) then will get relAZ=0
%NOTE - %with Joint files only get the relative AZ, but have kept the variable
%Solar_Azimuth_Angle for consistency with full L2 and have set the sensorAZ to be the
%relative AZ. But are the MODIS Relative_Azimuth_Angles with the 180
%already removed? I.e. do they follow the above convention? Check using L2?
%Yes, if calculate relAz as below using L2 then get the same answer as that
%given by Relative_Azimuth from the jointL2 product.
relAz_L2 = calc_relAz(sensorAz_L2,solarAz_L2);

%extrapolate the 5km fields to fill the whole swath
[Plat2_L2]=expand_swath_lin_extrap(Plat_L2,sampling_along,sampling_across,re);
[Plon2_L2]=expand_swath_lin_extrap(Plon_L2,sampling_along,sampling_across,re);
[t_top3]=expand_swath_lin_extrap(t_top2(row_L2_inds,col_L2_inds),sampling_along,sampling_across,re);

%go one to the side for the edges
[Plat3_L2]=Plat2_L2(sampling_across(1)-1:sampling_across(2),sampling_along(1)-1:sampling_along(2));
[Plon3_L2]=Plon2_L2(sampling_across(1)-1:sampling_across(2),sampling_along(1)-1:sampling_along(2));

Wflag='calc';
[N,H,W,k,Q,cw]=MODIS_N_H_func(tau,re*1e-6,Wflag,NaN,t_top2);

if exist('SST')==1 %==1 tests if it is a variable, since this also returns true if it's a directory, file, etc.
    [sst,dates,diff_days] = sst_for_date_range_McCoy_func(SST,sst_time,modis_date_time,'nearest');
    sst_L2 = GRIDDATA(lon_sst,lat_sst,(squeeze(sst)),Plon_L2_5km,Plat_L2_5km);
    %            ctt_L2 = GRIDDATA(Plon_L2_5km,Plat_L2_5km,(squeeze(t_top)),Plon_L2,Plat_L2);
    CTH_L2_5km = (273.15 + sst_L2 - t_top - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km        
    CTH_L2 = GRIDDATA(Plon_L2_5km,Plat_L2_5km,CTH_L2_5km,Plon2_L2,Plat2_L2);
end

%re_un and tau_un are given in % error (checked in the hdf attributes)
 percent_error_Nd = sqrt( (0.5*tau_un).^2 + (-5/2*re_un).^2 );
 re_un_abs = re_un.*re/100;
 tau_un_abs = tau_un.*tau/100;


%phase_flag = squeeze(qapq_1km(8,:,:));   %1621 phase flag
phase_flag = squeeze(qapp_1km(1,:,:));    %"Primary" phase flag
phase_retreival_outcome = squeeze(qapp_1km(2,:,:)); %Outcome of primary phase retrieval above. 0=not attempted
%(possibly clear, although can actually have a phase for these - but it
%should be ignored). 1=successful, is cloudy. All clear points have =0
%(&phase_flag==1)
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
