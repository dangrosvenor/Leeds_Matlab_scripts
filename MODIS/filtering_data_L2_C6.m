%Plat2_L2 is the lat array that matches the re21, etc. arrays
%Made using [Plat2_L2]=expand_swath_lin_extrap(Plat_L2,sampling_along,sampling_across,re21);
%Plat_L2 made using :-
%Plat_L2 = split_matrix_int(lat,5);

lat = inpaint_nans(lat); %extrapolates to fill in NaN values in case of missing rows, etc.
lon = inpaint_nans(lon); %extrapolates to fill in NaN values (inpaint_nans is located in peoples scripts folder)

Plat_L2 = split_matrix_int(lat,5);
Plon_L2 = split_matrix_int(lon,5);
Plat_L2_5km = lat;
Plon_L2_5km = lon;
cfL2 = split_matrix(cf,5);

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
% The above are not used in this script - are they used anywhere??

% Don't need t_top2 and p_top2 now we have 1km ctt, etc.
% t_top2 = NaN*ones(size(re21_un));
% t_top2(row_L2_inds,col_L2_inds) = split_matrix_int(t_top,5);
% 
% p_top2 = NaN*ones(size(re21_un));
% p_top2(row_L2_inds,col_L2_inds) = split_matrix_int(p_top,5);

solarZA_L2 = NaN*ones(size(re21_un));
solarZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(solar_zenith,5);

solarAz_L2 = NaN*ones(size(re21_un));
solarAz_L2(row_L2_inds,col_L2_inds) = split_matrix_int(solar_azimuth,5);

sensorZA_L2 = NaN*ones(size(re21_un));
sensorZA_L2(row_L2_inds,col_L2_inds) = split_matrix_int(sensor_zenith,5);

sensorAz_L2 = NaN*ones(size(re21_un));
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

%extrapolate the 5km fields to fill the whole swath - to match e.g. re21
%array
[Plat2_L2]=expand_swath_lin_extrap(Plat_L2,sampling_along,sampling_across,re21);
[Plon2_L2]=expand_swath_lin_extrap(Plon_L2,sampling_along,sampling_across,re21);
%[t_top3]=expand_swath_lin_extrap(t_top2(row_L2_inds,col_L2_inds),sampling_along,sampling_across,re);

%go one to the side for the edges
[Plat3_L2]=Plat2_L2(sampling_across(1)-1:sampling_across(2),sampling_along(1)-1:sampling_along(2));
[Plon3_L2]=Plon2_L2(sampling_across(1)-1:sampling_across(2),sampling_along(1)-1:sampling_along(2));

%Create new tau arrays for Nd calc in order to only use tau>5
tau16_Nd = tau16; tau16_Nd(tau16<5) = NaN;
tau21_Nd = tau21; tau21_Nd(tau21<5) = NaN;
tau37_Nd = tau37; tau37_Nd(tau37<5) = NaN;

Wflag='calc';
[N16,H16,W16,k,Q,cw]=MODIS_N_H_func(tau16,re16*1e-6,Wflag,NaN,t_top_1km);
[N21,H21,W21,k,Q,cw]=MODIS_N_H_func(tau21,re21*1e-6,Wflag,NaN,t_top_1km);
[N37,H37,W37,k,Q,cw]=MODIS_N_H_func(tau37,re37*1e-6,Wflag,NaN,t_top_1km);

[N16_pcl,H16_pcl,W16_pcl,k,Q,cw]=MODIS_N_H_func(tau16_pcl,re16_pcl*1e-6,Wflag,NaN,t_top_1km);
[N21_pcl,H21_pcl,W21_pcl,k,Q,cw]=MODIS_N_H_func(tau21_pcl,re21_pcl*1e-6,Wflag,NaN,t_top_1km);
[N37_pcl,H37_pcl,W37_pcl,k,Q,cw]=MODIS_N_H_func(tau37_pcl,re37_pcl*1e-6,Wflag,NaN,t_top_1km);

[N16_Nd,H16_Nd,W16_Nd,k,Q,cw]=MODIS_N_H_func(tau16_Nd,re16*1e-6,Wflag,NaN,t_top_1km);
[N21_Nd,H21_Nd,W21_Nd,k,Q,cw]=MODIS_N_H_func(tau21_Nd,re21*1e-6,Wflag,NaN,t_top_1km);
[N37_Nd,H37_Nd,W37_Nd,k,Q,cw]=MODIS_N_H_func(tau37_Nd,re37*1e-6,Wflag,NaN,t_top_1km);

N16(ind_not_liq)=NaN; %filter for non-liquid pixels
H16(ind_not_liq)=NaN; %filter for non-liquid pixels
W16(ind_not_liq)=NaN; %filter for non-liquid pixels
N21(ind_not_liq)=NaN; %filter for non-liquid pixels
H21(ind_not_liq)=NaN; %filter for non-liquid pixels
W21(ind_not_liq)=NaN; %filter for non-liquid pixels
N37(ind_not_liq)=NaN; %filter for non-liquid pixels
H37(ind_not_liq)=NaN; %filter for non-liquid pixels
W37(ind_not_liq)=NaN; %filter for non-liquid pixels

N16_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
H16_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
W16_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
N21_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
H21_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
W21_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
N37_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
H37_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels
W37_pcl(ind_not_liq)=NaN; %filter for non-liquid pixels

%More conservative Nd where also filter for sunglint, heavy aerosol, thin
%cirrus and shadow_flag
N16_cons = N16; N16_cons(ind_dubious) = NaN;
H16_cons = H16; H16_cons(ind_dubious) = NaN;
W16_cons = W16; W16_cons(ind_dubious) = NaN;
N21_cons = N21; N21_cons(ind_dubious) = NaN;
H21_cons = H21; H21_cons(ind_dubious) = NaN;
W21_cons = W21; W21_cons(ind_dubious) = NaN;
N37_cons = N37; N37_cons(ind_dubious) = NaN;
H37_cons = H37; H37_cons(ind_dubious) = NaN;
W37_cons = W37; W37_cons(ind_dubious) = NaN;

%Some less conervative filterings :-
%Just filter for sunglint, thin cirrus and shadow_flag
N16_cons_ignore_aerosol = N16; N16_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
H16_cons_ignore_aerosol = H16; H16_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
W16_cons_ignore_aerosol = W16; W16_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
N21_cons_ignore_aerosol = N21; N21_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
H21_cons_ignore_aerosol = H21; H21_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
W21_cons_ignore_aerosol = W21; W21_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
N37_cons_ignore_aerosol = N37; N37_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
H37_cons_ignore_aerosol = H37; H37_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;
W37_cons_ignore_aerosol = W37; W37_cons_ignore_aerosol(ind_dubious_ignore_aerosol) = NaN;

%Just  filter for thin cirrus and shadow_flag
N16_cons_just_cirrus_shadow = N16; N16_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
H16_cons_just_cirrus_shadow = H16; H16_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
W16_cons_just_cirrus_shadow = W16; W16_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
N21_cons_just_cirrus_shadow = N21; N21_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
H21_cons_just_cirrus_shadow = H21; H21_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
W21_cons_just_cirrus_shadow = W21; W21_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
N37_cons_just_cirrus_shadow = N37; N37_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
H37_cons_just_cirrus_shadow = H37; H37_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;
W37_cons_just_cirrus_shadow = W37; W37_cons_just_cirrus_shadow(ind_dubious_just_cirrus_shadow) = NaN;


if exist('SST')==1 %==1 tests if it is a variable, since this also returns true if it's a directory, file, etc.
    [sst,dates,diff_days] = sst_for_date_range_McCoy_func(SST,sst_time,modis_date_time,'nearest');
    sst_L2 = GRIDDATA(lon_sst,lat_sst,(squeeze(sst)),Plon_L2_5km,Plat_L2_5km);
    %            ctt_L2 = GRIDDATA(Plon_L2_5km,Plat_L2_5km,(squeeze(t_top)),Plon_L2,Plat_L2);
    CTH_L2_5km = (273.15 + sst_L2 - t_top - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km        
    CTH_L2 = GRIDDATA(Plon_L2_5km,Plat_L2_5km,CTH_L2_5km,Plon2_L2,Plat2_L2);
end

%re21_un and tau_un are given in % error (checked in the hdf attributes)
 percent_error_Nd = sqrt( (0.5*tau21_un).^2 + (-5/2*re21_un).^2 );
 re_un_abs = re21_un.*re21/100;
 tau_un_abs = tau21_un.*tau21/100;
 percent_error_Nd16 = sqrt( (0.5*tau16_un).^2 + (-5/2*re16_un).^2 );
 re16_un_abs = re16_un.*re16/100;
 tau16_un_abs = tau16_un.*tau16/100;
 percent_error_Nd37 = sqrt( (0.5*tau37_un).^2 + (-5/2*re37_un).^2 );
 re37_un_abs = re37_un.*re37/100;
 tau37_un_abs = tau37_un.*tau37/100;

%phase_flag = squeeze(qapq_1km(8,:,:));   %1621 phase flag
phase_flag = squeeze(qapp_1km(1,:,:));    %"Primary" phase flag
phase_retreival_outcome = squeeze(qapp_1km(2,:,:)); %Outcome of primary phase retrieval above. 0=not attempted
%(possibly clear, although can actually have a phase for these - but it
%should be ignored). 1=successful, is cloudy. All clear points have =0
%(&phase_flag==1)
surface_flag = squeeze(mask_1km(6,:,:));
tau_bounds = squeeze(qapq_1km(3,:,:)); %tau out of bounds flag


