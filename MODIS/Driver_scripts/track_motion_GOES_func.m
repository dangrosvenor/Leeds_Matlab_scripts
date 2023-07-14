function [lat_new,lon_new,lat_old,lon_old] = track_motion_GOES_func(LON_LAT_track,lat_track,lon_track,off_x_LWP,off_y_LWP,gcm_Plat2D_GOES,gcm_Plon2D_GOES)

% Fill in the NaNs - although will have to be careful with this... :-
off_x_LWP = inpaint_nans(off_x_LWP);
off_y_LWP = inpaint_nans(off_y_LWP);


%Interpolate the motion field for the exact location
off_x = griddata(lat_track,lon_track,off_x_LWP,LON_LAT_track(2),LON_LAT_track(1));
off_y = griddata(lat_track,lon_track,off_y_LWP,LON_LAT_track(2),LON_LAT_track(1));

%old position indices
[ilat,ilon]=getind_latlon_quick(gcm_Plat2D_GOES,gcm_Plon2D_GOES,LON_LAT_track(2),LON_LAT_track(1));

lat_old = gcm_Plat2D_GOES(ilat,ilon);
lon_old = gcm_Plon2D_GOES(ilat,ilon);

%ilat_new = round(off_y)+ilat;
%ilon_new = round(off_x)+ilon;

%% Will use the exact index change to estimate a new lat and lon by
%% interpolating assuming a base grid of 1:M and 1:N (see ineterp2 below)

%N.B. - it is correct that y is now for lat and x for lon - since y and x
%were swapped in call to track_motion_Dan in track_motion_GOES_12Nov2008_case.m
ilat_new = off_y+ilat;
ilon_new = off_x+ilon;

% -- use the actual starting lat and lon rather than that estimated from
% indices - just use the indices to get an estimate for change in lat and
% lon for given index change
%ilat_new = off_y + LON_LAT_track(2);
%ilon_new = off_x + LON_LAT_track(1);

%Now interpolate the amount of x and y index movement into an actual lat
%and lon value

%lat_new = gcm_Plat2D_GOES(ilat_new,ilon_new);
%lon_new = gcm_Plon2D_GOES(ilat_new,ilon_new);


%gcm_Plat2D_GOES(isnan(gcm_Plat2D_GOES)) = -20;

%ZI = INTERP2(Z,XI,YI) assumes X=1:N and Y=1:M where [M,N]=SIZE(Z).
% So, even though Z will be of size [550 1200] and M=550 N=1200 it makes X
% 1:N (1:1200) - so have to swap the order of ilat_new and ilon_new
lat_new_temp = interp2(gcm_Plat2D_GOES,ilon_new,ilat_new);
lon_new_temp = interp2(gcm_Plon2D_GOES,ilon_new,ilat_new);


dlat = lat_new_temp - lat_old;
dlon = lon_new_temp - lon_old;

lat_new = LON_LAT_track(2) + dlat;
lon_new = LON_LAT_track(1) + dlon;