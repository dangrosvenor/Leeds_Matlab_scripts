%find the (approximate) positions of the flight track on a vertical cross section
%This is done just based on the longitude of the flight track locations

ilim=[1:100:size(dat_flt19,1)];
ilim2=[1:3.161e5];

for i=1:length(ilim)   %size(dat_flt19,1)
    icross(i)=findheight_nearest(lon_slice,dat_flt19(ilim(i),3));
end

for i=1:length(lon_slice)   %size(dat_flt19,1)
    icross_slice(i)=findheight_nearest(dat_flt19(ilim2,3),lon_slice(i));
end
flt19_alt_cross = dat_flt19(icross_slice,11);

%time of the start of the 2250m mini-(temperature)inversion
it=findheight_nearest(time_flt19,19.475); 
lon_lower_warm=dat_flt19(it,3); 
lat_lower_warm=dat_flt19(it,2); 
icross_lower_warm=findheight_nearest(lon_slice,lon_lower_warm);
alt_lower_warm=dat_flt19(it,11);
it_lower_warm=it;

%time of the start of the 2 degree potemp drop
it=findheight_nearest(time_flt19,19.72); 
lon_upper_cool=dat_flt19(it,3); 
lat_upper_cool=dat_flt19(it,2); 
icross_upper_cool=findheight_nearest(lon_slice,lon_upper_cool);
alt_upper_cool=dat_flt19(it,11);
it_upper_cool=it;

%time of the start of the 0.5 degree potemp rise
it=findheight_nearest(time_flt19,19.715); 
lon_upper_warm=dat_flt19(it,3); 
lat_upper_warm=dat_flt19(it,2); 
icross_upper_warm=findheight_nearest(lon_slice,lon_upper_warm);
alt_upper_warm=dat_flt19(it,11);
it_upper_warm=it;

disp('Done find track position')