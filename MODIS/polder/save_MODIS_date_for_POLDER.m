%save MODIS Date_Time_Swath, daynum_timeseries3_MODIS and
%modisyear_timeseries3_MODIS, ilat_loaded and ilon_loaded

%Arctic box region 2007-2010
savefile_MODIS = '/home/disk/eos5/d.grosvenor/PARASOL/modis_DateTime_for_co_location.mat';

%VOCALS 0-40S, 140W to 68W, DJF 2006 + all of 2007
savefile_MODIS = '/home/disk/eos5/d.grosvenor/PARASOL/modis_DateTime_for_co_location_VOCALS_2007.mat';

%VOCALS 20N-50S, 160W to 60W, all of 2008
savefile_MODIS = '/home/disk/eos5/d.grosvenor/PARASOL/modis_DateTime_for_co_location_VOCALS_2008.mat';

modisyear_timeseries3_MODIS = modisyear_timeseries3;
ilat_loaded_MODIS = ilat_loaded;
ilon_loaded_MODIS = ilon_loaded;

%Date_Time_Swath is of size e.g. [3 51 2880] = [lat lon orbit]. The 2880 of orbit can
%be reshaped into [20 18 4 2] =[Nswaths Ndays Nyears Nsatellites]
sMod = size(Date_Time_Swath.timeseries3);
nswaths = 20;
nprod = sMod(3)/nswaths;
%ndays = length(days_L2L3);
%nyears= length(years);
%nsat = length(direcs);
nlat = size(Date_Time_Swath.timeseries3,1);
nlon = size(Date_Time_Swath.timeseries3,2);

%find the number of days for each year
uni_years = unique(modisyear_timeseries3_MODIS);
freq = ndhistc_run(modisyear_timeseries3_MODIS',[uni_years(1)-0.5:uni_years(end)+0.5]);
freq = freq / nswaths;


daynum_resized_MODIS = reshape(daynum_timeseries3,[nswaths nprod]);
%Date_Time_Swath_MODIS.timeseries3 = reshape(Date_Time_Swath.timeseries3,[nlat nlon nswaths ndays nyears nsat]);
Date_Time_Swath_MODIS.timeseries3 = reshape(Date_Time_Swath.timeseries3,[nlat nlon nswaths nprod]);
%just break down by nswath in case of varying numbers of days per year, or
%per satellite.



save(savefile_MODIS,'freq','daynum_resized_MODIS','Date_Time_Swath_MODIS','modisyear_timeseries3_MODIS','daynum_timeseries3_MODIS','ilat_loaded_MODIS','ilon_loaded_MODIS','-V7.3');

disp('Done save');