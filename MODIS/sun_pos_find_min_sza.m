function [sza_min, time_min]=sun_pos_find_min_sza(date_str,lat,lon)
% [sza_min, time_min]=sun_pos_find_min_sza(date_str,lat,lon)
% Finds the minimum SZA for the date supplied to date_str (in format
% '12-Nov-2008'. Also returns the time (UTC) of the min.

times=datenum(date_str)+[1:1/60:24]/24;
szas=sun_pos(times,-20,-76);

[a,b]=minALL(szas);

sza_min = a;
time_min = datestr(times(b(2)));

