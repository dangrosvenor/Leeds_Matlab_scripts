function [time,sza2] = get_time_for_SZA2(lat,lon,sza,day)

time=[day:1/(24):day+1];

sza2 = sun_pos(time,lat,lon);

[minval,imin]=min(abs(sza2-sza));


imin2=max([imin 2]);
imin3=min([imin length(time)-1]);

time2=[time(imin2-1):1/(24*60):time(imin3+1)];

sza2 = sun_pos(time2,lat,lon);

[minval,imin]=min(abs(sza2-sza));

time = time2(imin);
sza2 = sza2(imin);
