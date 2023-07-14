function x = sun_pos_fzero(time,lat,lon,sza)

%need an answer that is zero when we have the correct SZA
x = sun_pos(time,lat,lon) - sza;