%sets the correct column numbers for the aircraft data depending on the structure of the array
if size(dat_flt,2)==15
    %for flt_19
    col_temp=6;
    col_alt=11;
    col_lat=2;
    col_lon=3;
    col_press=4;
    col_wind=9;
    col_winddir=10;
else
    %for Feb2010 flights
    col_temp=5;
    col_alt=12;
    col_lat=2;
    col_lon=3;
    col_press=6;
    col_wind=9;
    col_winddir=10;
    col_frostpoint_hygro=7;
    col_frostpoint_humi=8;
    col_airspeed=4;
end