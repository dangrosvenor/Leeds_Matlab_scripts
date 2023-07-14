 LAT_extra = 24.08;
 LON_extra = 38.06;
 
[ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LAT,LON,0.1);
iloc=1;

fid=fopen('E:\wrfruns\Jaafar_Yanbu_nudging\pressure_for_nudging_13thMar09.txt','wt');

for time=1:12
    pres=WRFUserARW(nc,'p',time,ilat(iloc),ilon(iloc));
    Z=WRFUserARW(nc,'Z',time,ilat(iloc),ilon(iloc));

    heights=3000:-100:100;
    P=100*interp1(Z,pres,heights);
    
    for iz=1:length(heights)
        fprintf(fid,'%g \n',P(iz));
    end
    
    fprintf(fid,'\n',P(iz));
    
end

fclose(fid);
    

    