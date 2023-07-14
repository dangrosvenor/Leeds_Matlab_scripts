function [ilats,ilons,dmins]=am3_map_lat_lon_points(am3_lat,am3_lon,Plat,Plon)
%function [ilats,ilons,dmins]=am3_map_lat_lon_points(am3_lat,am3_lon,Plat,Plon)
 
max_ov=-9e99;
for i=1:length(Plat(:))

         dist = distlatlon(Plat(i),Plon(i),am3_lat,am3_lon);
         [minval,idist]=minALL(dist);
         
         ilats(i) = idist(1);
         ilons(i) = idist(2);
%         lats(i) = am3_lat(idist(1),idist(2));
%         lons(i) = am3_lon(idist(1),idist(2));            
         
         dmin = minALL(dist);
         if dmin>max_ov
             max_ov=dmin;
             imax=i;
             jmax=j;
         end
         
         dmins(i)=dmin;

 end
 
 