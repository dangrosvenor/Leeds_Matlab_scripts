function [out,lats_mid,Nout,std_out]=UM_calc_lat_mean(lat2d,lon2d,dat,lat_range,nlat)

if length(lat_range)==0
    lat_range(1)=minALL(lat2d);
    lat_range(2)=maxALL(lat2d);    
end

dlat = (lat_range(2)-lat_range(1) ) / nlat;

lats=[lat_range(1):dlat:lat_range(2)];
lats_mid=0.5*(lats(1:end-1)+lats(2:end));

for ilat=1:length(lats)-1
    ii = find(lat2d>lats(ilat) & lat2d<=lats(ilat+1));
    [out(ilat),Nout(ilat),std_out(ilat)] = meanNoNan(dat(ii),1); 
end



