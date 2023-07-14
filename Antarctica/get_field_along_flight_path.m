function [field,HGT,iz3000,ilat,ilon] = get_field_along_flight_path(nc,lat2d,lon2d,time,LATS,LONS,zfind,field_str,recalc_ilat,ilat_IN,ilon_IN,DX,DY)
% [field,HGT,iz3000] = get_field_along_path(nca,lat2d,lon2d,LATS,LONS,zfind,field_str)
%get values for a field along a particular lat lon path at a particular
%height

if recalc_ilat==1
    [ilat,ilon] = getind_latlon_quick(lat2d.var,lon2d.var,LATS,LONS,0.1);
else
    ilat=ilat_IN;
    ilon=ilon_IN;
end

nZ = length(WRFUserARW(nc,'Z',time,1,1));

field = NaN*ones([1 length(ilat)]);
iz3000 = NaN*ones([1 length(ilat)]);
HGT = NaN*ones([length(ilat) nZ]);

        for iloc=1:length(ilat)
           iloc;
           HGT(iloc,:)=WRFUserARW(nc,'Z',time,ilat(iloc),ilon(iloc));
           iz3000(iloc) = findheight_nearest(HGT(iloc,:),zfind);
           switch field_str
               case 'Pressure'
                   prof = ( nc{'P'}(time,:,ilat(iloc),ilon(iloc)) + nc{'PB'}(time,:,ilat(iloc),ilon(iloc)) ) / 100; %pressure                   
%                   field(iloc) = ( nc{'P'}(time,iz3000(iloc),ilat(iloc),ilon(iloc)) + nc{'PB'}(time,iz3000(iloc),ilat(iloc),ilon(iloc)) ) / 100; %pressure
                        %Using the above approach of finding the nearest
                        %level gives spurious results presumably because
                        %the level can jump. Interpolation (below) works
                        %better
               case 'Wind direction'
                   u=WRFUserARW(nc,'u',time,ilat(iloc),ilon(iloc));
                   v=WRFUserARW(nc,'v',time,ilat(iloc),ilon(iloc));
                   prof = wind_dir_compass_from_uv_wrf(u,v,lat2d,lon2d,ilat(iloc),ilon(iloc),DX,DY);

                   
                   
           end
           
           field(iloc) = interp1(HGT(iloc,:),prof,zfind);
                              
        end