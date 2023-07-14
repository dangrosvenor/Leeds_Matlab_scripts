function val=get_wrf_point_surface(nc,var,time,ilat,ilon)

if length(time)==1
    VAR=nc{var}(time,ilat,:);
    val=VAR(ilon);
else
    val=nc{var}(time,ilat,ilon);
end
