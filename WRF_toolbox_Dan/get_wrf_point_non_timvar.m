function val=get_wrf_point_non_timvar(nc,var,time,ilat,ilon)

VAR=nc{var}(time,ilat,:);
val=VAR(ilon)