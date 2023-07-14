%function height=get_met_height_all_points(nc)
clear pres_2300

lat2d.var=nc{'XLAT_M'}(:);
lon2d.var=nc{'XLON_M'}(:);

if prod(size(nc{'GHT'})) > 0
    height = nc{'GHT'}(1,:);
else

    for ilat=1:size(lat2d.var,1)
        fprintf(1,'ilat=%d\n',ilat);
        for ilon=1:size(lat2d.var,2)
            pres_2300(ilat,ilon)=get_met_height(nc,ilat,ilon,1,2300,500e2);
        end
    end

end


