

%wrf_dat = nc{'QVAPOR'}(time);

for ilat=1:size(wrf_dat,1);
    for ilon=1:size(wrf_dat,2);

        pdat(1).p(ilat,ilon) = interp1(Z(:,ilat,ilon),wrf_dat(:,ilat,ilon),h_wrf,[],'extrap');

    end
end

