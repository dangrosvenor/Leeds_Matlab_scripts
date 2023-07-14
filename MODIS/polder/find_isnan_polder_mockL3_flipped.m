%Find the instances where POLDER is NaN and return those indices into isnan

ilat_POL = find(LAT_MODIS>=LAT_val(1) & LAT_MODIS<LAT_val(end));
ilon_POL = find(LON_MODIS>=LON_val(1) & LON_MODIS<LON_val(end));
days = daynum_timeseries3(itime);
years = modisyear_timeseries3(itime);
itime_POL = NaN*ones([1 length(times)]);
for ifind=1:length(days)
    itime_POL(ifind) = find(daynum_timeseries3_MODIS==days(ifind) & modisyear_timeseries3_MODIS==years(ifind));
end
%            dat = flipdim(daymean_Par2_CDR,1);
dat = daymean_Par2_CDR(ilat_POL,ilon_POL,itime_POL);
%Now X and dat should be the same except that dat is flipped
%in the lat dimension (1st dim)
dat = flipdim(dat,1);
%X(isnan(dat))=NaN;
inan_POL = isnan(dat);
