%Load Ryan's CTh product and match in time and space to loaded L3 MODIS
%data

ryan_dir = '/home/disk/eos8/d.grosvenor/CTH_Ryan/';
%Just do daytime for now
ryan_file = 'Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_DG_hif_zb.nc.mat';

load([ryan_dir ryan_file]);

%Now match to the data in modisyear_timeseries3_MODIS,
%daynum_timeseries3_MODIS, MLAT and MLON

%Both MLAT and Ryans's lat run from high lats to low lats.
ilat=find(Ryan_day.lat<=MLAT(1) & Ryan_day.lat>=MLAT(end));
%MLON and Ryan's lon run from low to high
ilon=find(Ryan_day.lon>=MLON(1) & Ryan_day.lon<=MLON(end));

CTH_Ryan = NaN*ones([length(MLAT) length(MLON) length(daynum_timeseries3_MODIS)]);

for it=1:length(modisyear_timeseries3_MODIS)
   itf = find(Ryan_day.years==modisyear_timeseries3_MODIS(it) & Ryan_day.days==daynum_timeseries3_MODIS(it) );
   if length(itf)>0
       CTH_Ryan(:,:,it) = Ryan_day.cth(ilat,ilon,itf);
   end
end
