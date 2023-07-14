%finds the range of SZA values for each latitude for a give period

ilats=1:90; %30=60.5 degrees
day_str='21-Dec-2005';
iday2=datenum(day_str)-datenum('01-Jan-2005')+1;

clear sza_range_equals2

ndays=16;
ndays=1;
iday_count=0;
for iday=iday2:iday2+ndays-1
iday_count=iday_count+1;

%number of days to include
nday=1;
iday=iday:iday+nday-1;

minsza = NaN*ones([1 length(ilats)]);
maxsza = NaN*ones([1 length(ilats)]);
Np = NaN*ones([1 length(ilats)]);

for ilat=1:length(ilats)

mindat=Solar_Zenith_Minimum.timeseries3(ilat,:,[iday iday+365]); mindat=mindat(:);
maxdat=Solar_Zenith_Maximum.timeseries3(ilat,:,[iday iday+365]); maxdat=maxdat(:);
sensdat=Sensor_Zenith_Maximum.timeseries3(ilat,:,[iday iday+365]); sensdat=sensdat(:);

iyes=find(abs(mindat-maxdat)<=2 & sensdat<30);
Np(ilat)=length(iyes);

if Np(ilat)>0
    minsza(ilat)=min(mindat(iyes));
    maxsza(ilat)=max(maxdat(iyes));
end





%fprintf(1,'\nNp=%d, minsza=%f, maxsza=%f, LAT=%f\n',Np(ilat),minsza(ilat),maxsza(ilat),LAT(ilat));




end

fprintf(1,'\nDone\n');



sza_range = maxsza-minsza;
inot_nan=find(isnan(sza_range)==0);
sza_range2 = sza_range(inot_nan);

y=maxsza(inot_nan);
[x,I]=unique(sza_range2);
y=y(I);
sza_range_equals2(iday_count) = interp1(x,y,3);

x=y;
y=LAT(ilats(inot_nan));
y=y(I);
lat_val(iday_count) = interp1(x,y,sza_range_equals2(iday_count));

end


break

%%% ******        ******* %%%%


savedir='~/modis_work/plots/';

figure
plot(maxsza-minsza,maxsza,'ro-'); grid
xlabel_str = 'SZA range';
ylabel_str = 'Max SZA';
xlabel(xlabel_str);
ylabel(ylabel_str);
title(day_str);

savename=[savedir ylabel_str ' vs ' xlabel_str ' for ' day_str];



figure
plot(Np,LAT(ilats),'ro-'); grid;
xlabel_str = 'N datapoints';
ylabel_str = 'Max SZA';
xlabel(xlabel_str);
ylabel(ylabel_str);
title(day_str);

savename=[savedir ylabel_str ' vs ' xlabel_str ' for ' day_str];








