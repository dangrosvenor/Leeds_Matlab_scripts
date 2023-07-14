years=2005:2012;
Lyears=[365 365 365 366 365 365 365 366];
ii=1;
for i=1:length(years)
    daynum_timeseries3_MODIS(ii:ii+Lyears(i)-1)=[1:Lyears(i)];
    modisyear_timeseries3_MODIS(ii:ii+Lyears(i)-1)=years(i);    
    ii=ii+Lyears(i);    
end


