function[doy,ymd,dec_yr]=dayssince2doy(dayssince,basedate)
JULIAN_BASE=datenum(basedate(1),basedate(2),basedate(3)); 
%% gets days since 0 0 0 for dayssince
dayssince=dayssince+JULIAN_BASE;
%% dayssince is now a datenumber in matlab parlance
[Y, M, D] = datevec(dayssince);
BASE_DATE=datenum(Y,1,1);
doy=dayssince-BASE_DATE+1;
ymd=[Y M D];
dec_yr=Y+doy/366;
end
