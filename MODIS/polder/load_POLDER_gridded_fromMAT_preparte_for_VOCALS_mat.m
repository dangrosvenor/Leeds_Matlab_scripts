polfile='/home/disk/eos5/d.grosvenor/PARASOL/POLDER_Reff_2005-2012.mat';

%load(polfile);

%set the year to add data for
year_str = '2012'; 
% put that data in the daymean_Par2_CDR array
eval(['daymean_Par2_CDR = daymean_Par2_CDR_' year_str ';']);

%set the lat and lon to the polder ones
MLAT=LAT;
MLON=LON;

% find the averaage data for each day
d=eval( ['meanNoNan(daymean_Par2_MatlabTime_' year_str ',2);'] );
d2=meanNoNan(d,2);

[modisyear_timeseries3,M,D] = datevec(d2);
[daynum_timeseries3] = day_of_year_from_date_func(d2);
