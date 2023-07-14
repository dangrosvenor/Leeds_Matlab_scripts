%load John King's flight data

filename_ncdf = '/home/disk/eos1/d.grosvenor/Ant_flights/flight19/flight19.ncdf';
nc19 = netcdf(filename_ncdf);

time19 = nc19{'time'}(:); %milliseconds
   time19 = time19 /1e3 / 3600 /24; %convert to days
lon19 = nc19{'lon'}(:);
lat19 = nc19{'lat'}(:);
temp19 = nc19{'atn'}(:); %deg
p19 = nc19{'ps'}(:); %hPa
   p19 = p19 *1e2;
u19 = nc19{'u'}(:);
  u19(u19>100) = NaN;
v19 = nc19{'v'}(:);
w19 = nc19{'w'}(:);
alt19 = nc19{'altg'}(:); %GPS altitude (m)

pot19 = potemp(temp19+273.15,p19);

%find indices for the leg at 2900m at around 67S for comparison to WRF
[temp1,temp2] = findheight_nearest(time19,0.8138,0.8432);
inds_cs = temp1:temp2;


%average the flight data into 2km chunks since that is the size of the WRF
%gridboxes
mean_lat = meanNoNan(lat19(inds_cs),1);
dfs = distlatlon(mean_lat,lon19(inds_cs(1)),mean_lat,lon19(inds_cs));

chunk = 2; %averaging length (km)
nchunks = floor(max(dfs)/chunk);

pot19av = NaN * ones([nchunks 1]);
lon19av = NaN * ones([nchunks 1]);
u19av = NaN * ones([nchunks 1]);
w19av = NaN * ones([nchunks 1]);
for i=0:nchunks-1
   inds = find(dfs > i*chunk & dfs <= (i+1)*chunk);
   pot19av(i+1) = meanNoNan(pot19(inds_cs(inds)),1);  
   lon19av(i+1) = meanNoNan(lon19(inds_cs(inds)),1);  
   u19av(i+1) = meanNoNan(u19(inds_cs(inds)),1);  
   w19av(i+1) = meanNoNan(w19(inds_cs(inds)),1);     
end


% 
% netcdf flight19 {
% dimensions:
%         length = 787105 ;
%         one = 1 ;
% variables:
%         int time(length) ;
%                 time:units = "milliseconds" ;
%                 time:name = "Time in day" ;
%         float lat(length) ;
%                 lat:units = "degrees" ;
%                 lat:name = "Latitude" ;
%         float lon(length) ;
%                 lon:units = "degrees" ;
%                 lon:name = "Longitude" ;
%         float atn(length) ;
%                 atn:units = "Degrees C" ;
%                 atn:name = "Air Temperature" ;
%         float u(length) ;
%                 u:units = "ms-1" ;
%                 u:name = "Wind u-component" ;
%         float v(length) ;
%                 v:units = "ms-1" ;
%                 v:name = "Wind v-component" ;
%         float w(length) ;
%                 w:units = "ms-1" ;
%                 w:name = "Wind w-component" ;
%         float tf(length) ;
%                 tf:units = "Degrees C" ;
%                 tf:name = "Fast temperature" ;
%         float irt(length) ;
%                 irt:units = "Degrees C" ;
%                 irt:name = "Surface Temperature" ;
%         float swo(length) ;
%                 swo:units = "Wm-2" ;
%                 swo:name = "Short wave outgoing" ;
%         float swi(length) ;
%                 swi:units = "Wm-2" ;
%                 swi:name = "Short wave incoming" ;
%         float iri(length) ;
%                 iri:units = "Wm-2" ;
%                 iri:name = "IR incoming" ;
%         float iro(length) ;
%                 iro:units = "Wm-2" ;
%                 iro:name = "IR outgoing" ;
%         float ps(length) ;
%                 ps:units = "hPa" ;
%                 ps:name = "Static pressure" ;
%         float altg(length) ;
%                 altg:units = "m" ;
%                 altg:name = "GPS height" ;
%         float altr(length) ;
%                 altr:units = "m" ;
%                 altr:name = "Radar height" ;
%         float dpc(length) ;
%                 dpc:units = "degree C" ;
%                 dpc:name = "frost point" ;
%         float altl(length) ;
%                 altl:units = "m" ;
%                 altl:name = "Laser height" ;
%         float head(length) ;
%                 head:units = "deg true" ;
%                 head:name = "Aircraft heading" ;
%         float tas(length) ;
%                 tas:units = "m/s" ;
%                 tas:name = "True airspeed" ;
%         float qlic(length) ;
%                 qlic:units = "millimoles/mole" ;
%                 qlic:name = "LiCor mixing ratio" ;
%         float miss(one) ;
%                 miss:units = "none" ;
%                 miss:name = "missing data;" ;
% }
