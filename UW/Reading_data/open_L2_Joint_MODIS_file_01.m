%read an L2 joint file (reduced dataset from L2 at 5 km resolution)

region_flag='not epic';

if ~exist('comp')
    comp='uni';
    comp='laptop';
    comp='UWchallenger';
end

if exist('imodis_file_override') & imodis_file_override==1
    clear imodis_file_override  %don't set file_name_h5 and reset for next run
else
    filedir='/home/disk/eos1/d.grosvenor/';
 filedir='/home/disk/eos8/d.grosvenor/';    
 
    %    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0000.005.2006278184700.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0005.005.2006278184047.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2355.005.2006278230055.hdf';
%    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2220.005.2006278224532.hdf';
    file_name_h5='MPACE_JointL2_17days_all_lats/aqua/MODATML2.A2004277.0020.051.2010289113442.hdf';
end

filename_h5=[filedir file_name_h5];


iday=findstr(file_name_h5,'.A');
modis_year_str=file_name_h5(iday+2:iday+5);
modis_day_str=file_name_h5(iday+6:iday+8);
modis_hour_str=file_name_h5(iday+10:iday+11);
modis_min_str=file_name_h5(iday+12:iday+13);

date_str = datestr(datenum(['01-Jan-' modis_year_str])+str2num(modis_day_str)-1);
start_time = datenum([date_str ' ' modis_hour_str ':' modis_min_str ':00']);


switch comp
    case {'uni','UWchallenger'}

        INFO = hdfinfo(filename_h5);
end

SD_id = hdfsd('start',filename_h5,'read'); %open the file


[lat,dimsizes_lat]=get_hdf_data_dan('Latitude',SD_id,INFO,1);
[lon,dimsizes_lon]=get_hdf_data_dan('Longitude',SD_id,INFO,1);
[cf,dimsizes_cf]=get_hdf_data_dan('Cloud_Fraction',SD_id,INFO,1);
[t_top,dimsizes_ctt]=get_hdf_data_dan('Cloud_Top_Temperature',SD_id,INFO,1);
[p_top,dimsizes_ctp]=get_hdf_data_dan('Cloud_Top_Pressure',SD_id,INFO,1);
[tau,dimsizes_tau]=get_hdf_data_dan('Cloud_Optical_Thickness',SD_id,INFO,1);
[qa5,dimsizes_qa5]=get_hdf_data_dan('Cloud_Quality_Assurance',SD_id,INFO,1);
[mask5,dimsizes_mask5]=get_hdf_data_dan('Cloud_Mask',SD_id,INFO,1);
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO,1); %Rob had this as Effective_Particle_Radius?





qa5=flipdim(qa5,2);
mask5=flipdim(mask5,2);

cf=flipdim(cf,2);
tau=0.1*flipdim(tau,2);
re=0.1*flipdim(re,2);

lat=flipdim(lat,2);
lon=flipdim(lon,2);

t_top=0.01*flipdim(t_top+15000,2);  %L?  %the L in Rob's IDL script means long integer needed to make sure that the values
   %can get high enough. Shouldn't be a problem in Matlab because I convert them all to double - check!
p_top=0.1*flipdim(p_top,2);


%have saved the array scantime_L2 for an L2 swath - this should be the same
%for each swath - so here just add this on to the start time of the swath
scantime_matlab = start_time + scantime_L2;
sensor_zenith = sensor_zenith_L2;
%sun = sun_position(time,location);
%no. hours since start of the day
time_hours = (scantime_matlab-datenum(date_str))*24;
solar_zenith = zenith(day_arr,time_hours,lat,lon);



disp('Done read L2 MODIS');

