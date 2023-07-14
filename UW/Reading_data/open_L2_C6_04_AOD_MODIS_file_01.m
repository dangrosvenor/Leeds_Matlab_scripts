%N.B. MODIS files are hdf4 not hdf5, so need to use hdfinfo, etc rather
%than hdf5info

try %catch all errors later on - so that override_flags_openL2 flag can be
    %reset even if there is an error. The error is then rethrown
    %(re-issued).



if ~exist('override_flags_openL2') | override_flags_openL2==0
    day_only=0;
end

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
    filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/'
    

    %    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0000.005.2006278184700.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0005.005.2006278184047.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2355.005.2006278230055.hdf';
%    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2220.005.2006278224532.hdf';
    file_name_h5='MPACE_L2_17days/MOD06_L2.A2004277.0020.051.2010289113334.hdf';
%9th Oct MPACE 21:15 swath - high Nd point identified, but also some normal ones    
    file_name_h5='MPACE_L2_17days/MOD06_L2.A2004283.2115.051.2010289225041.hdf';  %the usual swath
    
   file_name_h5='MOD_L2/2012_06_05_hole_Australia/MYD06_L2.A2012157.0500.051.2012158011806.hdf';
   

   file_name_h5='MOD_L2/2012_06_05_hole_Australia/MYD06_L2.A2012157.0455.051.2012158011933.hdf'; %main image of the hole
%   file_name_h5='MOD_L2/2012_06_05_hole_Australia/MYD06_L2.A2012157.0450.051.2012158011035.hdf';

%    file_name_h5='MPACE_L2_17days_all_lats/aqua/MYD06_L2.A2004283.2130.051.2009014155700.hdf';
%    file_name_h5='MPACE_L2_17days_all_lats/aqua/MYD06_L2.A2004283.2135.051.2009014154943.hdf';
% filedir= [filedir 'AP_Dec_2009_daytime/AQUA/']; file_name_h5='MYD06_L2.A2009355.1900.051.2009357131942.hdf';
% filedir= [filedir 'summer_2004_17days_all_lats/aqua/']; file_name_h5='MYD06_L2.A2004181.2345.051.2009008183025.hdf'; 

% filedir= [filedir 'AP_Feb_2010_daytime/AQUA/']; file_name_h5='MYD06_L2.A2010043.1915.051.2010044195639.hdf'; 
% filedir= [filedir 'Arctic_summerL2_2007_20W-60E_70-80N/terra/']; file_name_h5='MOD06_L2.A2007175.1320.051.2010319152819.hdf';  
% file_name_h5='Arctic_summerL2_2007_20W-60E_70-80N/Joint_5km_files/terra/MODATML2.A2007166.0515.051.2010319070803.hdf';  %the usual swath    
 
%  filedir= [filedir 'MOD_L2/April2nd_2005_ZhangPlatnick_2011_L2_scene/']; file_name_h5='MOD06_L2.A2005092.1600.051.2010295045517.hdf';  
 
%filedir='/home/disk/eos8/d.grosvenor/';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1350.051.2010332165733.hdf';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1705.051.2010332163313.hdf';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1710.051.2010332162345.hdf';    

file_name_h5='AP_Feb_2010_daytime/06Feb2010_flight99/aqua/MYD06_L2.A2010037.1815.051.2010038223324.hdf'; %flight 99 Aqua 18:15 image
file_name_h5='MOD_L2/Finland_Sami/aqua/2006/MYD06_L2.A2006260.0935.051.2009069120820.hdf';   %17th Sep, 2006
file_name_h5='MOD_L2/Finland_Sami/terra/2006/MOD06_L2.A2006260.0920.051.2010310012604.hdf';   %17th Sep, 2006
%file_name_h5='MOD_L2/Finland_Sami/aqua/2006/MYD06_L2.A2006260.0755.051.2009069115237.hdf';   %17th Sep, 2006
%file_name_h5='MOD_L2/Finland_Sami/terra/2006/MOD06_L2.A2006260.1055.051.2010310011457.hdf';   %17th Sep, 2006, 10:55

%file_name_h5='MOD_L2/Finland_Sami/aqua/2006/MYD06_L2.A2006283.1120.051.2009071015900.hdf';   %10th Oct, 2006
%file_name_h5='MOD_L2/Finland_Sami/aqua/2006/MYD06_L2.A2006283.1115.051.2009071014057.hdf';   %10th Oct, 2006
%file_name_h5='MOD_L2/Finland_Sami/terra/2006/MOD06_L2.A2006278.0905.051.2010310121601.hdf';   %5th Oct, 2006
%file_name_h5='MOD_L2/Finland_Sami/MYD06_L2.A2006278.1045.051.2009070154100.hdf';   %5th Oct, 2006 - AQUA
%file_name_h5='MOD_L2/Finland_Sami/terra/2006/MOD06_L2.A2006275.1015.051.2010310111527.hdf';   %2nd Oct, 2006 
%file_name_h5='Arctic_summerL2_2007_20W-60E_70-80N/terra/MOD06_L2.A2007164.1025.051.2010319065119.hdf';
%file_name_h5='MOD_L2/Finland_Sami/MYD06_L2.A2006130.1045.051.2009059144147.hdf';   %10th May, 2006 - AQUA
file_name_h5='MOD_L2/Finland_Sami/MOD06_L2.A2006278.1045.051.2010310122952.hdf';   %5th Oct, 2006 - TERRA

%file_name_h5='/Leeds/MOD06_L2.A2014249.1155.006.2015077042234.hdf';   %
%file_name_h5='MYD06_L2.A2015282.1215.051.NRT.hdf';   %

file_name_h5='MOD_L2/C6/Hamish/MOD06_L2.A2016214.1040.006.2016215022430.hdf'; %Hamish CLARIFY test files
file_name_h5='MOD_L2/C6/Jesus/MYD06_L2.A2014343.1325.006.2014344210847.hdf'; %Hamish CLARIFY test files
file_name_h5='MOD_L2/NAtlantic/MYD06_L2.A2016214.0230.006.2016215074529.hdf'; %ACSIS 1st Aug 2016 case
file_name_h5='MOD_L2/ship_tracks_UK_Spain_NdReview/MYD06_L2.A2017357.1315.006.2017361230642.hdf'; %Ship track Nd cases
%file_name_h5='MOD_L2/ship_tracks_UK_Spain_NdReview/MYD06_L2.A2017357.1320.006.2017361230012.hdf'; %Ship track Nd cases

file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';

end

filename_h5=[filedir file_name_h5];


iday=findstr(file_name_h5,'.A');
modis_year_str=file_name_h5(iday+2:iday+5);
modis_day_str=file_name_h5(iday+6:iday+8);
modis_time_str=file_name_h5(iday+10:iday+13);
modis_hour_str=modis_time_str(1:2);
modis_min_str=modis_time_str(3:4);
aq_terr_str = file_name_h5(iday-8:iday-1);
date_str=datestr(datenum(['01-Jan-' modis_year_str])+str2num(modis_day_str)-1,1);
modis_date_time = datenum([date_str ' ' modis_hour_str ':' modis_min_str]);

if day_only==1
    if (str2num(modis_time_str)>1800 | str2num(modis_time_str)<800)
        night_time=1;
        return
    end
end


%% Start of read


switch comp
    case {'uni','UWchallenger'}

        INFO = hdfinfo(filename_h5,'eos');
        %N.B. this step (INFO) doesn't appear to be needed to read the data
        %anyway.

        %    INFO.GroupHierarchy.Datasets.Name %gives the info of all the dataset names

        %INFO.GroupHierarchy.Datasets(1)  %this gives all the info about dataset 1

        %    NVARS = length(INFO.GroupHierarchy.Datasets); %number of
        %    variables

        %    test = hdf5read(INFO.GroupHierarchy.Datasets(NVARS-1)); %this retrieves dataset 10



end

SD_id = hdfsd('start',filename_h5,'read'); %open the file


%to look at attributes do e.g. :
% [temp,dimsizes_temp]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO,0,1);

%typical attributes
% 1) _FillValue
% 2) long_name
% 3) units
% 4) scale_factor
% 5) add_offset
% 6) Parameter_Type
% 7) Cell_Along_Swath_Sampling
% 8) Cell_Across_Swath_Sampling
% 9) Geolocation_Pointer
% 10) Description (only for QA flags - useful info on how the bits work -see end of this script)

[lat,dimsizes_lat,sampling_along,sampling_across]=get_hdf_data_dan('Latitude',SD_id,INFO);
[lon,dimsizes_lon]=get_hdf_data_dan('Longitude',SD_id,INFO);

[AOD_550_Dark_Target_Deep_Blue_Combined,dimsizes_lon]=get_hdf_data_dan('AOD_550_Dark_Target_Deep_Blue_Combined',SD_id,INFO);
[AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag,dimsizes_lon]=get_hdf_data_dan('AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',SD_id,INFO);
%The above is a QA flag for the product: 3 = good, 0=bad.


[scantime,dimsizes_scantime]=get_hdf_data_dan('Scan_Start_Time',SD_id,INFO); %double
% (seconds since 1993-1-1 00:00:00.0) (!)
[solar_zenith,dimsizes_solar_zenith]=get_hdf_data_dan('Solar_Zenith',SD_id,INFO); %int16, min=-32767
[solar_azimuth,dimsizes_sensor_zenith]=get_hdf_data_dan('Solar_Azimuth',SD_id,INFO); %int16
[sensor_zenith,dimsizes_sensor_zenith]=get_hdf_data_dan('Sensor_Zenith',SD_id,INFO); %int16
[sensor_azimuth,dimsizes_sensor_zenith]=get_hdf_data_dan('Sensor_Azimuth',SD_id,INFO); %int16
%strange minima with solar and sensor ZA at the end of the swath


%[BTD,dimsizes_temp]=get_hdf_data_dan('Brightness_Temperature_Difference',SD_id,INFO);
%[BT,dimsizes_temp]=get_hdf_data_dan('Brightness_Temperature',SD_id,INFO);
%[rad_var,dimsizes_temp]=get_hdf_data_dan('Radiance_Variance',SD_id,INFO); %for same bands as BT (see below)
%Bands 29,31,32,33,34,35,36
%Wavelength 8.4-8.7 (29), 10.8-11.3 (31), 11.8-12.3 (32), 13.2-13.5 (33),
%13.5-13.8 (34) 13.8-14.1 (35), 14.1-14.4 (36) um

%% Close the HDF file - otherwise get problems if open too many files
hdfsd('end',SD_id);

%convert to Matlab date time (days since 01-Jan-0000)
scantime_matlab = scantime/3600/24 + datenum('01-Jan-1993');
%then can do datestr(scantime_matlab(ilat,ilon),31) to get the date string


catch openL2_error
   clear override_flags_openL2 %reset the flag if there was an error
   rethrow(openL2_error); %re-issue the error
end



    %Read 1km mask and flag data
disp('Done read M?D04_L2 MODIS');


