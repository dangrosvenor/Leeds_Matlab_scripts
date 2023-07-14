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
    %    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0000.005.2006278184700.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.0005.005.2006278184047.hdf';
    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2355.005.2006278230055.hdf';
%    file_name_h5='ladsweb.nascom.nasa.gov/allData/5/MOD06_L2/2000/355/MOD06_L2.A2000355.2220.005.2006278224532.hdf';
    file_name_h5='MPACE_L2_17days/MOD06_L2.A2004277.0020.051.2010289113334.hdf';
%9th Oct MPACE 21:15 swath - high Nd point identified, but also some normal ones    
    file_name_h5='MPACE_L2_17days/MOD06_L2.A2004283.2115.051.2010289225041.hdf';  %the usual swath
   
    
%    file_name_h5='MPACE_L2_17days_all_lats/aqua/MYD06_L2.A2004283.2130.051.2009014155700.hdf';
%    file_name_h5='MPACE_L2_17days_all_lats/aqua/MYD06_L2.A2004283.2135.051.2009014154943.hdf';
% filedir= [filedir 'AP_Dec_2009_daytime/AQUA/']; file_name_h5='MYD06_L2.A2009355.1900.051.2009357131942.hdf';
% filedir= [filedir 'summer_2004_17days_all_lats/aqua/']; file_name_h5='MYD06_L2.A2004181.2345.051.2009008183025.hdf'; 

% filedir= [filedir 'AP_Feb_2010_daytime/AQUA/']; file_name_h5='MYD06_L2.A2010043.1915.051.2010044195639.hdf'; 
% filedir= [filedir 'Arctic_summerL2_2007_20W-60E_70-80N/terra/']; file_name_h5='MOD06_L2.A2007175.1320.051.2010319152819.hdf';  
 file_name_h5='Arctic_summerL2_2007_20W-60E_70-80N/Joint_5km_files/terra/MODATML2.A2007166.0515.051.2010319070803.hdf';  %the usual swath    
 
%  filedir= [filedir 'MOD_L2/April2nd_2005_ZhangPlatnick_2011_L2_scene/']; file_name_h5='MOD06_L2.A2005092.1600.051.2010295045517.hdf';  
 
%filedir='/home/disk/eos8/d.grosvenor/';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1350.051.2010332165733.hdf';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1705.051.2010332163313.hdf';
%file_name_h5='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/MOD06_L2.A2008322.1710.051.2010332162345.hdf';    

%these two don't contain lat_10km and solar_zenith
file_name_h5='joint_L2/terra/2007/001/MODATML2.A2007001.0000.051.2010312202426.hdf';  
file_name_h5='joint_L2/terra/2007/001/MODATML2.A2007001.2355.051.2010312195516.hdf';  

%this one doesn't work
file_name_h5='joint_L2/terra/2006/306/MODATML2.A2006306.2355.051.2010311171442.hdf';  

%this from day 002 of 2007 works!
file_name_h5='joint_L2/terra/2007/002/MODATML2.A2007002.2245.051.2010313012329.hdf';  


%this works
file_name_h5='Arctic_summerL2_2007_20W-60E_70-80N/Joint_5km_files/terra/MODATML2.A2007181.2235.051.2010319175252.hdf';  

filedir='/home/disk/eos8/d.grosvenor/joint_L2/terra/'; 
file_name_h5='2007/032/MODATML2.A2007032.1755.051.2010316050338.hdf';
%as does this???
%file_name_h5='joint_L2/aqua/2007/001/MYDATML2.A2007001.0005.051.2009076122318.hdf';

%maybe the nighttime ones don't work for some reason??

%nighttime (all points are nighttime) - no solar_zenith data
%filedir='/home/disk/eos5/d.grosvenor/joint_L2/'; 
filedir='/home/disk/eos4/d.grosvenor/AMSRE_L2/'; 


file_name_h5='terra/2007/019/MODATML2.A2007019.2250.051.2010313181253.hdf';
%swath with some daytime (next swath) - does contain solar_zenith
file_name_h5='terra/2007/019/MODATML2.A2007019.2255.051.2010313214825.hdf';
file_name_h5='aqua/2010/179/MYDATML2.A2010179.0105.051.2010179211536.hdf';

%so if the swath contains no pixels with SZA<81 then the solar_zenith field
%doesn't exist (is NaN). lat and lon do exist though.

%used for checking the relative azimuth angle that is given by MODIS
file_name_h5='terra/2007/164/MODATML2.A2007164.1025.051.2010319065205.hdf';   %13th June, 2007, Terra

file_name_h5='Puijo/terra/2006/MODATML2.A2006260.1055.051.2010310011535.hdf';   %13th June, 2007, Terra

file_name_h5='Puijo/aqua/2006/MYDATML2.A2006333.0925.051.2009074073053.hdf';   %29th Nov, 2006, Aqua


file_name_h5='AE_Ocean.002/2007.01.01/AMSR_E_L2_Ocean_V06_200701010748_A.hdf';   %AMSRE test file

end

filename_h5=[filedir file_name_h5];


iday=findstr(file_name_h5,'_V')+3;
modis_year_str=file_name_h5(iday+2:iday+5);
modis_month_str = file_name_h5(iday+6:iday+7);
modis_dayofmonth_str = file_name_h5(iday+8:iday+9);
[temp,modis_day_str] = day_of_year_from_date_func([modis_dayofmonth_str '.' modis_month_str '.' modis_year_str],'dd.mm.yyyy');
modis_time_str=file_name_h5(iday+10:iday+13);
modis_hour_str=modis_time_str(1:2);
modis_min_str=modis_time_str(3:4);
aq_terr_str = 'AQUA'; %file_name_h5(iday-8:iday-1);
date_str=datestr(datenum(['01-Jan-' modis_year_str])+str2num(modis_day_str)-1,1);
modis_date_time = datenum([date_str ' ' modis_hour_str ':' modis_min_str]);

if day_only==1
    if (str2num(modis_time_str)>1800 | str2num(modis_time_str)<800)
        night_time=1;
        return
    end
end





switch comp
    case {'uni','UWchallenger'}

        INFO = hdfinfo(filename_h5);

        %    INFO.GroupHierarchy.Datasets.Name %gives the info of all the dataset names

        %INFO.GroupHierarchy.Datasets(1)  %this gives all the info about dataset 1

        %    NVARS = length(INFO.GroupHierarchy.Datasets); %number of
        %    variables

        %    test = hdf5read(INFO.GroupHierarchy.Datasets(NVARS-1)); %this retrieves dataset 10



end

SD_id = hdfsd('start',filename_h5,'read'); %open the file


%to look at attributes do e.g. :
% [temp,dimsizes_temp]=get_hdf_data_AMSRE('Cloud_Effective_Radius',SD_id,INFO,0,1);

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

[lat,dimsizes_lat,sampling_along,sampling_across]=get_hdf_data_AMSRE('Latitude',SD_id,INFO);
[lon,dimsizes_lon]=get_hdf_data_AMSRE('Longitude',SD_id,INFO);
%[lat_10km,dimsizes_lat,sampling_along_10km,sampling_across_10km,success]=get_hdf_data_AMSRE('Latitude_10km',SD_id,INFO);

[qa_L2amsre]=get_hdf_data_AMSRE('Ocean_products_quality_flag',SD_id,INFO); 
[sst_vlow_L2amsre]=get_hdf_data_AMSRE('Very_low_res_sst',SD_id,INFO); %
[sst_L2amsre]=get_hdf_data_AMSRE('Low_res_sst',SD_id,INFO); %
[wind_low_L2amsre]=get_hdf_data_AMSRE('Low_res_wind',SD_id,INFO); %
[wind_med_L2amsre]=get_hdf_data_AMSRE('Med_res_wind',SD_id,INFO); 
[vap_med_L2amsre]=get_hdf_data_AMSRE('Med_res_vapor',SD_id,INFO); %
[lwp_L2amsre]=get_hdf_data_AMSRE('High_res_cloud',SD_id,INFO); %

%if success==0 %if it failed to read Latitude_10km
        
%end

%[lon_10km,dimsizes_lon]=get_hdf_data_AMSRE('Longitude_10km',SD_id,INFO);
%[cf,dimsizes_cf]=get_hdf_data_AMSRE('Cloud_Fraction',SD_id,INFO); %int8 (0...127)
%[t_top,dimsizes_ctt]=get_hdf_data_AMSRE('Cloud_Top_Temperature',SD_id,INFO); %int16 (>0)
%[p_top,dimsizes_ctp]=get_hdf_data_AMSRE('Cloud_Top_Pressure',SD_id,INFO); %int16 (>0)
%[tau,dimsizes_tau,tmp,tmp,success_tau]=get_hdf_data_AMSRE('Cloud_Optical_Thickness',SD_id,INFO); %int16 (>0)
%if success_tau==0
%    return
%end

%[qa1,dimsizes_QA_1km]=get_hdf_data_AMSRE('Cloud_Quality_Assurance',SD_id,INFO,1); %int8 (-128 min) unsigned
%NOTE the qa array here is the same as the 1km L2 qa array (but at lower
%res, according to the hdf var list from the MODIS webiste)
%[qa5,dimsizes_qa5]=get_hdf_data_AMSRE('Quality_Assurance_5km',SD_id,INFO,1); %int8 (-93 min) unsigned
%[mask1,dimsizes_mask1]=get_hdf_data_AMSRE('Cloud_Mask_1km',SD_id,INFO,1); %int8 unsigned
%[mask5,dimsizes_mask5]=get_hdf_data_AMSRE('Cloud_Mask',SD_id,INFO,1); %int8 unsigned
%[re,dimsizes_re]=get_hdf_data_AMSRE('Cloud_Effective_Radius',SD_id,INFO); %int16, (401 min) sf=0.01 (microns)
%Rob had this as Effective_Particle_Radius?
%[re_diff,dimsizes_re]=get_hdf_data_AMSRE('Cloud_Effective_Radius_Difference',SD_id,INFO); %int16, (401 min) sf=0.01 (microns)
%Cloud Particle Effective Radius two-channel retrieval using band 6 and
%band 20 differenced from band 7 retrieval and either band 1, 2, or 5 
%(specified in Quality_Assurance_1km).
%Plane 1 of SDS is band 6 - band 7 retrieval (Re_1.6 - Re_2.1)
%Plane 2 is band 20 - band 7 retrieval.      (Re_3.7 - Re_2.1)
%size(re_diff) = [1354        2030           2]
%band 6 = 1.6um, band 7 = 2.1um (the usual one), band 20 = 3.7um
%band 1 = 0.62, band 2 = 0.84, band 5 = 1.2 um - don't think we need to
%worry about whether is band 1, 2 or 5?

%Uncertainties - are these high when have dodgy tau values?
%NOTE - they are in % not absolute values!
%[re_un,dimsizes_re_un]=get_hdf_data_AMSRE('Cloud_Effective_Radius_Uncertainty',SD_id,INFO);
%N.B. - don't think there are any uncertainties for 1.6 and 3.7 um Re - not
%listed in L2 variable list.
%[tau_un,dimsizes_tau_un]=get_hdf_data_AMSRE('Cloud_Optical_Thickness_Uncertainty',SD_id,INFO); %int16 (>0)

%[scantime,dimsizes_scantime]=get_hdf_data_AMSRE('Scan_Start_Time',SD_id,INFO); %double
% (seconds since 1993-1-1 00:00:00.0) (!)
%[solar_zenith,dimsizes_solar_zenith]=get_hdf_data_AMSRE('Solar_Zenith_10km',SD_id,INFO); %int16, min=-32767
%[solar_azimuth,dimsizes_sensor_zenith]=get_hdf_data_AMSRE('Solar_Azimuth_10km',SD_id,INFO); %int16
%[sensor_zenith,dimsizes_sensor_zenith]=get_hdf_data_AMSRE('Viewing_Zenith_10km',SD_id,INFO); %int16
%[relative_azimuth,dimsizes_sensor_zenith]=get_hdf_data_AMSRE('Relative_Azimuth_10km',SD_id,INFO); %int16
%strange minima with solar and sensor ZA at the end of the swath
%NOTE - joint L2 files don't have solar and sensor azimuth, just the
%relative azimuth

%[surface_temp,dimsizes_temp]=get_hdf_data_AMSRE('Surface_Temperature',SD_id,INFO);

%[BTD,dimsizes_temp]=get_hdf_data_AMSRE('Brightness_Temperature_Difference',SD_id,INFO);
%[BT,dimsizes_temp]=get_hdf_data_AMSRE('Brightness_Temperature',SD_id,INFO);
%[rad_var,dimsizes_temp]=get_hdf_data_AMSRE('Radiance_Variance',SD_id,INFO); %for same bands as BT (see below)
%Bands 29,31,32,33,34,35,36
%Wavelength 8.4-8.7 (29), 10.8-11.3 (31), 11.8-12.3 (32), 13.2-13.5 (33),
%13.5-13.8 (34) 13.8-14.1 (35), 14.1-14.4 (36) um

%convert to Matlab date time (days since 01-Jan-0000)
%scantime_matlab = scantime/3600/24 + datenum('01-Jan-1993');
%then can do datestr(scantime_matlab(ilat,ilon),31) to get the date string

% ****  NOTE - NEED to flip the dimension of new variables to be consistent
% with those that are flipped below

%set the QA NaN flags back to zero as it causes problems with the bit stripping - QA
%values of zero indicate a bad retrieval anyway,
%so can then NaN the data at the QA stage. But keep the NaN fill values in
%the other arrays

%qa1(isnan(qa1))=0;  
%mask5(isnan(mask5))=0;  


% % qa5=flipdim(qa5,3);
% % qa1=flipdim(qa1,3);
% % %qa1=flipdim(qa1,1);
% % mask1=flipdim(mask1,3);
% % mask5=flipdim(mask5,2);
% % cf=flipdim(cf,2);
% % tau=flipdim(tau,2);  %0.1 (but sf=0.01?)
% % re=flipdim(re,2); %0.1 (but sf=0.01?)
% % re_diff=flipdim(re_diff,2); %
% % tau_un=flipdim(tau_un,2);  %
% % re_un=flipdim(re_un,2); %
% % 
% % lat=flipdim(lat,2);
% % lon=flipdim(lon,2);
% % 
% % t_top=flipdim(t_top,2);  %L?  %the L in Rob's IDL script means long integer needed to make sure that the values
% %    %can get high enough. Shouldn't be a problem in Matlab because I convert them all to double - check!
% % p_top=flipdim(p_top,2);
% % %need to alter the arrays as Rob was applying the offset and scale factors
% % %here, but I had already applied them in get_hdf_data_AMSRE.m
% % %Should probably flip the dimension there too - or remove the dimension
% % %flipping altogether
% % 
% % scantime = flipdim(scantime,2);
% % scantime_matlab = flipdim(scantime_matlab,2);
% % solar_zenith = flipdim(solar_zenith,2);
% % sensor_zenith = flipdim(sensor_zenith,2);


%[mask_5km,qapq_5km,qapp_5km] = flagread_5km_Dan(mask5);
%[qapq_5km,qapp_5km] = flagread_1km_Dan(qa1,1);
%qapq_1km = qapq_5km; %set these for consistency with L2
%qapp_1km = qapp_5km;

status = hdfsd('end',SD_id); 
%filtering_data_joint_L2

%convert the qa flags (info contained in bits) into real numbers
[qaAMSRE] = flagread_qaAMSRE(qa_L2amsre);

%Screen LWP (high res) for bad data
i=find(qaAMSRE(10,:,:)>0);
lwp_L2amsre(i)=NaN;

%Screen for detected sea ice
i2=find(qaAMSRE(2,:,:)==1);
lwp_L2amsre(i2)=NaN;

clear override_flags_openL2 %reset the flag if no errors
    %Read 1km mask and flag data
    
    if ~exist('suppress_output') | suppress_output==0
        disp('Done read L2 MODIS');
    end

catch openL2_error
   status = hdfsd('end',SD_id);
   clear override_flags_openL2 %reset the flag if there was an error
   rethrow(openL2_error); %re-issue the error
end





