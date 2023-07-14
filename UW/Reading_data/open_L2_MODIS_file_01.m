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
[cf,dimsizes_cf]=get_hdf_data_dan('Cloud_Fraction',SD_id,INFO); %int8 (0...127)
[t_top,dimsizes_ctt]=get_hdf_data_dan('Cloud_Top_Temperature',SD_id,INFO); %int16 (>0)
[p_top,dimsizes_ctp]=get_hdf_data_dan('Cloud_Top_Pressure',SD_id,INFO); %int16 (>0)
[tau,dimsizes_tau]=get_hdf_data_dan('Cloud_Optical_Thickness',SD_id,INFO); %int16 (>0)
[qa1,dimsizes_QA_1km]=get_hdf_data_dan('Quality_Assurance_1km',SD_id,INFO,1); %int8 (-128 min) unsigned
[qa5,dimsizes_qa5]=get_hdf_data_dan('Quality_Assurance_5km',SD_id,INFO,1); %int8 (-93 min) unsigned
[mask1,dimsizes_mask1]=get_hdf_data_dan('Cloud_Mask_1km',SD_id,INFO,1); %int8 unsigned
[mask5,dimsizes_mask5]=get_hdf_data_dan('Cloud_Mask_5km',SD_id,INFO,1); %int8 unsigned
[re,dimsizes_re]=get_hdf_data_dan('Cloud_Effective_Radius',SD_id,INFO); %int16, (401 min) sf=0.01 (microns)
%Rob had this as Effective_Particle_Radius?
[re_diff,dimsizes_re]=get_hdf_data_dan('Effective_Radius_Difference',SD_id,INFO); %int16, (401 min) sf=0.01 (microns)
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
[re_un,dimsizes_re_un]=get_hdf_data_dan('Cloud_Effective_Radius_Uncertainty',SD_id,INFO);
%N.B. - don't think there are any uncertainties for 1.6 and 3.7 um Re - not
%listed in L2 variable list.
[tau_un,dimsizes_tau_un]=get_hdf_data_dan('Cloud_Optical_Thickness_Uncertainty',SD_id,INFO); %int16 (>0)

[scantime,dimsizes_scantime]=get_hdf_data_dan('Scan_Start_Time',SD_id,INFO); %double
% (seconds since 1993-1-1 00:00:00.0) (!)
[solar_zenith,dimsizes_solar_zenith]=get_hdf_data_dan('Solar_Zenith',SD_id,INFO); %int16, min=-32767
[solar_azimuth,dimsizes_sensor_zenith]=get_hdf_data_dan('Solar_Azimuth',SD_id,INFO); %int16
[sensor_zenith,dimsizes_sensor_zenith]=get_hdf_data_dan('Sensor_Zenith',SD_id,INFO); %int16
[sensor_azimuth,dimsizes_sensor_zenith]=get_hdf_data_dan('Sensor_Azimuth',SD_id,INFO); %int16
%strange minima with solar and sensor ZA at the end of the swath

[surface_temp,dimsizes_temp]=get_hdf_data_dan('Surface_Temperature',SD_id,INFO);

[BTD,dimsizes_temp]=get_hdf_data_dan('Brightness_Temperature_Difference',SD_id,INFO);
[BT,dimsizes_temp]=get_hdf_data_dan('Brightness_Temperature',SD_id,INFO);
[rad_var,dimsizes_temp]=get_hdf_data_dan('Radiance_Variance',SD_id,INFO); %for same bands as BT (see below)
%Bands 29,31,32,33,34,35,36
%Wavelength 8.4-8.7 (29), 10.8-11.3 (31), 11.8-12.3 (32), 13.2-13.5 (33),
%13.5-13.8 (34) 13.8-14.1 (35), 14.1-14.4 (36) um

%convert to Matlab date time (days since 01-Jan-0000)
scantime_matlab = scantime/3600/24 + datenum('01-Jan-1993');
%then can do datestr(scantime_matlab(ilat,ilon),31) to get the date string

% ****  NOTE - NEED to flip the dimension of new variables to be consistent
% with those that are flipped below

%set the QA NaN flags back to zero as it causes problems with the bit stripping - QA
%values of zero indicate a bad retrieval anyway,
%so can then NaN the data at the QA stage. But keep the NaN fill values in
%the other arrays

qa1(isnan(qa1))=0;  
qa5(isnan(qa5))=0;  
mask1(isnan(mask1))=0;  
mask5(isnan(mask5))=0;  


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
% % %here, but I had already applied them in get_hdf_data_dan.m
% % %Should probably flip the dimension there too - or remove the dimension
% % %flipping altogether
% % 
% % scantime = flipdim(scantime,2);
% % scantime_matlab = flipdim(scantime_matlab,2);
% % solar_zenith = flipdim(solar_zenith,2);
% % sensor_zenith = flipdim(sensor_zenith,2);

%if s_dims(0) eq 1354 then begin
make_subset=0;
if make_subset==1
    if size(re,1)==1354
        %FULL HDF FILE case
        mask1=mask1(:,160:1189,:);
        qa1=qa1(:,160:1189,:);
        mask5=mask5(33:238,:);
        qa5=qa5(:,33:238,:);
        re=re(160:1189,:);
        tau=tau(160:1189,:);

        if strcmp(region_flag,'epic')==1
            lwp=flipdim(lwp,2);
            lwp=lwp(160:1189,:);
        end


        lon=lon(33:238,:);
        lat=lat(33:238,:);
        cf=cf(33:238,:);

        t_top=t_top(33:238,:);
        p_top=p_top(33:238,:);
    end
end

[mask_5km,qapq_5km,qapp_5km] = flagread_5km_Dan(mask5,qa5);
[qapq_1km,qapp_1km] = flagread_1km_Dan(qa1,1);
%qapp_1km is a continuation of the 1km QA as documented in the
%MODIS_Quality_Assurance_etc....pdf
%plan document. Starts at the 3rd byte ("Primary Cloud Retrieval Phase
%Flag").
[mask_1km] = maskread_1km_Dan(mask1);

% --- Do some screening based on the confidence of the retrieval ----

%  mask(1,:,:) - Cloud Mask Status Flag  
% 0 = Undetermined, 1 = Determined. Although prob not any undetermined cloud points with good confidence 

%  mask(2,:,:) - Cloud Mask Confidence  0: confident cloudy

%  qapq_1km(7) : WATER PATH CONFIDENCE
%                    0: bad% 1 marginal% 2 good% 3 very good

if ~exist('i_apply_water_path_confidence') | i_apply_water_path_confidence==1
    %Only use points where the water path confidence is very good
    % Was found that this restriction (when combined with a restriction to
    % liquid only clouds) causes the re2.1 values to be truncated at 20um,
    % which is likely to have a large impact.
    ind=find( mask_1km(1,:,:) < 1 | mask_1km(2,:,:) > 0 | qapq_1km(7,:,:) < 3);       
else
    %Use all cases except bad
     ind=find( mask_1km(1,:,:) < 1 | mask_1km(2,:,:) > 0 | qapq_1km(7,:,:) < 1);      
end

    %mask_1km(2,:,:) has changed cf Rob's script - the (2,:,:) here is actually a
    %(1,:,:) one in Rob's script due to the zero-based indices quoted there
    %(checked in the .pdf below, 15th Feb, 2012).
    % Also, 0 now means good! See
    %http://modis-atmos.gsfc.nasa.gov/_docs/QA_Plan_2011_01_26.pdf -
 
%  mask_1km(3)  :  SUN GLINT FLAG  - so is actually mask_1km(4,:,:)
%			0: yes% 1: no

%  **** N.B. - need to add one to these indices for Matlab 1-based indices
%  ************************************************************************
%===========================================================================
%  QUALITY ASSURANCE FLAGS for 1km DATA - PRODUCT QUALITY (BYTES 0 and 1)
%  INFORMATION: (http://modis-atmos.gsfc.nasa.gov/reference_atbd.html
%
%  qapq_1km(0) : OPTICAL THICKNESS RETRIEVAL 
%                    0: not useful% 1 useful
%  qapq_1km(1) : OPTICAL THICKNESS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_1km(2) : OPTICAL THICKNESS OUT OF BOUNDS 
%                    0: within bounds (tau<150)
%                    1: (100<tau<150)
%                    2: (tau>150)
%                    3: surface reflectance too large
%  qapq_1km(3) : EFFECTIVE RADIUS GENERAL QA
%		     0: not useful
% 		     1: useful
%  qapq_1km(4) : EFFECTIVE RADIUS CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good			
%  qapq_1km(5) : WATER PATH GENERAL QA
%		     0: not useful
% 		     1: useful
%  qapq_1km(6) : WATER PATH CONFIDENCE 
%                    0: bad% 1 marginal% 2 good% 3 very good
%  qapq_1km(7) : CLOUD PHASE DETERMINATION
%		     0: SWIR algorithm not run
%		     1: CLEAR		
%		     2: WATER
%		     3: ICE
%		     4: MIXED PHASE or UNDETERMINED	

%  qapq_1km(8) : 1621 cloud retrieval outcome
%		     0: not attempted or unsuccessful
%		     1: successful		

%  
%--------------------------------------------------------------------------
%-
    
    
    %  lwp(ind)=0.00
    %if are not confident cloudy then we NaN the cloud variables like re
    %and tau
    re(ind)=NaN;   %set to zero or to NaN??
    
    re_diff_plane1 = re_diff(:,:,1);
    re_diff_plane1(ind)=NaN;
    re_diff_plane2 = re_diff(:,:,2);
    re_diff_plane2(ind)=NaN;
    re_diff(:,:,1)=re_diff_plane1;  %the two planes refer to different bands (see above)
    re_diff(:,:,2)=re_diff_plane2;
        
%    W=5/9*rhoL.*tau.*reff; %this is the one that takes into account
%    adiabaticity - kg/m2  = 0.555*1e3*tau*reff
    
%5/9 = 0.555...
%    lwp=1.0e3*555.5*1.0e-6.*re.*tau; %from Rob - so is in g per m2 - same
%    as below
lwp = 5/9*1000*tau.*re*1.0e-6 *1e3; %the last 1e3 is to convert to g/m2

    lwp(ind)=NaN;
    tau(ind)=NaN;
    
    re_un(ind)=NaN;
    tau_un(ind)=NaN;



%sub-sample to 5km
%need to check that these are the correct locations for 5 km subsampling
%Perhaps can do this with the L2 joint data?
re_5km = re(5:5:end,5:5:end);
tau_5km = tau(5:5:end,5:5:end);
lwp_5km = lwp(5:5:end,5:5:end);
re_un_5km = re_un(5:5:end,5:5:end);
tau_un_5km = tau_un(5:5:end,5:5:end);

%   segsize=256L
%   overlap=128L

   segsize=256;
   overlap=128;


   % This sorted out a problem with the different sensors - think this
   % should have been removed for later collections, but should check by
   % plotting
%   remove_lines,lwp,npts,offset,re
%   remove_lines_block,lwp,re

%   ind=where(re le 0.0)
   ind2=find(re < 0.0);
%   tau=1.5*lwp./re;
   tau(ind2)=0.0;



   s_dims=size(lwp);
   s_dims_5km=size(lon);
   
   
%This part is to put the swath data into a grid where each point represents
%an area that is the same size (1x1km). So if the pixel is larger (as they are at
%the swath edges) then it is split into N 1x1 km areas.
%We won't use this for now (Dan).
%    restore,'~/MODIS/modis_pixel_size.idl' 
%    restore,'~/MODIS/modis_pixel_size_5km.idl'
%    nxn=1280;
%    nyn=1920;
%    nxn_5km=256;
%    nyn_5km=384;
% 
%    newpix=newpix+s_dims(1)/2;
%    newpix5km=newpix5km+s_dims_5km(1)/2;
% 
% %    lwp=lwp(newpix,*)
% %    tau=tau(newpix,*)
% %    re=re(newpix,*)
% % 
% %    mask_1km=mask_1km(*,newpix,*)
% %    qapq_1km=qapq_1km(*,newpix,*)
% %    qapp_5km=qapp_5km(*,newpix5km,*)
% %    qapq_5km=qapq_5km(*,newpix5km,*)
% % 
% % 
% %    cf=cf(newpix5km,*)
% %    lat=lat(newpix5km,*)
% %    lon=lon(newpix5km,*)
% %    p_top=p_top(newpix5km,*)
% %    t_top=t_top(newpix5km,*)
%    
%    
%    lwp=lwp(newpix,:);
%    tau=tau(newpix,:);
%    re=re(newpix,*:);
% 
%    mask_1km=mask_1km(:,newpix,:);
%    qapq_1km=qapq_1km(:,newpix,:);
%    qapp_5km=qapp_5km(:,newpix5km,:);
%    qapq_5km=qapq_5km(:,newpix5km,:);
% 
% 
%    cf=cf(newpix5km,:);
%    lat=lat(newpix5km,:);
%    lon=lon(newpix5km,:);
%    p_top=p_top(newpix5km,:);
%    t_top=t_top(newpix5km,:);

   

catch openL2_error
   clear override_flags_openL2 %reset the flag if there was an error
   rethrow(openL2_error); %re-issue the error
end

clear override_flags_openL2 %reset the flag if no errors


filtering_data_L2

    %Read 1km mask and flag data
disp('Done read L2 MODIS');




% Quality Assurance 1km reports on Cloud Optical Properties                           
% algorithm performance.  Refer to MOD_PR06OD User Documentation and the              
% MODIS atmosphere QA plan for complete descriptions and coding examples.             
%                                                                                     
%  Bit fields within each byte are numbered from the left:                            
%  7, 6, 5, 4, 3, 2, 1, 0.                                                            
%  The left-most bit (bit 7) is the most significant bit.                             
%  The right-most bit (bit 0) is the least significant bit.                           
%                                                                                     
%                                                                                     
%  Byte 0                                                                             
%  ------                                                                             
%                                                                                     
%  bit field       Description                             Key                        
%  ---------       -----------                             ---                        
%                                                                                     
%  Byte 0 -----------------------------------------------------------------           
%   0              Optical Thickness General QA       0 = Not Useful                  
%                                                     1 = Useful                      
%   2,1            Optical Thickness Confidence QA    00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%   4,3            Optical Thickness out-of-bounds    00 = OT < 100                   
%                                                     01 = 100 < OT < 150             
%                                                     10 = OT > 150                   
%                                                     11 = Albedo  too high           
%   5              Effective Radius General QA        0 = Not Useful                  
%                                                     1 = Useful                      
%   7,6            Effective Radius Confidence QA     00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%  Byte 1 -----------------------------------------------------------------           
%   0              Liquid Water Path General QA       0 = Not Useful                  
%                                                     1 = Useful                      
%   2,1            Liquid Water Path Confidence QA    00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%   5,4,3          1621 Retrieval processing path     000 = No Cloud Mask             
%                                                     001 = No Cloud                  
%                                                     010 = Water Cloud               
%                                                     011 = Ice Cloud                 
%                                                     100 = Unknown Cloud             
%   6              1621 Retrieval Outcome             0 = Failed/No  attempt          
%                                                     1 = Successful                  
%  Byte 2 -----------------------------------------------------------------           
%   2,1,0          Primary retrieval processing path  000 = No Cloud Mask             
%                                                     001 = No Cloud                  
%                                                     010 = Water Cloud               
%                                                     011 = Ice Cloud                 
%                                                     100 = Unknown Cloud             
%   3              Retrieval Outcome                  0 = Failed/No  attempt          
%                                                     1 = Successful                  
%   4              Rayleigh Correction                0 = No Correction               
%                                                     1 = Correction                  
%   5              Water Vapor Correction             0 = No Correction               
%                                                     1 = Correction                  
%   7,6            Band Used for Optical Thickness Retrieval                          
%                                                     00 = No attempt                 
%                                                     01 = .645 micron                
%                                                     10 = .858 micron                
%                                                     11 = 1.24 micron                
%  Byte 3 -----------------------------------------------------------------           
%   0              Optical Thickness 1621 General QA  0 = Not Useful                  
%                                                     1 = Useful                      
%   2,1            Optical Thickness 1621 Condifence QA                               
%                                                     00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%   3              Effective Radius 1621 General QA   0 = Not Useful                  
%                                                     1 = Useful                      
%   5,4            Effective Radius 1621 Confidence QA                                
%                                                     00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%   6,7            Clear Sky Restoral Type QA                                         
%                                        00 = Not Restored                            
%                                        01 = Restored Via Edge detection             
%                                        10 = Restored Via Spatial  Variance          
%                                        11 = Restored Via 250m Tests                 
%  Byte 4 -----------------------------------------------------------------           
%   0              Water Path 1621 General QA         0 = Not Useful                  
%                                                     1 = Useful                      
%   2,1            Water Path 1621 Confidence QA      00 = No confidence              
%                                                     01 = Marginal                   
%                                                     10 = Good                       
%                                                     11 = Very Good                  
%   5,4,3          Multi Layer Cloud Flag      000 = Cloud Mask Undet                 
%                                              001 = Decision tree stop               
%                                              010 = single layer: water              
%                                              011 = multi layer: water               
%                                              100 = single layer: ice                
%                                              101 = multi layer: ice                 
%                                              110 = single layer:  unknown           
%                                              111 = multi layer: unknown             
% 
