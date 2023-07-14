function Hawaii_read_save_all_Nd_MODIS_L2()

filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/MOD06_L2/'    
%file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';

save_name = 'Hawaii_all_L2_06_data.mat';

%List of variables to save for all times (plus, all of the data for the last swath will be saved)
var_list={'lat','lon','Plon2_L2','Plat2_L2','N37','re37','tau37','W37','N37_pcl','re37_pcl','tau37_pcl','W37_pcl','cf','cth_1km'...
    ,'cloud_phase'... %,'mask_5km','qapq_5km','qapp_5km','mask_1km','qapq_1km','qapp_1km'...
    ,'cloud_phase_flag','cloud_mask_status_flag','cloud_mask_confidence_flag'...
    ,'scantime','scantime_matlab','solar_zenith','solar_azimuth','sensor_zenith','sensor_azimuth','modis_date_time','iday'};

%cloud_phase % 0 = cloud mask unavailable, missing data,etc. (no phase result), 1= no phase result due to clear sky, etc.; 
  %2=liquid water; 3=ice; 4=undetermined
  % Is NaN for fill values where there was no daylight, etc. So will want
  % to ignore those. And also ignore =0 values where the cloud mask was
  % unavailable.
  % Can work out a liquid cloud fraction as Nliq/(Nliq+Nice+Nundet+Nclear)                                                                                                                                                                                                                                                                                                                      


files=dir([filedir 'M*D06*.hdf']);

count=0
for i=1:length(files)
    count=count+1;
    %i
    file_L2_06{i} = files(i).name;
    file_name_h5 = files(i).name
    imodis_file_override=1;
    open_L2_C6_MODIS_file_01    
     %run the read function
    %if count>1 & (size(lat,1) ~= size(lat_old,1) | size(lat,2) ~= size(lat_old,2))
    %    siz = size(lat); siz2 = size(lat_old);
    %    fprintf(1,['\n***Error - size(lat)=[' num2str(siz(1)) ' ' num2str(siz(2)) '] , but size(lat_old)=[' num2str(siz2(1)) ' ' num2str(siz2(2)) ']' ]);
    %end
    %lat_old = lat;
    for ivar=1:length(var_list)
       var = var_list{ivar};
       eval_str = [var '_all_times{i} = ' var ';']; eval(eval_str);
    end

end

save([filedir save_name],'-V7.3');