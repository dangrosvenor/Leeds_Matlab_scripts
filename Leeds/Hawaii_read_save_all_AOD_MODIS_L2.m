function Hawaii_read_save_all_AOD_MODIS_L2()

filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/'    
%file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';

save_name = 'Hawaii_all_L2_04_data.mat';

var_list={'lat','lon','AOD_550_Dark_Target_Deep_Blue_Combined','AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag','scantime','scantime_matlab',...
    'solar_zenith','solar_azimuth','sensor_zenith','sensor_azimuth','modis_date_time','iday'};


files=dir([filedir 'M*D04*.hdf']);

count=0
for i=1:1:length(files)
    count=count+1;
    %i
    file_L2_04{i} = files(i).name;
    file_name_h5 = files(i).name
    imodis_file_override=1;
    open_L2_C6_04_AOD_MODIS_file_01 %run the read function
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