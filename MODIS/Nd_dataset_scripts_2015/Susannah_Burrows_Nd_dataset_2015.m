% Will for now just convert the data sent to Dan McCoy from the .mat files
% to .nc

data_dir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Southern_Ocean_all_lons_no_water_path_confidence/CF_0.8_minCTT_173_ice_allowed_SZA_65/aqua/';

just_dates='no';
switch just_dates
    case 'no'

        files = dir([data_dir '*mockL3_saved_data_*_no_confidence_screening*.mat'])

        for i=1:length(files)
            year = str2num(files(i).name(19:22));
            
                [date_str,date_num_213] = date_from_day_of_year_func(213,year);
                [date_str,date_num_365] = date_from_day_of_year_func(365,year);
                [date_str,date_num_1] = date_from_day_of_year_func(1,year+1);
                [date_str,date_num_121] = date_from_day_of_year_func(121,year+1);

                times = [date_num_213:date_num_365 date_num_1:date_num_121];

                array_name = ['Days_since_01Jan0000'];
%                file_name = 'Dates_for_Susannah_Burrows_Jan2016.mat';

                eval_str=[array_name ' = times;']; eval(eval_str);
                
                save([data_dir files(i).name],array_name,'-APPEND');
            
            mat2nc_Dan([data_dir files(i).name],[data_dir files(i).name '.nc']);
        end

end

% In Jan 2016 she asked about the dates - realised that had just been using
% the day of the year without regard for the leap year. So, since the data
% in file Y2007 runs from Aug 2007 to Apr 2008 this will be affected along
% with the 2008 file. The code below produces a time axis for her using
% Matlab time convention (days since 01-Jan-0000)
% Data ran from Y2006 to Y2013 and used days [213:365] for Y and [1:121]
% for Y+1 in each array.

% years=[2006:2013];
% for iyear=1:length(years)
%     year = years(iyear);
% 
%     [date_str,date_num_213] = date_from_day_of_year_func(213,year);
%     [date_str,date_num_365] = date_from_day_of_year_func(365,year);
%     [date_str,date_num_1] = date_from_day_of_year_func(1,year+1);
%     [date_str,date_num_121] = date_from_day_of_year_func(121,year+1);
% 
%     times = [date_num_213:date_num_365 date_num_1:date_num_121];
% 
%     array_name = ['Y' num2str(year) '_days_since_01Jan0000'];
%     file_name = 'Dates_for_Susannah_Burrows_Jan2016.mat';
% 
%     eval_str=[array_name ' = times;']; eval(eval_str);
% 
%     if iyear==1
%         save([data_dir file_name],array_name);
%     else
%         save([data_dir file_name],array_name,'-APPEND');
%     end
% 
% end
% 
% mat2nc_Dan([data_dir file_name],[data_dir file_name '.nc']);
