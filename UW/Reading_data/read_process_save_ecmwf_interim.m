%Process ECMWF .nc files into .mat files
%Read in ECMWF LTS, etc from Manuels's data (converted to NetCDF)
%and calculates daily mean temps and monthly temps for each hour of the day
%(6 hourly)
% Manuel's data is in this folder:-
% /home/disk/cresta2/mzuluaga/Data/ERA_interim/pressure/
% and then in the q_press and t_press folders
% But I copied it to /home/disk/eos5/d.grosvenor/ERA_Interim/ManuelZ/
%This calls   -- read_process_save_ecmwf_interim_func --   to do the processing
%       Use   --  load_saved_ecmwf_interim --    to read in the .mat files

%adds this to the filename of the .mat files


keyword_EC='daily_LTS_'; %for LTS calculated from the daily (instantaneous) temperatures


Man_dir = '/home/disk/eos5/d.grosvenor/ERA_Interim/ManuelZ/';

years = [1998:2012];
years = [2008];


var_strs = {'t_press','q_press'};
var_strs = {'t_press'};

for ivar=1:length(var_strs)

    var_str = var_strs{ivar};


    for iyear=1:length(years)
        iyear
        year_str = num2str(years(iyear));

        files = dir([Man_dir var_str '*' year_str '.nc']);
        nc_file=[Man_dir files(1).name];
        save_file = [nc_file keyword_EC '.mat'];


        switch var_str
            case 'q_press'
                var_str_ec = 'Q_GDS0_ISBL';
                [ecqv_daily_Man,time_daily,ecqv_mon_Man,ecLat_Man,ecLon_Man,levs_Man_L6,ecLTS_mon_Man] = read_process_save_ecmwf_interim_func(nc_file,var_str_ec);
                save(save_file,'ecqv_daily_Man','time_daily','ecqv_mon_Man','ecLat_Man','ecLon_Man','levs_Man_L6','-V7.3');
            case 't_press'
                var_str_ec = 'T_GDS0_ISBL';
                [ecT_daily_Man,time_daily,ecT_mon_Man,ecLat_Man,ecLon_Man,levs_Man_L6,ecLTS_mon_Man] = read_process_save_ecmwf_interim_func(nc_file,var_str_ec);
                %N.B. ecLTS_mon_Man is now the monthly LTS calculated from
                %daily LTS values (not from monthly mean T values)
                save(save_file,'ecLTS_mon_Man','ecT_daily_Man','time_daily','ecT_mon_Man','ecLat_Man','ecLon_Man','levs_Man_L6','-V7.3');
        end





    end


end







