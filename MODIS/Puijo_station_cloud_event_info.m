CDP_datatype = 'New - 5 min averages';
%CDP_datatype = 'Old hourly averages';


switch CDP_datatype
    case 'New - 5 min averages'
        filedir_vocals = '/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/';
        filename_vocals = 'CDP_5min_MT_D_Nd_LWC_stdNd.mat';
        %
        %         Name                   Size             Bytes  Class     Attributes
        %
        %   Final_5min_CDP      8896x10            711680  double

        %Columns
        %  1  = Matlab Date
        %  2 - 7 = year, month, day, hour, min, sec
        %  8  = Nd
        %  9  = LWC
        % 10  = std_dev


        load([filedir_vocals filename_vocals]);



        switch assume_UTC
            case 'no'
                MatlabTime_Puijo = Final_5min_CDP(:,1) - 2/24 - 2.5/60/24; %CDP data is in UTC+2, so convert to UTC - for this datafile we should
                %have the mid-point (in time) of the 5 mins of the average. Have subtracted
                %2.5 mins since then the time of the CDP average will correspond to the
                %MODIS swath where the centre of the swath was sampled at teh same time as
                %the CDP time mid-point. E.g. CDP measurement at 10:02:30 and MODIS swath
                %at 10:00 (which runs from 10:00 to 10:05).

            case 'yes'  %here assume that the CDP is in UTC (not UTC + 2) for a test
                MatlabTime_Puijo = Final_5min_CDP(:,1) - 0/24 - 2.5/60/24;
                
            case 'UTC_minus_10mins'  %here assume that the CDP is in UTC (not UTC + 2) for a test
                MatlabTime_Puijo = Final_5min_CDP(:,1) - 0/24 - 2.5/60/24 + 10/60/24;                

        end



        Nd_Puijo = Final_5min_CDP(:,8);
        LWC_Puijo = Final_5min_CDP(:,9);
        stdNd_Puijo = Final_5min_CDP(:,10);

    case 'Old hourly averages'

        filedir_vocals = '/home/disk/eos1/d.grosvenor/modis_work/Irshad_data/';
        filename_vocals = 'Corrected CDP 2006-2012_CSV.csv';

        %fid_voc=fopen(filename_vocals,'rt');
        % fgetl(fid_voc)
        %[dat_VOCALS]=fscanf(fid_voc,'%f',[6 Inf]);
        [MatlabTime_Puijo,temp,temp,temp,temp,Nd_Puijo] = textread([filedir_vocals filename_vocals],'','delimiter',',');

        MatlabTime_Puijo = MatlabTime_Puijo - 1.5/24; %CDP data is in UTC+2, so convert to UTC and add 30 mins
        %because the data is averaged between ay 01:00 and 02:00 - i.e. will use the midpoint

        %[month_voc,day_voc,hours_voc,mins_voc,secs_voc,lat_voc,lon_voc]=read_vocals_info_files_func(filedir_vocals,filename_vocals);
        %year_voc=2008*ones(size(month_voc));

end