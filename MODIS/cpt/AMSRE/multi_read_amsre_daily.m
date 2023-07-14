%read in monthly amsre data, or process daily files into monthly averages
%/home/disk/eos5/d.grosvenor/AMSRE
%AMSRE is onboard AQUA!
% Also have SSMI (several satellites f13, f15, f16, f17) and Windsat in
% /home/disk/eos5/d.grosvenor/SSMI      /f13   /windsat, etc.

% use save_amsre_loaded_data for saving
% This script just loads it all in without saving. Woudl be useful to have
% sst files for each year of L3 that can be automatially loaded in
% load_modis_saved_vars


try

    clear years_select months_select years months

    if ~exist('ioverride_read_amsre') | ioverride_read_amsre==0
        
        sat = 'AMSRE';

        clear chosen_files

        amsre_action='process'; %process the daily fields into monthly averages
        amsre_action='read'; %read in a selection of monthly averages
        amsre_action='read_daily'; %read in a selection of daily files and store in a big array
        %amsre_action='read_saved'; %read in a specific saved dataset


        case_study='MODIS_L2 comparisons';
        case_study='CPT'; %is all data 2006-2010 (but with Dec 2005 and Jan 2011)
        %case_study='2008'; %just 2008
        %case_study='2007-2008 DanMcCoy';
%        case_study='2002-2005'; %Covering the time period not covered by the CPT dataset - would putting it all in one file be too much?
        %case_study = 'Arctic Box MODIS_L2 comparisons';
%        case_study = 'One month'; %One month of data to retrieve, e.g. lat lon
%        case_study = 'Choose days'; %Just seleect one day specified below in chosen_files:-
%         case_study = 'ACSIS global 2009-2010';

        iday=1;
        %chosen_files(iday).name = 'amsre_20081025v7'; year_single = 'y2008'; month_single = 'm10';  iday=iday+1; %26th Oct, 2008 for UM POC study - select 'Choose days' above
        %chosen_files(iday).name = 'amsre_20081026v7'; year_single = 'y2008'; month_single = 'm10';  iday=iday+1;%26th Oct, 2008 for UM POC study - select 'Choose days' above
        %chosen_files(iday).name = 'amsre_20081027v7'; year_single = 'y2008'; month_single = 'm10'; iday=iday+1; %26th Oct, 2008 for UM POC study - select 'Choose days' above
        %chosen_files(iday).name = 'amsre_20081115v7'; year_single = 'y2008'; month_single = 'm11';  iday=iday+1;%12th Nov, 2008 for UM Boutle case - select 'Choose days' above
        %chosen_files(iday).name = 'amsre_20081114v7'; year_single = 'y2008'; month_single = 'm11';  iday=iday+1;%12th Nov, 2008 for UM Boutle case - select 'Choose days' above
        chosen_files(iday).name = 'amsre_20081113v7'; year_single = 'y2008'; month_single = 'm11';  iday=iday+1;%12th Nov, 2008 for UM Boutle case - select 'Choose days' above
        %chosen_files(iday).name = 'amsre_20081112v7'; year_single = 'y2008'; month_single = 'm11';  iday=iday+1;%12th Nov, 2008 for UM Boutle case - select 'Choose days' above


        %% Location of the AMSRE stored files

        %filedir = '/home/disk/eos10/robwood/AMSR/bmaps_v05/y2007/m01/';
        filedir = '/home/disk/eos5/d.grosvenor/AMSRE/';
        % %file_name2 = 'amsre_200701v5.gz';
        % file_name2 = 'amsre_200701v7'; %these are monthly values - should be 6 maps - however they've averaged ascending and descending
        % file_name2 = 'amsre_20070131v5'; %thess are the daily maps (should be 7 maps - the extra one is time UTC,
        % %which is the first one)
        % %file_name2 = 'amsre_20070131v5_d3d'; %d3d are the 3 day averages
        % %file_name2 = 'amsre_20111001v7'; %weekly files - same as monthly (i.e. have also averaged day and night)
        % file_name2 = 'amsre_20020612v7';
        
        var_names={'sst','lwp','time','rain','wspdLF','wspdMF','vapor'};  %time is for each location
        
             %   [mingmt,sst,wspdLF,wspdMF,vapor,cloud,rain,wspdAW,wdir]
                                        %   mingmt is gmt time in hours
                                        %   sst  is surface water temperature at depth of about 1 mm in deg C
                                        %   wspdLF is 10 meter surface wind in m/s made using 10.7 GHz channel and above, low frequency channels
                                        %   wspdMF is 10 meter surface wind in m/s made using 18.7 GHz channel and above, medium frequency channels
                                        %   vapor is atmospheric water vapor in millimeters
                                        %   cloud is liquid cloud water in millimeters
                                        %   rain  is rain rate in millimeters/hour
                                        %   wspdAW is 10 meter surface wind for all weather conditions made using 3 algorithms
                                        %   wdir is wind direction
                                        %   oceanographic convention, blowing North = 0 in degrees
% HOWEVER, only these ones are output from the read function                                        
%[amsre_time,sst,amsre_wspdLF,amsre_wspdMF,amsre_vapor,amsre_lwp,amsre_rain]                                        


    end



    gcm_str = 'AMSRE';
    gcm_str_select = gcm_str;

    switch case_study
        case 'Choose days'
            Nsub_av=1;
            ireduce_to_1deg=0;  %Flag to tell it to do this - have changed to zero for case_study = 'Choose days' option.
        otherwise
            ireduce_to_1deg=1;
            Nsub_av = 4; %average over 4x4 blocks to make 1x1 degree for 'read' action
    end

    time_series_type = 'AMSRE';


    %% Select what to read
    switch case_study
        case 'MODIS_L2 comparisons'
            years_select = {'y2006','y2007'}; months_select = { {'m11', 'm12'},{'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'} };
            %months_select = { {'m12'},{'m01','m02','m06','m07','m08'} };
            %         years_select = {'y2007'}; months_select = { {'m08','m09','m10','m11'} };

            %         years_select = {'y2007'}; months_select = { {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11'} };

        case 'Arctic Box MODIS_L2 comparisons'
            years_select = {'y2007','y2008','y2009','y2010'};
            months_select = { {'m06'},{'m06'},{'m06'},{'m06'} };

        case 'CPT'
            years_select = {'y2005','y2006','y2007','y2008','y2009','y2010','y2011'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m12'},all_months,all_months,all_months,all_months,all_months,{'m01'}};

        case '2002-2005'
            years_select = {'y2002','y2003','y2004','y2005','y2006'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m06','m07','m08','m09','m10','m11','m12'},all_months,all_months,all_months,{'m01'}};
            %Only have from June for 2002 since this is when the satellite was
            %launched. So will prob miss a few days in June 2002.

        case {'One month','Choose days'}
            %        years_select = {'y2008'};
            %        all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            %        months_select = {{'m10'}};

            years_select = {year_single};
            months_select = {{month_single}};

        case '2008'
            %load some data either side to pad for the smoothing
            years_select = {'y2007','y2008','y2009'};
            %        years_select = {'y2007'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m12'},all_months,{'m01'}};
            %        months_select = {{'m12'}};

        case '2005'
            %load some data either side to pad for the smoothing
            years_select = {'y2004','y2005','y2006'};
            %        years_select = {'y2007'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m12'},all_months,{'m01'}};
            %        months_select = {{'m12'}};


        case '2007-2008 DanMcCoy'
            %load some data either side to pad for the smoothing
            years_select = {'y2006','y2007','y2008','y2009'};
            %        years_select = {'y2007'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m12'},all_months,all_months,{'m01'}};
            %        months_select = {{'m12'}};

        case 'ACSIS global 2009-2010'
            years_select = {'y2009','y2010'};
            all_months = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'}
            months_select = {{'m04','m05','m06','m07','m08','m09','m10','m11','m12'},{'m01','m02','m03'}};



    end




    if strcmp(months_select{1}{1},'all')==1
        for iyear=1:length(years_select)
            months_select{iyear} = {'m01','m02','m03','m04','m05','m06','m07','m08','m09','m10','m11','m12'};
        end
    end


    for iyear=1:length(years_select)
        years_requested(iyear) = str2num(years_select{iyear}(2:end));
    end


    years_amsre_str='';
    if max(diff(years_requested)==1)
        years_amsre_str=[num2str(years_requested(1)) ' to ' num2str(years_requested(end))];
    else
        for i=1:length(years_requested)
            years_amsre_str=[years_amsre_str ' ' num2str(years_requested(i))];
        end
    end

    gcm_years_loaded_str = years_amsre_str;

    years_required_for_mean = years_requested;



    %these are the cell centres
    ILAT=[1:720];
    ILON=[1:1440];

    XLAT=0.25*ILAT-90.125;
    XLON=0.25*ILON-0.125;

    XLAT_edges=[0.25*ILAT-90.25 90];
    XLON_edges=[0.25*ILON-0.25 360];

    if ireduce_to_1deg

        sub_lat = [1:4:720];
        sub_lon = [1:4:1440];

        sub_lat2 = [1:4:721];
        sub_lon2 = [1:4:1441];

        %the new centres are halfway between the old edges
        lat_centres_sub = 0.5 * (XLAT_edges(sub_lat2(1:end-1)) + XLAT_edges(sub_lat2(2:end)));
        lon_centres_sub = 0.5 * (XLON_edges(sub_lon2(1:end-1)) + XLON_edges(sub_lon2(2:end)));

        lat_edges_sub = [lat_centres_sub-0.5 90];
        lon_edges_sub = [lon_centres_sub-0.5 360];

    else
        lat_centres_sub = XLAT;
        lat_edges_sub = XLAT_edges;
        lon_centres_sub = XLON;
        lon_edges_sub = XLON_edges;

        sub_lat = lat_centres_sub;
        sub_lon = lon_centres_sub;

    end

    nlon = size(lon_centres_sub,2);
    dlon = lon_centres_sub(2) - lon_centres_sub(1);

    %will manipulate the data so that 0deg is in the centre, so do the same here
    lon_centres_sub2(1:nlon/2) = lon_centres_sub(nlon/2+1:nlon);
    lon_centres_sub2(nlon/2+1:nlon) = lon_centres_sub(1:nlon/2);
    i180 = find(lon_centres_sub2>180);
    lon_centres_sub2(i180) = lon_centres_sub2(i180) - 360;

    lon_edges_sub2 = [lon_centres_sub2-dlon/2 180];

    %flip the latitudes to be consistent with MODIS (12th Feb, 2013)
    lat_edges_sub = flipdim(lat_edges_sub,2);
    lat_centres_sub = flipdim(lat_centres_sub,2);

    [gcm_Plon2D_edges_AMSRE,gcm_Plat2D_edges_AMSRE]=meshgrid(lon_edges_sub2,lat_edges_sub);
    [gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE]=meshgrid(lon_centres_sub2,lat_centres_sub);




    switch amsre_action
        case 'process2'
            years = dir([filedir 'y20*']);
        case {'read','read_daily','process'}
            nmonths=0;
            for iy=1:length(years_select)
                years(iy).name = [years_select{iy}];
                nmonths = nmonths + length(months_select{iy});
            end

            switch amsre_action
                case {'read'}
                    %set up the blank NaN arrays for each variable
                    for ivar=1:length(var_names)
                        eval([var_names{ivar} '_amsre = NaN * ones([1440/Nsub_av 720/Nsub_av nmonths 2]);']);
                        %                    sst_amsre = NaN * ones([1440 720 nmonths]);
                        if ivar==1
                            daynum_timeseries3_AMSRE=[1:nmonths];
                            gcm_time_UTC_AMSRE=zeros([1 nmonths]);
                        end

                    end
                case {'read_daily'}
                    %set up the blank NaN arrays for each variable
                    for ivar=1:length(var_names)
                        %for sst will just do the average of the day and night values whereas for LWP will want
                        %to keep them separate
                        if strcmp(var_names{ivar},'sst')==1
                            nascdesc=1;
                        else
                            nascdesc=2;
                        end

                        switch case_study
                            case 'Choose days'
                                eval([var_names{ivar} '_amsre = NaN * ones([length(sub_lon) length(sub_lat) length(chosen_files) nascdesc]);']);
                            otherwise
                                eval([var_names{ivar} '_amsre = NaN * ones([length(sub_lon) length(sub_lat) nmonths*31 nascdesc]);']);
                        end
                        %                sst_amsre = NaN * ones([length(sub_lon) length(sub_lat) nmonths*31]);
                    end
                    %                eval(['lwp2_amsre = NaN * ones([length(sub_lon) length(sub_lat) nmonths*31]);']);
            end

        case 'read_saved'
            load([filedir 'sst_DJF06_JJA07']);
            return
            disp('Done load saved AMSRE data');
    end


    %% Start the actual processing

    idat=0;
    idat2=0;
    for iyear=1:length(years)

        year_dir = [filedir years(iyear).name '/'];

        switch amsre_action
            case 'process2'
                months = dir([year_dir 'm*']);
            case {'read','read_daily','process'}
                clear months
                for im=1:length(months_select{iyear})
                    months(im).name = [months_select{iyear}{im}];
                end
        end

        for imonth=1:length(months)
            idat=idat+1;
            imonth
            mon_dir = [year_dir months(imonth).name];
            switch amsre_action
                case 'read'
                    filepath = [mon_dir '/monthly_mean.mat'];
                    if exist(filepath)==2
                        load(filepath);
                    end

                    month_amsre(idat) = str2num(months(imonth).name(2:end));
                    year_amsre(idat) = str2num(years(iyear).name(2:end));

                    for ivar=1:length(var_names)
                        dat = eval(['reduce_matrix_subsample_mean(' var_names{ivar} '_monav(:,:,1),Nsub_av,Nsub_av);']);
                        clear dat2
                        dat2(1:180,:)=dat(181:360,:);
                        dat2(181:360,:)=dat(1:180,:);
                        eval([var_names{ivar} '_amsre(:,:,idat,1)=dat2;']);

                        dat = eval(['reduce_matrix_subsample_mean(' var_names{ivar} '_monav(:,:,2),Nsub_av,Nsub_av);']);
                        clear dat2
                        dat2(1:180,:)=dat(181:360,:);
                        dat2(181:360,:)=dat(1:180,:);
                        eval([var_names{ivar} '_amsre(:,:,idat,2)=dat2;']);
                        %                    sst_amsre(:,:,idat) = mean(sst_monav,3);
                    end

                case {'process','read_daily'}
                    switch case_study
                        case 'Choose days'
                            %                        files(1).name = [filedir year_single '/' month_single '/' single_file];
                            %                        files(1).name = [single_file];
                            files = chosen_files;
                        otherwise
                            files = dir([filedir years(iyear).name '/' months(imonth).name '/amsre_*v7*']);
                            %                        files = files(3:end); %Remove the . and .. entries
                    end

                    switch amsre_action
                        case 'process'
                            %make arrays for the monthly averages
                            sst_monav   = zeros([1440 720 2]); Nsst_monav   = zeros([1440 720 2]);
                            uLF_monav   = zeros([1440 720 2]); NuLF_monav   = zeros([1440 720 2]);
                            uMF_monav   = zeros([1440 720 2]); NuMF_monav   = zeros([1440 720 2]);
                            vapor_monav = zeros([1440 720 2]); Nvapor_monav = zeros([1440 720 2]);
                            lwp_monav   = zeros([1440 720 2]); Nlwp_monav   = zeros([1440 720 2]);
                            rain_monav  = zeros([1440 720 2]); Nrain_monav  = zeros([1440 720 2]);
                    end

                    nfile=0;

                    for ifile=1:length(files) %loop through all of the files in the dir for this month
                        file_dir2 = [filedir years(iyear).name '/' months(imonth).name '/'];
                        file_name = [file_dir2 files(ifile).name];
                        file_name2 = files(ifile).name;
                        if length(strfind(files(ifile).name,'d3d'))==0 %if not a 3-day file
                            if strfind(files(ifile).name,'.gz')
                                %unzip the file
                                eval(['!gzip -d ' file_name]);
                                file_name = remove_character(file_name,'.gz','');
                                file_name2 = remove_character(file_name2,'.gz','');
                            end

                            %Normal files look like amsre_20081124v7, whereas
                            %the monthly mean files look like:-
                            %amsre_200811v7
                            iu=findstr(file_name2,'_');
                            iv=findstr(file_name2,'v'); %for v7. v7.0.1, etc.


                            if iv(end)-iu(1)~=7 | strcmp(sat,'MAC3')==1 %if not the monthly av file

                                nfile=nfile+1;
                                idat2=idat2+1;

                                %                            day_file = str2num(files(ifile).name(13:14));
                                day_file = str2num(files(ifile).name(iu(1)+7:iu(1)+8));

                                %read the daily file
                                switch sat
                                    case {'AMSRE','TMI'}
                                        [amsre_time,sst,amsre_wspdLF,amsre_wspdMF,amsre_vapor,amsre_lwp,amsre_rain]=read_amsr_day_v7_Dan(file_name);
                                        %   [mingmt,sst,wspdLF,wspdMF,vapor,cloud,rain,wspdAW,wdir]
                                        %   mingmt is gmt time in hours
                                        %   sst  is surface water temperature at depth of about 1 mm in deg C
                                        %   wspdLF is 10 meter surface wind in m/s made using 10.7 GHz channel and above, low frequency channels
                                        %   wspdMF is 10 meter surface wind in m/s made using 18.7 GHz channel and above, medium frequency channels
                                        %   vapor is atmospheric water vapor in millimeters
                                        %   cloud is liquid cloud water in millimeters
                                        %   rain  is rain rate in millimeters/hour
                                        %   wspdAW is 10 meter surface wind for all weather conditions made using 3 algorithms
                                        %   wdir is wind direction oceanographic convention, blowing North = 0 in degrees
                                        %
                                        %  The center of the first cell of the 1440 column and 720 row map is at 0.125 E longitude and -89.875 latitude.
                                        %  The center of the second cell is 0.375 E longitude, -89.875 latitude.
                                        % 		XLAT=0.25*ILAT-90.125
                                        %		XLON=0.25*ILON-0.125

                                        %convert LWP from mm to g/m2. Mass_water/area = depth*rhoW
                                        %amsre_lwp = amsre_lwp * 1e-3 * 1e3 * 1e3; %i.e. convert mm to m, multiply by rhoW and convert kg into g
                                        %                amsre_lwp = amsre_lwp * 1e3; %leave this until later
                                        % so mm are equivalent to kg/m2  -- have left as kg/m2 here

                                    case {'SSMI-f13','SSMI-f15','SSMI-f16','SSMI-f17'}
                                        [amsre_time,wspdMF,vapor,amsre_lwp,rain]=read_ssmi_day_v7_Dan(file_name);
                                        
                                    case {'Windsat'}
                                        [amsre_time,sst,wspdLF,wspdMF,vapor,amsre_lwp,rain,windAW,wdir]=read_windsat_daily_v7_Dan(file_name);
                                        
                                    case {'MAC'}
                                        [amsre_time,wspdMF,vapor,amsre_lwp,rain]=read_MAC_one_month_daily(file_name);
                                    otherwise
                                        fprintf(1,'\n\n*** This type of data has not been added yet!!\nAdd to multi_read_amsre_daily ***\n\n');
                                        return

                                end


                                switch amsre_action
                                    case 'process'

                                        [sst_monav,Nsst_monav] = add_to_running_average(sst,sst_monav,Nsst_monav);
                                        [uLF_monav,NuLF_monav] = add_to_running_average(wspdLF,uLF_monav,NuLF_monav);
                                        [uMF_monav,NuMF_monav] = add_to_running_average(wspdMF,uMF_monav,NuMF_monav);
                                        [vapor_monav,Nvapor_monav] = add_to_running_average(vapor,vapor_monav,Nvapor_monav);
                                        [lwp_monav,Nlwp_monav] = add_to_running_average(amsre_lwp,lwp_monav,Nlwp_monav);
                                        [rain_monav,Nrain_monav] = add_to_running_average(rain,rain_monav,Nrain_monav);

                                    case 'read_daily'
                                        %average day and night - arrays are
                                        %currently [1440 720 2], i.e. lon lat
                                        %ascending/descending

                                        for ivar=1:length(var_names)
                                            var_name = var_names{ivar};
                                            if strcmp(var_name,'sst')==1
                                                %Average over the ascending and
                                                %descending nodes
                                                sst_read = meanNoNan(sst,3);

                                                if ireduce_to_1deg==1
                                                    %reduce resolution 4x by averaging
                                                    tmp = reduce_matrix_subsample_mean(sst_read,4,4);
                                                else
                                                    tmp = sst_read;
                                                end

                                                nlon = size(tmp,1);

                                                %rearrange to put lon=0 in the middle as with MODIS
                                                %N.B. at this stage AMSRE data is
                                                %ordered as (lon,lat) - permute is
                                                %performed later.
                                                clear tmp2
                                                %                                            tmp2(1:180,:)=tmp(181:360,:);
                                                %                                            tmp2(181:360,:)=tmp(1:180,:);
                                                tmp2(1:nlon/2,:)=tmp(nlon/2+1:nlon,:);
                                                tmp2(nlon/2+1:nlon,:)=tmp(1:nlon/2,:);
                                                sst_amsre(:,:,idat2)=tmp2;

                                            else
                                                %for other data keep ascending and descending separate
                                                for iasc=1:2
                                                    eval(['var_read = squeeze(amsre_' var_name '(:,:,iasc));']);

                                                    if ireduce_to_1deg==1
                                                        %reduce resolution 4x by averaging
                                                        tmp_var = reduce_matrix_subsample_mean(var_read,4,4);
                                                    else
                                                        tmp_var = var_read;
                                                    end

                                                    nlon = size(tmp_var,1);

                                                    clear tmp2_var
                                                    %move the longitude halves around
                                                    tmp2_var(1:nlon/2,:)=tmp_var(nlon/2+1:nlon,:);
                                                    tmp2_var(nlon/2+1:nlon,:)=tmp_var(1:nlon/2,:);


                                                    eval([var_name '_amsre(:,:,idat2,iasc)=tmp2_var;']);
                                                end

                                            end

                                        end


                                        month_amsre(idat2) = str2num(months(imonth).name(2:end));
                                        year_amsre(idat2) = str2num(years(iyear).name(2:end));
                                        day_amsre(idat2) = day_file;
                                        %                                    time_amsre(idat2) = time;


                                end

                            end

                        end %if
                    end %ifile - all files in the dir for this month

                    switch amsre_action
                        case 'process'
                            if nfile>0
                                sst_monav = sst_monav ./ Nsst_monav; sst_monav(Nsst_monav==0)=NaN;
                                uLF_monav = uLF_monav ./ NuLF_monav; uLF_monav(NuLF_monav==0)=NaN;
                                uMF_monav = uMF_monav ./ NuMF_monav; uMF_monav(NuMF_monav==0)=NaN;
                                vapor_monav = vapor_monav ./ Nvapor_monav; vapor_monav(Nvapor_monav==0)=NaN;
                                lwp_monav = lwp_monav ./ Nlwp_monav; lwp_monav(Nlwp_monav==0)=NaN;
                                rain_monav = rain_monav ./ Nrain_monav; rain_monav(Nrain_monav==0)=NaN;
                                %now output the data
                                save([file_dir2 'monthly_mean.mat'],'-V7.3','sst_monav','Nsst_monav','uLF_monav','NuLF_monav','uMF_monav','NuMF_monav','vapor_monav','Nvapor_monav','lwp_monav','Nlwp_monav','rain_monav','Nrain_monav');
                            end

                    end

            end

        end %imonth - loop through all the months

    end %iyear






    switch amsre_action
        case {'read_daily'}
            for ivar=1:length(var_names)
                var_name = var_names{ivar};
                if strcmp(var_name,'sst')==1
                    sst_amsre = permute(sst_amsre,[2 1 3]);
                    %remove the extra indices not required
                    sst_amsre(:,:,idat2+1:end)=[];
                    %flip the lat dimension to be consistent with MODIS
                    sst_amsre = flipdim(sst_amsre,1);
                    switch case_study
                        case 'Choose days'
                        otherwise
                            %smooth over 5 days to help remove NaNs from swath gaps
                            sst_amsre_smooth = amsre_block_av_time_function2(sst_amsre,5,year_amsre,month_amsre,day_amsre);
                            %but not for LWP as we want to col-locate with MODIS Terra
                    end
                else
                    eval([var_name '_amsre = permute(' var_name '_amsre,[2 1 3 4]);']);
                    %remove the extra indices not required
                    eval([var_name '_amsre(:,:,idat2+1:end,:)=[];']);
                    %flip the lat dimension to be consistent with MODIS
                    eval([var_name '_amsre = flipdim(' var_name '_amsre,1);']);
                end
            end

            switch case_study
                case 'Choose days'
                otherwise
                    if exist('modisyear_timeseries3')
                        %make an array the same size as the MODIS array
                        amsre_block_average_time
                    end
            end

    end




    clear ioverride_read_amsre
catch amsre_error

    clear ioverride_read_amsre
    rethrow(amsre_error);

end