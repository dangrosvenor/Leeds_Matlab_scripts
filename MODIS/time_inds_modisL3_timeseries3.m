try
    
    if ~exist('ioverride_years_time_screen') | ioverride_years_time_screen==0        
        years_required_for_mean=[1990:2020]; %cover lots of years to make sure that all are used.
    end
    

    if ~exist('ioverride_time_selection') | ioverride_time_selection==0

        days_required_for_mean = [1:366]; time_mean_str = 'ANNUAL';

        time_mean_str = 'ALL'; %all days - same as for annual except different string for the title
        
%        time_mean_str = 'choose_years'; years_required_for_mean = [2007];

        %for seasonal means might want to do "solar-centric" months - i.e. centred
        %around 21st of June and 21st Dec - so for winter that would actually
        %be 5th Nov to 5th Feb, so more like NDJ than DJF. Although I guess
        %that the tempeature response is not really changing like that? Maybe
        %should do 21st Dec - 21st Feb? Or 21st Nov - 21st Jan, but call it
        %DJF? Or just do DJF??

        %straight DJF etc.
%                 days_required_for_mean = [336:366 1:60]; time_mean_str = 'DJF';
        %    %straight MAM
        %        days_required_for_mean = [61:152]; time_mean_str = 'MAM';
        %    %straight JJA
%                days_required_for_mean = [153:244]; time_mean_str = 'JJA';
        %    %straight SON
%                days_required_for_mean = [245:335]; time_mean_str = 'SON';

        %solar DJF etc.
        %         days_required_for_mean = [354-45:366 1:mod(354+45,366)]; time_mean_str = 'solar DJF';
        %    %straight MAM
        %        days_required_for_mean = [79-45:79+45]; time_mean_str = 'solar MAM';
        %    %straight JJA
        %        days_required_for_mean = [171-46:171+45]; time_mean_str = 'solar JJA';
        %    %straight SON
        %        days_required_for_mean = [263-46:263+45]; time_mean_str = 'solar SON';


        %
%             days_required_for_mean = [1:31]; time_mean_str = 'Jan';
        %     days_required_for_mean = [32:60]; time_mean_str = 'Feb';
%              days_required_for_mean = [61:91]; time_mean_str = 'Mar';
%         days_required_for_mean = [92:121]; time_mean_str = 'Apr';
%        days_required_for_mean = [122:152]; time_mean_str = 'May';
        %     days_required_for_mean = [153:182]; time_mean_str = 'Jun';
        % days_required_for_mean = [183:213]; time_mean_str = 'Jul';
%             days_required_for_mean = [214:244]; time_mean_str = 'Aug';
%        days_required_for_mean = [245:274]; time_mean_str = 'Sep';
        %     days_required_for_mean = [275:305]; time_mean_str = 'Oct';
        % days_required_for_mean = [306:335]; time_mean_str = 'Nov';
%            days_required_for_mean = [336:366]; time_mean_str = 'Dec';
        %    days_required_for_mean = [289:320]; time_mean_str = 'VOCALS';

        %    days_required_for_mean = [2]; time_mean_str = ['t='
        %    num2str(days_required_for_mean)];
        %    days_required_for_mean = [183]; time_mean_str = ['t='
        %    num2str(days_required_for_mean)];
        %    days_required_for_mean = [92]; time_mean_str = ['t=' num2str(days_required_for_mean)];
        %    days_required_for_mean = [12:12:48 1:12:48 2:12:48]; time_mean_str = ['t=' num2str(days_required_for_mean)];

        %    days_required_for_mean = [354-15:354+15]; time_mean_str = 'Dec'; %time_mean_str = ['1 month centred on 20th Dec'];  %
        %    days_required_for_mean = [79-15:79+15]; time_mean_str = 'Mar'; %time_mean_str = ['1 month centred on 20th Jun'];  %
        %    days_required_for_mean = [171-15:171+15]; time_mean_str = 'Jun'; %time_mean_str = ['1 month centred on 20th Jun'];  %
        %    days_required_for_mean = [263-15:263+15]; time_mean_str = 'Sep'; %time_mean_str = ['1 month centred on 20th Jun'];  %

        %    days_required_for_mean = [79]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %20th March
        %    days_required_for_mean = [171]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %20th June
        %    days_required_for_mean = [263]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %20th Sep
        %    days_required_for_mean = [354]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %20th Dec
        %    days_required_for_mean = [336]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %1st Dec 2008
        %    days_required_for_mean = [60]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %29th Feb 2008
        %    days_required_for_mean = [19]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %19th Jan
        %    days_required_for_mean = [165]; time_mean_str = ['t=' num2str(days_required_for_mean)];  %14th June
%            days_required_for_mean = [210]; time_mean_str = ['t=' num2str(days_required_for_mean)];  % arbitrary single day       
%            days_required_for_mean = [1]; time_mean_str = ['t=' num2str(days_required_for_mean)];
        %    days_required_for_mean = [1 2 12 13 14 24 25 26 36 37 38 48]; time_mean_str = ['DJF'];  %DJF Calipso
        %    days_required_for_mean = [6 7 8  6+12 7+12 8+12  6+24 7+24 8+24  6+36 7+36 8+36]; time_mean_str = ['JJA'];  %DJF Calipso
%            days_required_for_mean = [199]; time_mean_str = ['17th July, 2008'];  %17th July,2008
        %    days_required_for_mean = [1:size(Solar_Zenith_Maximum_Daily.timeseries3,3)]; time_mean_str = ['ALL days (from daily) '];
        %special indices for data that has been consolidated into daily data

        switch time_series_type
            case {'monthly calipso precip','CALIPSO','AMSRE'}
                years_required_for_mean = [2007:2008];
                years_required_for_mean = [2008];                
%                years_required_for_mean = [2006:2010];

                if ~exist('ioverride_time_selection') | ioverride_time_selection==0
                    matt_year_str = [num2str(years_required_for_mean(1)) '-' num2str(years_required_for_mean(end))];
                    %        months_required_for_mean = [1:12]; time_mean_str = ['month=' num2str(months_required_for_mean)];
                    months_required_for_mean = [1:12]; time_mean_str = ['ANNUAL, ' matt_year_str];
                    %            months_required_for_mean = [12 1 2]; time_mean_str = ['DJF, ' matt_year_str];
%                                months_required_for_mean = [3 4 5]; time_mean_str = ['MAM, ' matt_year_str];
%                                months_required_for_mean = [6 7 8]; time_mean_str = ['JJA, ' matt_year_str];
                    %            months_required_for_mean = [9 10 11]; time_mean_str = ['SON, ' matt_year_str];
                    
                    


%                    years_required_for_mean = unique(year_amsre); %if want
%                    all the available AMSRE years
                    
                    

                end

                %         case {'CALIPSO','AMSRE'}
                %            months_required_for_mean = [12 1 2]; time_mean_str = ['DJF'];

        end




    else
        %    clear ioverride_time_selection
    end

if exist('time_mean_str')
    time_mean_str_orig = time_mean_str;
end

    switch time_series_type
        case 'daily'
            
            
            


            switch time_mean_str
                case {'ANNUAL','ALL'}
                    time_inds_average = [1:length(daynum_timeseries3)];
                case 'choose_years'
                    time_inds_average=[];
                    for iyear_mean = 1:length(years_required_for_mean)
                        iyears_find = find(modisyear_timeseries3==years_required_for_mean(iyear_mean));
                        %                    time_inds_average = [time_inds_average idays_find];
                        if length(iyears_find)>0
                            time_inds_average = cat(1,time_inds_average,iyears_find(:));
                        end
                    end

                    if size(time_inds_average,1)~=1
                        time_inds_average=time_inds_average';
                    end
                    
                otherwise

                    time_inds_average=[];
                    
                    if exist('days_required_for_mean') & exist('modisyear_timeseries3')
                        for iyear_mean = 1:length(years_required_for_mean)
                            for iday_mean = 1:length(days_required_for_mean)
                                idays_find = find(daynum_timeseries3==days_required_for_mean(iday_mean) & modisyear_timeseries3==years_required_for_mean(iyear_mean));
                                %                    time_inds_average = [time_inds_average idays_find];
                                if length(idays_find)>0
                                    time_inds_average = cat(1,time_inds_average,idays_find(:));
                                end
                            end
                        end
                    else
                        fprintf(1,'\nWARNING - time_inds_average not set in time_inds_modisL3_timeseries3\n');
                    end

                    if size(time_inds_average,1)~=1
                        time_inds_average=time_inds_average';
                    end

            end
            



            if ~exist('gcm_time_of_day_select')
                gcm_time_of_day_select=0;
            end

            switch gcm_time_of_day_select
                case 0
                    time_inds_average2 = time_inds_average;
                case 1
                    times_required = [0 3 6 9 12 15 18 21];
                    %                times_required = [15 18 21];

                    time_inds_average2=[];
                    time_UTC_str=' ';
                    for itime_mean = 1:length(times_required)
                        itimes_find = find(gcm_time_UTC(time_inds_average)==times_required(itime_mean));
                        time_inds_average2 = [time_inds_average2 time_inds_average(itimes_find)];
                        time_UTC_str=[time_UTC_str num2str(times_required(itime_mean))  ','];
                    end

                    time_inds_average = time_inds_average2;


                        time_mean_str = [time_mean_str time_UTC_str(1:end-1) ' UTC'];


                case 2 %longitude dependent local time selection

                    if ~exist('ioverride_time_selection') | ioverride_time_selection==0

                        %remember that we only have 0,3,6 etc. UTC hours from the
                        %model output - although this will shifted depending on
                        %longitude. But if we just choose one UTC time then

                        %create a 3D array [time lat lon] of indices to make NaN
                        %(i.e. those not at the required times)
                        times_required = [0:24]; %need to specify all hours (in 1 hour increments)
                        %24 not needed, but 0-24 sounds better for time_UTC_str
                        %                times_required = [15:21];
                        %                times_required = [11:16];


                        %               times_required = [10:11]; %stricter Terra daytime (more suitable for AMSRE with its narrower swath?)
                        %but need to choose a range of at least 3 hours tomake
                        %sure that each time zone has a possible output time
                        %(since only have output every 3 hours).
%                                       times_required = [9:12]; %Terra daytime
                        %                times_required = [21:23 0]; %Terra nighttime

%                        times_required = [12:15];  %Aqua daytime
%                        times_required = [0:3]; %Aqua nighttime

%EMCWF
%                         times_required = [20:23 0:6];
%                         times_required = [7:19];
%                         times_required = [9:15];
%                          times_required = [21:23 0:3];



                        %                 times_required = [9:15];  %Aqua/Terra daytime
                        %                 times_required = [21:23 0:3];  %Aqua/Terra nighttime

                        %                times_required = [0:3 12:15];

                    end

                    %first we need a 3D array of the local times based on the longitude
                    SIZ = size(Plon2D);
                    %                SIZ_time = length(gcm_time_UTC(time_inds_average));
                    SIZ_time = length(gcm_time_UTC);
                    times_UTC_3D = repmat(gcm_time_UTC,[1 SIZ(1) SIZ(2)]);
                    days_3D = repmat(daynum_timeseries3,[1 SIZ(1) SIZ(2)]);
                    Plon_3D = repmat(Plon2D,[1 1 SIZ_time]);
                    Plon_3D = permute(Plon_3D,[3 1 2]);
                    %round to the nearest hour for simplicity (15 gegrees per
                    %hour). The mod prevents hours <0 or >24
                    times_local_3D = mod( round(times_UTC_3D + Plon_3D/15) , 24 );

                    siz3D = size(times_local_3D);
                    LT = length(time_inds_average);

                    %taking the approach of a 3D matrix that contains all NaNs
                    %except at the required times when it is zero. Then this can be added to
                    %dat_modis, which will make the data that we don't want to include NaN
                    time_inds_average2=NaN*ones(size(times_local_3D));
                    time_UTC_str=' ';

                    %                for iday_mean = 1:length(days_required_for_mean)
                    %create a set of lat and lon inds for each time value
                    if size(time_inds_average,1)==1
                        IT=repmat(time_inds_average',[1 siz3D(2) siz3D(3)]);
                    else
                        IT=repmat(time_inds_average,[1 siz3D(2) siz3D(3)]);
                    end


                    [ILAT2,ILON2] = meshgrid([1:siz3D(2)],[1:siz3D(3)]);

                    %ILAT=repmat([1:siz3D(2)], [siz3D(3) 1 LT]);
                    ILAT=repmat(ILAT2, [1 1 LT]);
                    ILAT=permute(ILAT,[3 1 2]);

                    %ILON=repmat([1:siz3D(3)], [siz3D(2) 1 LT]);
                    ILON=repmat(ILON2, [1 1 LT]);
                    ILON=permute(ILON,[3 1 2]);

                    itdays = sub2ind(siz3D,IT(:),ILAT(:),ILON(:));
                    %these are now linear indices in times_local_3D just for the days we want


                    for itime_mean = 1:length(times_required)
                        itimes_find = find(times_local_3D(itdays)==times_required(itime_mean) );
                        %make the days that we want zero (replace the NaNs)
                        time_inds_average2(itdays(itimes_find))=0;
                        %                    time_UTC_str=[time_UTC_str num2str(times_required(itime_mean))  ','];
                    end



                    time_UTC_str=[' ' num2str(times_required(1)) '-' num2str(times_required(end)) ' LOCAL TIME ']

                    %                time_inds_average = time_inds_average2;
                    

                        time_mean_str = [time_mean_str time_UTC_str];


                case 22 %OLD longitude dependent local time selection
                    %create a 3D array [time lat lon] of indices to make NaN
                    %(i.e. those not at the required times)
                    times_required = [0:24]; %need to specify all hours (in 1 hour increments)
                    %24 not needed, but 0-24 sounds better for time_UTC_str
                    %                times_required = [15:21];
                    %                times_required = [11:16];

                    %                times_required = [9:12]; %Terra daytime
                    %                times_required = [21:23 0]; %Terra nighttime

                    %                times_required = [12:15];  %Aqua daytime
                    %                times_required = [0:3]; %Aqua nighttime


                    %                 times_required = [9:15];  %Aqua/Terra daytime
                    %                 times_required = [21:23 0:3];  %Aqua/Terra nighttime

                    %first we need a 3D array of the local times based on the longitude
                    SIZ = size(Plon2D);
                    SIZ_time = length(gcm_time_UTC(time_inds_average));
                    times_UTC_3D = repmat(gcm_time_UTC(time_inds_average),[1 SIZ(1) SIZ(2)]);
                    Plon_3D = repmat(Plon2D,[1 1 SIZ_time]);
                    Plon_3D = permute(Plon_3D,[3 1 2]);
                    %round to the nearest hour for simplicity (15 gegrees per
                    %hour). The mod prevents hours <0 or >24
                    times_local_3D = mod( round(times_UTC_3D + Plon_3D/15) , 24 );

                    %taking the approach of a 3D matrix that contains all NaNs
                    %except at the required times when it is zero. Then this can be added to
                    %dat_modis, which will make the data that we don't want to include NaN
                    time_inds_average2=NaN*ones(size(times_local_3D));
                    time_UTC_str=' ';
                    for itime_mean = 1:length(times_required)
                        itimes_find = find(times_local_3D==times_required(itime_mean));
                        time_inds_average2(itimes_find)=0;
                        %                    time_UTC_str=[time_UTC_str num2str(times_required(itime_mean))  ','];
                    end

                    time_UTC_str=[' ' num2str(times_required(1)) '-' num2str(times_required(end)) ' LOCAL TIME ']

                    %                time_inds_average = time_inds_average2;

                    

                        time_mean_str = [time_mean_str time_UTC_str];



            end

    
        case {'monthly calipso precip','CALIPSO','AMSRE'}
         %for switch time_series_type

            switch time_series_type
                case {'AMSRE'}
                    month_calipso_matt = month_amsre;
                    year_calipso_matt = year_amsre;
                    iasc_desc=0;
%                    years_required_for_mean = unique(year_amsre);
                case {'CALIPSO'}
                    month_calipso_matt = month_calipso_cf;
                    year_calipso_matt = year_calipso_cf;
                    iasc_desc=0;
                    years_required_for_mean = unique(year_calipso_cf);
                case {'ERAInt'}
                    month_calipso_matt = month_ERAInt;
                    year_calipso_matt = year_ERAInt;
                    iasc_desc=0;
                    years_required_for_mean = unique(year_ERAInt)
            end


            switch time_mean_str
                case {'ANNUAL','ALL'}
                    time_inds_average = [1:length(daynum_timeseries3)];
                otherwise
                    time_inds_average=[];
                    for iyear_mean = 1:length(years_required_for_mean)
                        for imonth_mean = 1:length(months_required_for_mean)
                            imonths_find = find(month_calipso_matt==months_required_for_mean(imonth_mean) & year_calipso_matt==years_required_for_mean(iyear_mean));
                            if iasc_desc==1
                                imonths_find=imonths_find(1:2:end);
                                if imonth_mean==1
                                    time_mean_str = [time_mean_str ' ASCENDING ONLY'];
                                end
                            elseif iasc_desc==2
                                imonths_find=imonths_find(2:2:end);
                                if imonth_mean==1
                                    time_mean_str = [time_mean_str ' DESCENDING ONLY'];
                                end
                            end %otherwise take all of them (e.g. CALIPSO, AMSRE)
                            time_inds_average = [time_inds_average imonths_find];
                        end
                    end

            end






    end




    clear ioverride_time_selection ioverride_years_time_screen
catch time_error
    clear ioverride_time_selection ioverride_years_time_screen
    rethrow(time_error);
end