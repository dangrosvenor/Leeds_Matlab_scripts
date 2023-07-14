function [dat_out] = time_match_data(time01,time02,dat02,exception_yrs,exception_days)
%[dat_out] = time_match_data(time01,time02,dat02,exception_yrs,exception_days)
%match the data in dat02 to the times for the data in dat01 (times are time01 and time02)
%data arrays are in 3D (lat,lon,time).
%Time arrays are either corresponding 3D arrays in "Matlab" time.
%OR they are daynum_timeseries3 and modisyear_timeseries3 (i.e. vectors of
%the day-of-the-year and year). If inputting like this then enter as a cell
%array with the day first. E.g.
%{daynum_timeseries3_MODIS,modisyear_timeseries3}

if ~exist('exception_yrs')
    %these days are missing from the database, so add them as a known
    %exception
    exceptions_yrs = []; %[2006 2007];
    exceptions_days = []; %[322 332];
end

%Doing this since each time entry will be for the specific overpass time,
%but for a given day they should globally be the same day (based on UTC
%time). So am just creating a 1D vector here.
if iscell(time01)
    daynum01 = time01{1}(:);
    year01 = time01{2}(:);
else

    %time_vec01 = meanNoNan(meanNoNan(time01,1),1);
    
    %get the day of the year from the Matlab time.
    for it=1:length(time01)
        %daynum01(it,1) = day_of_year_from_date_func(time_vec01(it));
        daynum01(1,it) = day_of_year_from_date_func(time01(it));        
    end


    %Get the year
    %[year01,M,D] = datevec(time_vec01);
    [year01,M,D] = datevec(time01);    


end

if iscell(time02)
    daynum02 = time02{1}(:);
    year02 = time02{2}(:);
else

    %time_vec02 = meanNoNan(meanNoNan(time02,1),1);

    %get the day of the year from the Matlab time.
    for it=1:length(time02)
        %daynum02(it,1) = day_of_year_from_date_func(time_vec02(it));
        daynum02(1,it) = day_of_year_from_date_func(time02(it));        
    end

    %Get the year
    %[year02,M,D] = datevec(time_vec02);
    [year02,M,D] = datevec(time02);
end


%data will be time matched to time01, so time dimension will be the same as
%dat01
dat_out = NaN*ones([size(dat02,1) size(dat02,2) length(year01)]);

ain='stop';
for it=1:length(year01)
    i = find(year02==year01(it) & daynum02==daynum01(it));
    if length(i)>0
        dat_out(:,:,it) = dat02(:,:,i);
    else %watch out in case can't find the data for a particular day
        %but allow these known exceptions
        iex = find( exceptions_yrs==year01(it) & exceptions_days==daynum01(it));
        if strcmp(ain,'stop')==1 & length(iex)==0    
            fprintf(1,'\n*** WARNING - cannot find seaice data for day %d.\n ***',daynum01(it));
%            ain=input('Press enter to ignore all future warnings of this');
             fprintf(1,'\n*** continuing....\n ***');
             ain='';
        end

    end

end