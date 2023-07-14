%match seaice arrays in time
% Creates seaice_time3 and seaice_max_time3 in order to match the times in
% daynum_timeseries3 and modisyear_timeseries3
% Is created to be the same size as Cloud_Fraction_Liquid.timeseries3

%these days are missing from the database, so add them as a known
%exception
exceptions_yrs = []; %[2006 2007];
exceptions_days = []; %[322 332];

clear daynum_seaice
for it=1:length(seaice_array_1deg_datenum)
    daynum_seaice(it) = day_of_year_from_date_func(seaice_array_1deg_datenum(it));
end
[year_seaice,M,D] = datevec(seaice_array_1deg_datenum);

seaice_time3 = NaN*ones(size(Cloud_Fraction_Liquid.timeseries3));
seaice_max_time3 = NaN*ones(size(Cloud_Fraction_Liquid.timeseries3));

for it=1:length(modisyear_timeseries3)
    i = find(year_seaice==modisyear_timeseries3(it) & daynum_seaice==daynum_timeseries3(it));
    if length(i)>0
        inds_save(it) = i;
    end
    
    ain='stop';
    
    if length(i)>0
       seaice_time3(:,:,it) = seaice_array_1deg(:,:,i);
       seaice_max_time3(:,:,it) = seaice_array_1deg_max(:,:,i);
    else %watch out in case can't find the data for a particular day
        %but allow these known exceptions
        iex = find( exceptions_yrs==modisyear_timeseries3(it) & exceptions_days==daynum_timeseries3(it));
        if strcmp(ain,'stop')==1 & length(iex)==0
            fprintf(1,'\n*** WARNING - cannot find seaice data for day %d.\n ***',daynum_timeseries3(it));
            ain=input('Press enter to ignore all future warnings of this');
        end
        
    end
    
end