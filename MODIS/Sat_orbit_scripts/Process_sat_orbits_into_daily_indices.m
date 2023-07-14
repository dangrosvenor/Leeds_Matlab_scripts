%%  Process into daily indexed arrays (from arrays with an index for each
%%  TLE)
for itle=1:size(sat_times,1)
    for ilat=1:size(sat_times,3) %only one index for lat/lon
        %sat_times are in days starting from tspan(1)
        %        sat_times = sat_times + epoch_days_thisyear(itle) + datenum(years(iyear),1,1)-1;
        %        sat_times = sat_times + epoch_days_thisyear(itle);
        D = floor(squeeze(sat_times(itle,:,ilat)));
        %make into days from the start of the current year
        %N.B. epoch_days are decimal days (i.e they include the time too)
        %        [Y,MO,D,H,MI,S] = datevec(sat_times);
        %days since the start of the current year
        %        sat_times = sat_times - datenum(years(iyear),1,1) - 1;
        D(isnan(D))='';
        day_uni = unique(D);

        for idaysat=1:length(day_uni)
            dsat = day_uni(idaysat);
            idaysati=find(D==dsat);
            sat_times_daily(dsat,1:length(idaysati),ilat) = sat_times(itle,idaysati,ilat);
            sat_sza_daily(dsat,1:length(idaysati),ilat) = sat_sza(itle,idaysati,ilat);
            sat_sensor_daily(dsat,1:length(idaysati),ilat) = sat_sensor(itle,idaysati,ilat);                   
        end

    end
end