LAT_val = [72 75]; LON_val = [-3:48]; %Arctic

ilat = find(LAT>=LAT_val(1) & LAT<LAT_val(end));
ilon = find(LON>=LON_val(1) & LON<LON_val(end));


%find the individual components of the date
dat = Date_Time_Swath.timeseries3(ilat,ilon,:);
sza_dat = Solar_Zenith_Mean.timeseries3(ilat,ilon,:);
dat=dat(:);
[Y,MO,D,H,MI,S] = datevec(dat);
%Just use the day (not hours and mins)
all_days = datenum(Y,MO,D);
%            all_days = reshape(all_days,size(N_time3));
           
%find the unique days
[B,I,J]=unique(all_days);

%put the required variables in modis_var
save_or_load='save';
%make_mockL3_variables



isw_max=-9e99;
%loop through all the days we have and find the max number of swaths from
%all the days for array creation
for ib=1:length(B)
    isw=find(all_days==B(ib));
    isw_max=max(isw_max,length(isw));
end

Solar_Zenith_Daily_pdfs.timeseries3 = NaN * ones([isw_max length(B)]);

%now loop again and store the values for each day
            %loop through all of the days

                    for ib=1:length(B)
                        isw=find(all_days==B(ib));
                        %the daily max SZA
                        Solar_Zenith_Daily_pdfs.timeseries3(1:length(isw),ib) = sza_dat(isw);  
                        
                        G{ib} = datestr(B(ib));
                        
                    end
                    

