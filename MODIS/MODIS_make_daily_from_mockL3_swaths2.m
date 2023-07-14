  %find the individual components of the date
            [Y,MO,D,H,MI,S] = datevec(Date_Time_Swath.timeseries3(:));
            %Just use the day (not hours and mins)
            all_days = datenum(Y,MO,D);   
            all_days = reshape(all_days,size(N_time3));
           
%find the unique days
[B,I,J]=unique(all_days);

%put the required variables in modis_var
save_or_load='save';
make_mockL3_variables

Solar_Zenith_Maximum_Daily.timeseries3 = NaN * ones([size(all_days,1) size(all_days,2) length(B)]);
Solar_Zenith_Minimum_Daily.timeseries3 = NaN * ones([size(all_days,1) size(all_days,2) length(B)]);
Solar_Zenith_CosMean_Daily.timeseries3 = NaN * ones([size(all_days,1) size(all_days,2) length(B)]);


for ivar=1:length(modis_var)
    eval_str = [modis_var{ivar} '_Daily.timeseries3 = NaN * ones([size(all_days,1) size(all_days,2) length(B)]);'];
%Solar_Zenith_Mean_Daily.timeseries3 = NaN * ones([size(all_days,1) size(all_days,2) length(B)]);
    eval(eval_str);
end




            %loop through all of the days
            for ii=1:size(all_days,1)
                ii
                for jj=1:size(all_days,2)                                  
                    %loop through all individual days
                    for kk=1:length(B)
                        id=find(all_days(ii,jj,:)==B(kk));
                        %the daily max SZA
                        Solar_Zenith_Maximum_Daily.timeseries3(ii,jj,kk) = max(Solar_Zenith_Mean.timeseries3(ii,jj,id));
                        Solar_Zenith_Minimum_Daily.timeseries3(ii,jj,kk) = min(Solar_Zenith_Mean.timeseries3(ii,jj,id));
                        Solar_Zenith_MeanCos_Daily.timeseries3(ii,jj,kk) = 180/pi*acos(meanNoNan(cos(pi/180*Solar_Zenith_Mean.timeseries3(ii,jj,id)),3));
                        
                        for ivar=1:length(modis_var)
                            eval_str=[modis_var{ivar} '_Daily.timeseries3(ii,jj,kk) = meanNoNan(' modis_var{ivar} '.timeseries3(ii,jj,id),3);'];
                            %                        Solar_Zenith_Mean_Daily.timeseries3(ii,jj,kk) = meanNoNan(Solar_Zenith_Mean.timeseries3(ii,jj,id),3);
                            eval(eval_str);
                        end
                        
                    end
                                            
                end
            end