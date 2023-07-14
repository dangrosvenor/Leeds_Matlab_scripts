%		? Data loaded using
% 			§ MODIS_multi_DAY_processL3L2
% 				? For variables used (multiL2L3_case) used :-
% 					® make_mockL3_variables_mockL3_basic 
% 				? For screening case used (daily_averaged_files_loc2) used :-
% 					® '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/no_screening/'


savedir01 = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/';                
savefile=[savedir01 'max_SZA_vs_lat_from_mockL3_TERRA.mat'];


set_month = {1,2,3,4,5,6,7,8,9,10,11,12};




clear max_SZA_vs_lat
for im=1:length(set_month)
                month_modis = set_month{im};                                
                              
                %Pick out the correct day nunbers of the month in
                %question
                [days_required_for_mean,time_mean_str] = days_of_month(month_modis);

                clear idays
                icount=0;
                for i=1:length(days_required_for_mean)
                    idays_temp = find(daynum_timeseries3_MODIS==days_required_for_mean(i));                
                    if length(idays_temp)>0
                        icount=icount+1;
                        idays(icount)=idays_temp;
                    end
                end
                
                max_SZA_vs_lat_temp = max(Solar_Zenith_Maximum.timeseries3(:,:,idays),[],3); %Max over the month
                    %N.B. the Maximum2 variable has fewer gaps than Maximum
                    %(can't remember the reason now!)
                    %But only seems to be available for the first day??.
                max_SZA_vs_lat(:,im) = max(max_SZA_vs_lat_temp(:,:),[],2); %Max over lon    
                
end


save(savefile,'max_SZA_vs_lat','MLAT','-V7.3');

                