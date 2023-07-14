  %find the individual components of the date
            [Y,MO,D,H,MI,S] = datevec(Date_Time_Swath.timeseries3(:));
            %Just use the day (not hours and mins)
            all_days = datenum(Y,MO,D);   
            all_days = reshape(all_days,size(N_time3));
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
           %remove the days for the points that are being screened
%            all_days(ihtot) = NaN;
            all_days(isnan(dat_modis))=NaN; %to also remove the points that were originally NaN
            %constrain by required times
            all_days = all_days(:,:,time_inds_average);
                        
            a=zeros(size(all_days));
            max_sza=zeros(size(all_days(:,:,1)));

            %loop through all of the days
            for ii=1:size(all_days,1)
                for jj=1:size(all_days,2)
                    dat_ij = all_days(ii,jj,:);
                    inot_nan = find(isnan(dat_ij)==0);

                    %for each lat/lon find the frequencies of each unique
                    %day (J below)
                    %[B,I,J] = UNIQUE(...) also returns index vectors I and J such
                    %that B = A(I) and A = B(J) (or B = A(I,:) and A = B(J,:)).
                    [B,I,J]=unique(dat_ij(inot_nan)); %unique doesn't work for NaNs - all are still included
                    %
                    
                    SZA=NaN; %then min(SZA) will result in NaN below (i.e. no data)
                    %loop through all individual days
                    for kk=1:length(B)
                        id=find(dat_ij==B(kk));
                        %the daily max SZA
                        SZA(kk) = max(Solar_Zenith_Mean.timeseries3(ii,jj,time_inds_average(id)));
                    end
                        
%                    a(ii,jj,time_inds_average(inot_nan(I)))=1;
%                    [mode_day,nmode_freq] = mode(J); %this gets the most frequent day and the number of occurences
                    %so we can know the max no. swaths per day
                    
                    %the minimum for all days of the daily max SZA
                    max_sza(ii,jj) = min(SZA);
                end
            end