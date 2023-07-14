function [time_driver,time_out,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM,it] ...
            = read_time_UM(nc_out,tstr,time_in,time_tol)
        
    version_mat=version;
    
    vfind = strfind(version_mat,'R2007b');
    
%if length(vfind)>0
%    case '7.5.0.338 (R2007b)'                
        time=nc_out{tstr}(:);
        t0_str=nc_out{tstr}.time_origin{1};
        cal_str = nc_out{tstr}.calendar{1};
%else
%         time = double(ncread(nc_out,tstr));
%         t0_str = ncreadatt(nc_out,['/' tstr],'time_origin');
%         try
%             cal_str = ncreadatt(nc_out,['/' tstr],'calendar');
%         catch
%             cal_str='';
%         end
% end     
        
        
    
    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];    
    
    
    
    if length(cal_str)>=7 & strcmp( cal_str(1:7),'360_day' )==1
        %360 day caldendar.
        %time is no. days since origin but with only 30 days per month
        
        %Origin time
        [Y,M,D] = datevec(datenum(t0_str2));
        %Just work in no. months knowing that each month is 30 days.
        months = time/30;
        months_whole = floor(months);
        months_rem = months - months_whole;
        
        M2 = M + months_whole; %N.B. can just specify say month 13, which will give
         %Jan of the next year.
         
        % Add the days onto the origin based on 30 days per month 
        days = D + months_rem*30;
        
        %But have to deal with when go over 30 days
        months_whole2 = floor(days/30);
        D2 = days - months_whole2*30;
        M2 = M2 + months_whole2;
        
        time_driver = datenum(Y,M2,D2);
        
        
        
        
        
    else
        time_driver = datenum(t0_str2) + time;
    end
    
    time_out = time_driver;

    it = 1:length(time_driver);
    
%    time_tol = 1/3600/24; %set as 1 minute  *** passed into function now
    if length(time_in)>0 %This is there in case want all times (time-in set to [])
        it=[];
        for i=1:length(time_in)

                ii = find(abs(time_driver-time_in(i)) < time_tol );
                if length(ii)>1
                    error('***More than one matched time ***');
                elseif length(ii)>0 %to avoid when no time is matched
                    it(i) = ii;  
                else
                    error('***No matched time ***');
                end

            time_out(i) = time_driver(it(i));
        end
    end     
    
    daynum_timeseries3_UM = [1:length(time_out)]';
    modisyear_timeseries3_UM = [1:length(time_out)]';
    gcm_time_UTC_UM = [1:length(time_out)]';
    gcm_time_matlab_UM = time_out;
