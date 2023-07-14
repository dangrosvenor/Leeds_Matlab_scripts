function [F_daily,F_mon,time_daily] = ecmwf_averages(F,time) 

ecDateVec_Man = datevec(time);    %[Y,MO,D,H,MI,S]




%average into daily values
        sMan = size(F);
        colon_str='';
        nan_str = '';
        for i=1:length(sMan)-1
            %colon_str = ',:,:,:';
            colon_str = [colon_str ',:'];
            nan_str = [nan_str ' sMan(' num2str(i+1) ')'];
        end

        istart=1; stride=4;
        n = sMan(1)/stride;
%        F_daily = NaN*ones([n-1 sMan(2) sMan(3) sMan(4)]);
        eval_str = ['F_daily = NaN*ones([n-1 ' nan_str ']);'];
        eval(eval_str);
        
       for i=1:n-1
          % F_daily(i,:,:,:) = meanNoNan(F(istart:istart+stride-1,:,:,:),1);
          eval_str =['F_daily(i' colon_str ') = meanNoNan(F(istart:istart+stride-1' colon_str '),1);'];
          eval(eval_str);
          istart = istart + stride;
       end

        time_daily = time(1:stride:n*stride);
        
        %make monthly averages for each time of day (0, 6 , 12 and 18 UTC).
%        F_mon = NaN*ones([4 12 sMan(2) sMan(3) sMan(4)]);
        eval_str = ['F_mon = NaN*ones([4 12 ' nan_str ']);'];
        eval(eval_str);
        
        hours=[0 6 12 18];
for i=1:12
    imon = find(ecDateVec_Man(:,2)==i);
    for j=1:4
        itime = find(ecDateVec_Man(imon,4)==hours(j));
        %F_mon(j,i,:,:,:) = meanNoNan(F(imon(itime),:,:,:),1);
        eval_str = ['F_mon(j,i' colon_str ') = meanNoNan(F(imon(itime)' colon_str '),1);'];
        eval(eval_str);
        
    end
end


        
        
        
        