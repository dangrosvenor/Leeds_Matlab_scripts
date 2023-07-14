%declare the arrays to start with
for ivar=1:length(vars_PAR)
    eval_str = [vars_PAR{ivar} '_MODIS = NaN*ones(size(Date_Time_Swath.timeseries3));']; eval(eval_str);

end

                                
%loop through all of the days
for it=1:length(modisyear_timeseries3_MODIS)
    
    i = find(modisyear_timeseries3_POLDER==modisyear_timeseries3_block(it) & daynum_timeseries3_POLDER==daynum_timeseries3_MODIS(it));
    if length(i)>0
       sst_amsre_time3(:,:,it) = sst_amsre_smooth(ilat,ilon,i);
       %lwp_amsre is larger in the 3rd dimension than sst_asmre_smooth, but
       %the end of it will just be padded with zeros as is created as
       %nmonths*31, whereas smooth one is only the size of the actual data
       %since it is based on year_amsre, which grows each loop
       if exist('lwp_amsre')
           lwp_amsre_time3(:,:,it,:) = lwp_amsre(ilat,ilon,i,:);
       end
    end
    
end
