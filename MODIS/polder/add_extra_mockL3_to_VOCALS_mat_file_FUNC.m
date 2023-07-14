function add_extra_mockL3_to_VOCALS_mat_file_FUNC(dat,file,varname,varname2,ilat,ilon,it)
%varname2 is e.g. '.timeseries3' or just '' if no extension

load(file,'modisyear_timeseries3_MODIS')


eval_str = ['load(file,''' varname ''');'];
eval(eval_str);

if ~exist(varname)
    if strcmp(varname,'lwp_amsre_time3')==1
        eval_str = [varname varname2 ' = NaN*ones([length(ilat) length(ilon) length(modisyear_timeseries3_MODIS) 2]);'];        
    else
        eval_str = [varname varname2 ' = NaN*ones([length(ilat) length(ilon) length(modisyear_timeseries3_MODIS)]);'];
    end
    
    
    eval(eval_str);
end

eval_str = [varname varname2 '(ilat,ilon,it,:) = squeeze(dat' varname2 ');'];
eval(eval_str);

save(file,varname,'-V7.3','-APPEND');



