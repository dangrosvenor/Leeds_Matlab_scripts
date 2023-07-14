fprintf(1,'\n saving Nd_timeseries etc.');


if ~exist('savevarname') %this is set in averages_of_MODIS_files - so don't worry about files being overwritten!
    savedir_var='/home/disk/eos8/d.grosvenor/saved_data_L3/'
%    savevarname = [savedir_var 'timeseries3_TProfiles_CTT_CTP_data_' dataset_modis '_' datestr(now,30)];
    savevarname = [savedir_var '_' dataset_modis '_' datestr(now,30)];    
end



app_str=')';

%first save the non-standard ones that were made
MODIS_varname='Nd_timeseries';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

MODIS_varname='W_timeseries';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

MODIS_varname='H_timeseries';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

MODIS_varname='LWC_timeseries';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

MODIS_varname='daynum_timeseries3';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

MODIS_varname='modisyear_timeseries3';
if exist(MODIS_varname)
    eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
    app_str=',''-APPEND'')';
end

%now save all of the standard timeseries3 variables as defined in modis_var
for iread_modis=1:length(modis_var)
    
    fprintf(1,'\n saving %d out of %d',iread_modis,length(modis_var));

    MODIS_varname = modis_var{iread_modis};

    eval(['save(''' savevarname ''',''' MODIS_varname ''',''-v7.3''' app_str]);
    
    app_str=',''-APPEND'')';

end


save(savevarname,'MLAT','-APPEND');
save(savevarname,'MLON','-APPEND');

