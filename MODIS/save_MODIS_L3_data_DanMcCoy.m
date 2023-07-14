% runs make_mockL3_variables_DanMcCoy_saved first to define what to save
% Saves the data in savevarname

% -------------------------------------------------------------
%    --- Some previously saved files ---
% -------------------------------------------------------------
%saved 40-60S for 2007 in this file :-
file_load=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/'...
    'CF_0.8_minCTT_173_ice_allowed_SZA_65/terra/mockL3_saved_data_20140424T015323.mat'];
%There is also this file that previously sent him - global data for days Dec 2006 to Nov 2007 :-
file_load2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/'...
    'CF_0.8_minCTT_173_ice_allowed_SZA_65/mockL3_saved_data_20121017T115935.mat'];


% Data with SZA<90 and CF>80
file_load_sza90 = ['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/'...
    'CF_0.8_minCTT_173_ice_allowed_SZA_0-90/terra/mockL3_saved_data_20140506T084406.mat'];

file_load_CF0_sza65 =['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/'...
    'CF_0.0_minCTT_173_ice_allowed_SZA_65/terra/mockL3_saved_data_20140507T091712.mat'];
% ----------------------------------------------------------------------------

fprintf(1,'\n saving Nd_timeseries etc.');

% --- Running script here to define which variables to save :-
make_mockL3_variables_DanMcCoy_saved
% -----------------------------------

pre_str = 'mockL3_saved_data_2013_no_confidence_screening';

%savedir_var specified in MODIS_multi_DAY_processL3L2
    savevarname = [savedir_var pre_str '_' datestr(now,30) '.mat'];    




app_str=')';
% 
% %first save the non-standard ones that were made
% MODIS_varname='N_time3';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end
% 
% MODIS_varname='W_timeseries';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end
% 
% MODIS_varname='H_timeseries';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end
% 
% MODIS_varname='LWC_timeseries';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end
% 
% MODIS_varname='daynum_timeseries3';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end
% 
% MODIS_varname='modisyear_timeseries3';
% if exist(MODIS_varname)
%     eval(['save(''' savevarname ''',''' MODIS_varname '''' app_str]);
%     app_str=',''-APPEND'')';
% end

%now save all of the standard timeseries3 variables as defined in modis_var
for iread_modis=1:length(modis_var)
    
    fprintf(1,'\n saving %d out of %d',iread_modis,length(modis_var));

    MODIS_varname = modis_var{iread_modis};

    eval(['save(''' savevarname ''',''' MODIS_varname ''',''-v7.3''' app_str]);
    
    app_str=',''-APPEND'')';

end


save(savevarname,'MLAT','-APPEND');
save(savevarname,'MLON','-APPEND');

