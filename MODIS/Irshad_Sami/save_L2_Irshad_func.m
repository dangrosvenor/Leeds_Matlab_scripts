function save_L2_Irshad_func(dat,savevarname,MODIS_varname,app_str,icell)

% if strcmp(MODIS_varname,'MODIS_swath_filename')
%     '';
% end

eval_str = [MODIS_varname '=squeeze(dat);'];
eval(eval_str);

eval_str = ['save(''' savevarname ''',''' MODIS_varname '''' app_str]; 
eval(eval_str);