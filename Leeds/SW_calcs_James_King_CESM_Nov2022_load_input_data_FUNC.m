%Called from SW_calcs_James_King_CESM_Nov2022_load_input_data_FUNC02

%file = dir([data_dir UM_run '/' file_prefix '*' var '_yrs*.nc']);
file = dir([data_dir UM_run '/' file_prefix '*' var '.nc']);

if length(file)==0
    file = dir([data_dir UM_run '/' file_prefix '*' var '_*.nc']);
    
    
    for i=1:length(file)
        if length(strfind(file(i).name,'_d1'))==0
            filename = file(i).name;
        end
    end
    
else
    
    filename = file(1).name;
    
end


%load_file = [data_dir UM_run '/' file_prefix stash '.nc'];
nc = netcdf([data_dir UM_run '/' filename]);
dat = nc{var}(:);
output.gcm_Plat1D_UM = nc{'lat'}(:); %Make gcm_Plat2D_UM?
output.gcm_Plon1D_UM = nc{'lon'}(:); %Make gcm_Plat2D_UM?
output.plevs = nc{'lev'}(:);
output.plevs_bnds = nc{'lev_bnds'}(:);
output.pilevs = nc{'ilev'}(:);
output.pilevs_bnds = nc{'ilev_bnds'}(:);
%[gcm_Plat2D_UM
nc = close(nc);

[output.gcm_Plat2D_UM,output.gcm_Plon2D_UM,output.gcm_Plat2D_edges_UM,output.gcm_Plon2D_edges_UM,output.gcm_area_UM] =...
    UM_lat_lon_intialise_FUNC(output.gcm_Plat1D_UM,output.gcm_Plon1D_UM);