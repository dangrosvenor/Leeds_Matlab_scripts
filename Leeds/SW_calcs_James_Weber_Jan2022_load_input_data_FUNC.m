%Called from SW_calcs_James_Weber_Jan2022_load_input_data_FUNC02

%stash = '01217'; var='UM_m01s01i217_vn1105';
%UM_run = 'u-ck598'; file_suffix='_8_years'; %control run

switch UM_run
    case {'u-co022','u-cr012','u-co023','u-cr013','u-co112','u-cr140','u-co113','u-cr141'}
        %stash = stash(2:end); %remove the first '0'   
        stash_dir ='13_years';
    otherwise
        stash_dir = stash;
end

if length(strfind('0',stash(1)))==1
    stash = stash(2:end);
end
    


load_file = [data_dir UM_run '/' stash_dir '/' stash file_suffix '.nc3'];
nc = netcdf(load_file);
dat = nc{var}(:);
output.gcm_Plat1D_UM = nc{'latitude'}(:); %Make gcm_Plat2D_UM?
output.gcm_Plon1D_UM = nc{'longitude'}(:); %Make gcm_Plat2D_UM?
%[gcm_Plat2D_UM
nc = close(nc);

[output.gcm_Plat2D_UM,output.gcm_Plon2D_UM,output.gcm_Plat2D_edges_UM,output.gcm_Plon2D_edges_UM,output.gcm_area_UM] =...
    UM_lat_lon_intialise_FUNC(output.gcm_Plat1D_UM,output.gcm_Plon1D_UM);