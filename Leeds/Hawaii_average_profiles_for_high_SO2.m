function [potential_temp_3d_mean_col_PI,potential_temp_3d_mean_col_PD] = Hawaii_average_profiles_for_high_SO2(UM_base_dir,um_case_PI,um_case_PD,varname_str,SO2_col_PD_ALL)

%SO2 thresholds
thresh_SO2 = -1;
thresh_SO2 = 1e-5;
%thresh_SO2 = 1e-4;

%var_UM = 'potential_temp_3d';
var_UM = varname_str;

dir_UM = [UM_base_dir um_case_PI];
dir_UM2 = [dir_UM '/' var_UM];
filename = [dir_UM2 '/merged.nc'];
nc_PI=netcdf(filename);
       
dir_UM = [UM_base_dir um_case_PD];
dir_UM2 = [dir_UM '/' var_UM];
filename = [dir_UM2 '/merged.nc'];
nc_PD=netcdf(filename);

all_cols_PI=0;
all_cols_PD=0;
N_PI = 0;
N_PD = 0;
for it=1:size(SO2_col_PD_ALL,3)
    fprintf(1,'%d ',it);

    potential_temp_3d_it_PI = nc_PI{var_UM}(it,:,:,:);
    potential_temp_3d_it_PD = nc_PD{var_UM}(it,:,:,:);
    inds = find(SO2_col_PD_ALL(:,:,it) > thresh_SO2);
    [arr_out_PI] = indices_apply_2D_to_3D_all_column(inds,potential_temp_3d_it_PI);
    [arr_out_PD] = indices_apply_2D_to_3D_all_column(inds,potential_temp_3d_it_PD);
    
    %If might have NaN data then probably better to save all the columns
    %and avearage after.
    %all_cols_PI = cat(2,all_cols_PI,arr_out_PI);
    %all_cols_PD = cat(2,all_cols_PD,arr_out_PD);
    
    %But if not then can just sum to save memory
    all_cols_PI = all_cols_PI + sum(arr_out_PI,2);
    all_cols_PD = all_cols_PD + sum(arr_out_PD,2);
    N_PI = N_PI + size(arr_out_PI,2);
    N_PD = N_PD + size(arr_out_PD,2);
    
end

nc_PI = close(nc_PI);
nc_PD = close(nc_PD);

potential_temp_3d_mean_col_PI = all_cols_PI/N_PI;
potential_temp_3d_mean_col_PD = all_cols_PD/N_PD;

% 
% %% Get height info - prob a better way than using the slices...
% var_UM = 'SO2_perkg_lon_height_slice_at_ilat=264_lat=19.50';
% 
% um_case=um_case_PI; pole_lat=71.0; pole_lon=17.0; run_type = run_type_DRIVER; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
% opts.lat_var='height'; opts.lon_var='Longitude';
% var_name_out = 'SO2_slice_PI_ALL';
% UM_load_var_commands
% z2 = dat_global.gcm_Plat2D_UM(:,1);
% 
% %% Plot profiles
% figure('color','w')
% hold on
% plot(potential_temp_3d_mean_col_PI,z2/1e3,'b');
% plot(potential_temp_3d_mean_col_PD,z2/1e3,'r');
% legend({'Volc OFF','Volc ON'});
% set(gca,'ylim',[0 6]);






