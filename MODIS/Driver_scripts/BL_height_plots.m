cont_col_str_DRIVER='k';
iplot_mgrid_lines_DRIVER=1;
irestrict_domain_DRIVER=0;
icoarse_grain=0;
ioverride_ticks_DRIVER=0;

var_UM = 'BL_type';    
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[BL_type_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
%Stash 3-476

var_UM = 'BL_height';    
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PI emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[BL_height_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
% This is stash 3-025; Looks like this only shows the top of the lowermost
% layer for decoupled BL types (types 2,4 and 5). Will need stash 3-304 to
% get the top of the decoupled cloud above.

var_UM = 'BL_mixing_height';    
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[BL_mixing_height_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
% This is stash 3-304; Looks like this only shows the top of the lowermost
% layer for decoupled BL types (types 2,4 and 5). Will need stash 3-304 to
% get the top of the decoupled cloud above.

var_UM = 'BL_parcel_height';    
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[BL_parcel_height_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
% This is stash 3-359; The height of the diagnostic parcel - I'm hoping
% this might give the top of the "2nd inversion" for cumulus boundary
% layers.

var_UM = 'BL_height_dthdz';    
um_case=um_case_PD; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
[BL_height_dthdz_ALL,nT] = UM_get_time_data_mat_nc(dat_global,time_inds,load_type,gcm_Plat2D_UM);
% BL height based on max height of the theta gradient since CAM doesn't
% have the same metrics for BL height.

tit_str_clean='';

%% Create the "combined" product where use the mixing height for all except the cumulus type 6 where use the parcel height.
BL_height_combined_ALL = BL_mixing_height_ALL;
i = find(BL_type_ALL == 6);
BL_height_combined_ALL(i) = BL_parcel_height_ALL(i);

%% Plot BL height for theta gradient technique
var_UM = 'BL height dthdz';
subtitle_str=var_UM;
add_str='';
dat_modis = meanNoNan(BL_height_dthdz_ALL,3);    
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
caxis([0 2500]);

%% Plot BL combined technique
var_UM = 'BL height from BL scheme';
subtitle_str=var_UM;
add_str='';
dat_modis = meanNoNan(BL_height_combined_ALL,3);    
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
caxis([0 2500]);


%% Plot fraction of time and BL height for all BL types

BL_plot_type = 'BL height combined';
BL_plot_type = 'BL height theta gradient';

for ibl=1:7
    var_UM = ['Fraction for BL type ' num2str(ibl)];
    i = find(BL_type_ALL == ibl); %Can treat as integers.
    inot = find(BL_type_ALL ~= ibl); %Can treat as integers.
    n_bl = zeros(size(BL_type_ALL));
    n_bl(i) = 1;
    n_all = ones(size(BL_type_ALL));
    
    dat_modis = sum(n_bl,3) ./ sum(n_all,3);
     
    %figure
    %UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    %caxis([0 1]);
    
    
    switch BL_plot_type
        case 'BL height combined'
            %BL_height_i = BL_height_ALL;
            %BL_height_i = BL_mixing_height_ALL;
            %BL_height_i = BL_parcel_height_ALL;
            BL_height_i = BL_height_combined_ALL;
        case 'BL height theta gradient'
            BL_height_i = BL_height_dthdz_ALL;
    end
    
    var_UM = [BL_plot_type ' for type ' num2str(ibl)];
    
    BL_height_i(inot) = NaN;
    dat_modis = meanNoNan(BL_height_i,3);
    
    subtitle_str=var_UM;
    add_str='';
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    caxis([0 2500]);
    
end




