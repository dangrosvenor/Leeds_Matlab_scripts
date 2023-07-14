clear region_choice_multi region_biases
i=1;

region_choice_multi{i} = 'Global'; i=i+1;
region_choice_multi{i} = 'North Atlantic'; i=i+1;
region_choice_multi{i} = 'Southern NA'; i=i+1;
region_choice_multi{i} = 'Northern NA'; i=i+1;
%region_choice_multi{i} = 'Rosenfeld VOCALS'; i=i+1;
%region_choice_multi{i} = 'Rosenfeld ALL'; i=i+1;
region_choice_multi{i} = 'VOCALS CPT'; i=i+1;
region_choice_multi{i} = 'VOCALS coastal'; i=i+1;
region_choice_multi{i} = 'East of US SW trend'; i=i+1;
region_choice_multi{i} = 'North Atlantic SW paper'; i=i+1;


land_ocean = 'land and ocean';
land_ocean = 'ocean only';
%land_ocean = 'land only';

load_type = 'merged netCDF';
var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);

land_mask_mean = squeeze(ones(size(um_data)));
switch land_ocean
    case 'land and ocean'
        %don't need to do anything
    case 'ocean only'
        land_mask_mean(dat_global.dat==1)=NaN;
    case 'land only'
        land_mask_mean(dat_global.dat==0)=NaN;
end

um_data2 = um_data.*land_mask_mean;


[gcm_area_UM] = calc_area_lat_lon2d(gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
%gcm_area_UM(:)=1;

for i=1:length(region_choice_multi)       
    
    % Run the script to pick the region
    [LAT_val_DRIVER2, LON_val_DRIVER2, region_shortname] = UM_ACSIS_choose_region(region_choice_multi{i});    
    
    if iplot_box_on_map==1 & i>1 %not plotting the box for the full NA region
        UM_ACSIS_run_plot_box_commands
    end
    
%weight by gridbox area
    region_biases{i} = calc_biases_for_regions(0.1,sat_data,um_data2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_area_UM);
    region_biases{i}.region_name = region_choice_multi{i};
    region_biases{i}.region_shortname = region_shortname;
end

iappend=0;
for i=1:length(region_biases)
    savename2 = [savename land_ocean];
    savename2 = remove_character(savename2,' ','_');
    savename2 = remove_character(savename2,'%','_');
    savename2 = remove_character(savename2,':','_');
    %This saves the varaible to the file (savename2)
    
    latex_newcommand_from_structure(region_biases{i},[region_biases{i}.region_shortname var_Latex],savename2,iappend);
    iappend=1;
end