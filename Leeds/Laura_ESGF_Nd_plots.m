%%

load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_ukesm_all_' var_ukesm '.mat'];
expt_str2 = '';
dat = load(load_file,'dat_ens_mean');
eval_str2 = ['[me'  expt_str2 ',N' expt_str2 ',std' expt_str2 '] = meanNoNan(dat.dat_ens_mean,1);'];
eval(eval_str2);
    
%%

models={'','_GFDL_ESM4','_CESM2'}; % '' = UKESM
model_titles={'UKESM','GFDL ESM4','CESM2'};

for im=1:length(models)
im
eval_str = ['dat_modis = N' models{im} ';']; eval(eval_str);
var_UM='';

eval_str = ['gcm_Plat2D_UM = dat_ukesm' models{im} '.gcm_Plat2D_UM;']; eval(eval_str);
eval_str = ['gcm_Plon2D_UM = dat_ukesm' models{im} '.gcm_Plon2D_UM;']; eval(eval_str);
eval_str = ['gcm_Plat2D_edges_UM = dat_ukesm' models{im} '.gcm_Plat2D_edges_UM;']; eval(eval_str);
eval_str = ['gcm_Plon2D_edges_UM = dat_ukesm' models{im} '.gcm_Plon2D_edges_UM;']; eval(eval_str);




proj_type_DRIVER='ortho'; plot_region_str='NA'; %Spherical globe projection - "angle of view" is chosen in plot_global_maps at present 
proj_type_DRIVER='other'; plot_region_str='global'; %Full global map in miller projection

LAT_val_DRIVER_override = [-1e9 1e9]; LON_val_DRIVER_override = [-1e9 1e9];


  %var_UM = [model_str_map ' ' var_str ' trend of ensemble mean between ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and ' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
    %tit_str_clean = ['UKESM, iens=' num2str(iens) ' ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];
    %subtitle_str = tit_str_clean;
    %add_str = [' ' units_str_trend];
    %subtitle_str=num2str(iens); add_str='';
    subtitle_str=''; add_str='';
    clims=[0 1980];
    units_str = 'Number of months';
    
    %run plotting script
    figure
    ioverride_proj_type=1;
    ioverride_LAT_plots=0;
    %proj_type_DRIVER='ortho'; %set at top of script.
    %proj_type_DRIVER='other';
    irestrict_domain_DRIVER=0;
    igeneric_plot=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims);
    %increase_font_size_map_figures;   %This creates a gap between map and
        %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?
    %caxis([-0.3 0.3]);
    xlabel(hc,units_str); %label the colour bar
    title(model_titles{im});
    
    
end