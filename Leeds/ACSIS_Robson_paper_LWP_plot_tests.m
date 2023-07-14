
ioverride_LAT_plots=0;
ioverride_ticks_DRIVER=0;
irestrict_domain_DRIVER=1; %set to zero for a global map - otherwise shows the area specfied in UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global

%%
year=2010;


iyear = find(dat_ukesm.years_ukesm_1d == 2010);

var_UM = ['UKESM ' var_ukesm ' ' num2str(year) ' (g m^{-2})'];
tit_str_clean = var_UM; subtitle_str = var_UM;
icoarse_grain=0; M_coarse_grain=3; N_coarse_grain=3;
dat_modis = squeeze(dat_ukesm.dat_annual(iyear,:,:));
um_data = dat_modis;
sat_data = dat_modis; %for the stats (don't need sat ones).
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
caxis([0 150]);

%%
iyear = find(dat_ukesm_amip.years_ukesm_1d == 2010);

var_UM = ['UKESM-AMIP ' var_ukesm ' ' num2str(year) ' (g m^{-2})'];
tit_str_clean = var_UM; subtitle_str = var_UM;
icoarse_grain=0; M_coarse_grain=3; N_coarse_grain=3;
dat_modis = squeeze(dat_ukesm_amip.dat_annual(iyear,:,:));
um_data = dat_modis;
sat_data = dat_modis; %for the stats (don't need sat ones).
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
caxis([0 150]);

%% MAC
iyear = find(MAC_dat.years == 2010);

var_UM = ['MAC LWP ' num2str(year) ' (g m^{-2})'];
tit_str_clean = var_UM; subtitle_str = var_UM;
icoarse_grain=0; M_coarse_grain=3; N_coarse_grain=3;
dat_modis = squeeze(obs_annual_map_MAC(iyear,:,:));
um_data = dat_modis;
sat_data = dat_modis; %for the stats (don't need sat ones).
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
caxis([0 150]);