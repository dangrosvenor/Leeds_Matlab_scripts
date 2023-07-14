%iyear = find(years_sw_calc==1965);
%iy=iyear-5:iyear+5;
%iy=156:165;
%iy=1:11;
iy=2;

PI_lab = remove_character(UM_run_PI,'_',' ');
PD_lab = remove_character(UM_run_PD,'_',' ');


gcm_Plat2D_edges_UM = dat_PI.gcm_Plat2D_edges_UM;
gcm_Plon2D_edges_UM = dat_PI.gcm_Plon2D_edges_UM;
gcm_Plat2D_UM = dat_PI.gcm_Plat2D_UM;
gcm_Plon2D_UM = dat_PI.gcm_Plon2D_UM;



irestrict_domain_DRIVER=0;
ioverride_proj_type=0;
proj_type_DRIVER='ortho';
ioverride_LAT_plots=0;
iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;
icoarse_grain=0;
icontour_DRIVER=0; cont_col_str_DRIVER='';
time_round=0; time_format_str='';
isave_plot=0;
iplot_wind_arrows=0;



SW_calc_TOA = meanNoNan(SWTOA_calc_model_annual(iy,:,:),1);
SW_actual_TOA = meanNoNan(SWTOA_model_annual_PI(iy,:,:),1);

SW_calc_TOA_PD = meanNoNan(SWTOA_calc_model_annual_PD(iy,:,:),1);
SW_actual_TOA_PD = meanNoNan(SWTOA_model_annual_PD(iy,:,:),1);



%% SW clear sky

ioverride_LAT_plots=0;
dat_modis = meanNoNan(dat_PI.SW_up_TOA_clear_sky,1);
%var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
var_UM = 'SWTOA Clear-sky'; subtitle_str = [var_UM];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 250]);

%% CF

ioverride_LAT_plots=0;
dat_modis = meanNoNan(dat_PI.totCF,1);
%var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
var_UM = 'Tot CF'; subtitle_str = [var_UM];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 1]);

%% LWPic

ioverride_LAT_plots=0;
dat_modis = 1e3*meanNoNan(dat_PI.lwp_ic,1);
%var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
var_UM = 'LWPic (g m^{-3})'; subtitle_str = [var_UM];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);

%% Nd

ioverride_LAT_plots=0;
dat_modis = meanNoNan(dat_PI.Nd,1);
%var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
var_UM = 'Nd (cm^{-3})'; subtitle_str = [var_UM];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 200]);


%% SW biases

ioverride_LAT_plots=0;
dat_modis = SW_calc_TOA;
%var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' iy=' num2str(iy)];
var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 250]);

ioverride_LAT_plots=0;
dat_modis = SW_actual_TOA;
var_UM = 'Actual SWTOA'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 250]);

%Bias plot
ioverride_LAT_plots=0;
dat_modis = SW_calc_TOA - SW_actual_TOA;
var_UM = 'SWTOA Bias'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

% Percentage bias plot
ioverride_LAT_plots=0;
dat_modis = 100*(SW_calc_TOA - SW_actual_TOA)./SW_actual_TOA;
var_UM = 'SWTOA % Bias'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

%Bias plot
ioverride_LAT_plots=0;
dat_modis = SW_calc_TOA_PD - SW_actual_TOA_PD;
var_UM = 'SWTOA PD Bias'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

% Percentage bias plot
ioverride_LAT_plots=0;
dat_modis = 100*(SW_calc_TOA_PD - SW_actual_TOA_PD)./SW_actual_TOA_PD;
var_UM = 'SWTOA PD % Bias'; subtitle_str = [var_UM ' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-100 100]);

% No point plotting clear-sky comparison if using the model clear-sky as
% input since they are the same.

% ioverride_LAT_plots=0;
% dat_modis = meanNoNan(SWTOA_clear_sky_calc_model_annual(iy,:,:),1);
% var_UM = 'Calculated SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];
% %run plotting script
% figure
% UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
% lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
% caxis([0 75]);
% 
% ioverride_LAT_plots=0;
% dat_modis = meanNoNan(SWTOA_clear_sky_model_annual_PI(iy,:,:),1);
% var_UM = 'Actual SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];
% %run plotting script
% figure
% UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
% lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
% caxis([0 75]);

% ioverride_LAT_plots=0;
% dat_modis = meanNoNan(SWTOA_clear_sky_calc_model_annual(iy,:,:),1) - meanNoNan(sw_TOA_up_clear_dat_annual(iy,:,:),1);
% var_UM = 'Bias SWTOA clear-sky'; subtitle_str = [var_UM 'trans_clear_sky=' num2str(trans_clear_sky) ' m=' num2str(m) ' c=' num2str(c)];
% %run plotting script
% figure
% UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
% lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
% caxis([- 20 20]);

%                 ioverride_LAT_plots=0;
%                 dat_modis = meanNoNan(surf_albedo_annual(iy,:,:),1); var_UM = 'Annual mean surface albedo'; subtitle_str = [var_UM ' iy=' num2str(iy)];
%                 %run plotting script
%                 figure
%                 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                 caxis([0 0.2]);

%                 ioverride_LAT_plots=0;
%                 dat_modis = meanNoNan(trans_orig_ann(iy,:,:),1); var_UM = 'Annual mean transmissivity'; subtitle_str = [var_UM ' iy=' num2str(iy)];
%                 %run plotting script
%                 figure
%                 UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                 lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                 caxis([0.6 0.9]);

%montly plots
%

%ifind = find(years_sw_calc_monthly==1970);
%ims=ifind:ifind+11;

%
%                 for im=ims
%
%                     imon = im-ims(1)+1
%                     %
%                     %         ioverride_LAT_plots=0;
%                     %         dat_modis = squeeze(SWTOA_calc_model(im,:,:)); var_UM = 'Calculated SWTOA'; subtitle_str = [var_UM ',im=' num2str(im)];
%                     %         %run plotting script
%                     %         figure
%                     %         UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     %         lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     %         caxis([0 250]);
%                     %
%                     %         ioverride_LAT_plots=0;
%                     %         dat_modis = squeeze(SW_up_TOA_dat_ens_monthly(im,:,:)); var_UM = 'Actual SWTOA'; subtitle_str = [var_UM ',im=' num2str(im)];
%                     %         %run plotting script
%                     %         figure
%                     %         UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     %         lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     %         caxis([0 250]);
%
%                     ioverride_LAT_plots=0;
%                     dat_modis = meanNoNan(surf_albedo(im,:,:),1); var_UM = 'Monthly mean surface albedo'; subtitle_str = [var_UM ' im=' num2str(imon)];
%                     %run plotting script
%                     figure
%                     UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
%                     lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%                     caxis([0 0.2]);
%                     %
%                 end


%                 %For calculation of map of trends for checking against actual trend map
%                 dat_ukesm.years_ukesm_1d = years_sw_calc;
%                 dat_ukesm.dat_annual = SWTOA_calc_model_annual;
%                 MIP = 'DAMIP hist-aer SW calculated'; %Be careful with
%                 this as it screws up the MIP value needed in this
%                 calculation script...
%                 var_ukesm = 'SW_up_TOA';


            
            

%% Change in Nd

%year iy
ioverride_LAT_plots=0;
dat_modis = squeeze(dNd_annual(iy,:,:));
var_UM = ['Vertically averaged \DeltaN_d (' UM_run_PD ' minus ' PI_lab '; cm^{-3}), year=' num2str(iy)]; subtitle_str = [var_UM]; %' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

%average over both years
ioverride_LAT_plots=0;
dat_modis = meanNoNan(dNd_annual(:,:,:),1);
var_UM = ['Vertically averaged \DeltaN_d (' UM_run_PD ' minus ' PI_lab '; cm^{-3}), average of all years']; subtitle_str = [var_UM]; %' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-10 10]);

%% Forcing due to Nd change (instantaneous radiative forcing, IRF).

ioverride_LAT_plots=0;
dat_modis = squeeze(aci(iy,:,:));
var_UM = ['IRF_{\DeltaNd} (' UM_run_PD ' minus ' PI_lab '; W m^{-2}), year=' num2str(iy)]; subtitle_str = [var_UM]; %' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-1 1]);

ioverride_LAT_plots=0;
dat_modis = meanNoNan(aci(:,:,:),1);
var_UM = ['IRF_{\DeltaNd} (' UM_run_PD ' minus ' PI_lab '; W m^{-2}), average of all years']; subtitle_str = [var_UM]; %' tr_atmos_constant_val=' num2str(tr_atmos_constant_val)];
%run plotting script
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-1 1]);



