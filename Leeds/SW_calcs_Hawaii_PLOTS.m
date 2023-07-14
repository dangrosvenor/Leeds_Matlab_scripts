thresh_SO2(1) = 1e-5; thresh_SO2(2) = 1e15; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));

var_str = '\DeltaSW_{TOA} Estimated';
dat_SO2_PI = SW_PI_TOA; dat_SO2_PI(inan) = NaN;
dat_SO2_PD = SW_PD_TOA; dat_SO2_PD(inan) = NaN;

thresh_dNd = -1;
thresh_N = -1;
thresh_str='';

thresh_N = 0; %Overwritten below
icoarse_grain=0;
M_coarse_grain=10; N_coarse_grain=10;
var_UM = '';
%dat_modis = mean_dat_SO2_dNd_PD - mean_dat_SO2_dNd_PI;
%thresh_N = 5; dat_modis(N_mean_dat_SO2_dNd_PI<thresh_N) = NaN;
dat_modis = meanNoNan(dat_SO2_PD,3) - meanNoNan(dat_SO2_PI,3);
[mean_overall_est,N_overall std_dev_overall] = meanNoNan(dat_modis(:),1)
%thresh_N = 5; dat_modis(N_mean_d_dat<thresh_N) = NaN;
subtitle_str = ['Time mean ' var_str ' for SO_2>' num2str(thresh_SO2,'%1.0e') ', ' thresh_str ...
    ', \DeltaN_d>' num2str(thresh_dNd) ', N>=' num2str(thresh_N)];
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-150 150]);
icoarse_grain=0;

%%
thresh_SO2(1) = 1e-5; thresh_SO2(2) = 1e15; inan = find(dSO2 < thresh_SO2(1) | dSO2 >= thresh_SO2(2));

var_str = '\DeltaSW_{TOA} Actual';
dat_SO2_PI = SW_up_TOA_PI_ALL; dat_SO2_PI(inan) = NaN;
dat_SO2_PD = SW_up_TOA_PD_ALL; dat_SO2_PD(inan) = NaN;

thresh_dNd = -1;
thresh_N = -1;
thresh_str='';

thresh_N = 0; %Overwritten below
icoarse_grain=0;
M_coarse_grain=10; N_coarse_grain=10;
var_UM = '';
%dat_modis = mean_dat_SO2_dNd_PD - mean_dat_SO2_dNd_PI;
%thresh_N = 5; dat_modis(N_mean_dat_SO2_dNd_PI<thresh_N) = NaN;
dat_modis = meanNoNan(dat_SO2_PD,3) - meanNoNan(dat_SO2_PI,3);
[mean_overall_actual,N_overall std_dev_overall] = meanNoNan(dat_modis(:),1)
%thresh_N = 5; dat_modis(N_mean_d_dat<thresh_N) = NaN;
subtitle_str = ['Time mean ' var_str ' for SO_2>' num2str(thresh_SO2,'%1.0e') ', ' thresh_str ...
    ', \DeltaN_d>' num2str(thresh_dNd) ', N>=' num2str(thresh_N)];
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-150 150]);
icoarse_grain=0;




