%% Plot UKESM ensemble mean values for May for comparison with Leighton (used May 2015).

load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
dat_Leighton=load(load_file,'dat_ens_mean','years_ukesm_1d','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM'...
    , 'months_ukesm','time_ukesm');

%%

time_ukesm = dat_Leighton.time_ukesm'; %is ordered in the opposite way to as in dat_ens_mean...

%N.B. - data is not in the correct time order prior to around 1969, but
%after that it looks ok. Since we are plotting for the satellite era it
%should therefore be ok.

%year_match=2014;
month_match_str='May-2014';
day_match_str=['01-' month_match_str];
day_match=datenum(day_match_str);
it=find(time_ukesm==day_match);
%iyear = find(dat_Leighton.years_ukesm_1d == year_match);

dat_modis = squeeze(dat_Leighton.dat_ens_mean(it,:,:));

var_UM = ['UKESM SW TOA up ' month_match_str];
tit_str_clean = ['UKESM SW TOA up ' month_match_str];
%run plotting script
figure
ioverride_proj_type=1;
proj_type_DRIVER='ortho';
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis(clims);
caxis([0 250]);
%xlabel(hc,units_str_trend); %label the colour bar

%%
day_match_str=['15-' month_match_str];
day_match=datenum(day_match_str);
it=find(time_matlab_CERES==day_match);

dat_modis2 = squeeze(SW_TOA_ceres(it,:,:));
dat_modis =  griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);

var_UM = ['CERES SW TOA up ' month_match_str];
tit_str_clean = ['CERES SW TOA up ' month_match_str];

figure
ioverride_proj_type=1;
proj_type_DRIVER='ortho';
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 250]);
xlabel(hc,'W m^{-2}'); %label the colour bar


