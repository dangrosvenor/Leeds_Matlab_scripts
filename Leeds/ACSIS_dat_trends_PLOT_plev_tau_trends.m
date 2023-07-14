%% Choose region and land/ocean screening
box_region = '4';
box_region = '18'; %Northern NA as used in ACP paper
%box_region = '19'; %Southern NA as used in ACP paper
ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.

inew_folder=0;
if inew_folder==1
    savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
    eval(['!mkdir ' savedir_date]);
else
    savedir_date='/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
end

%%
region_str=box_region_str;

region_str=box_region_str;
var_UM='clisccp';

%Load dir
loadfile_pre_str = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_' region_str '_trend_ukesm'];
loadfile = [loadfile_pre_str '_' var_UM '.mat'];
dat_load = load(loadfile);

%From the NetCDF file :-
plevs = [90000, 74000, 62000, 50000, 37500, 24500, 9000]/100; 
plev_bnds = [100000 80000 68000 56000 44000 31000 18000 0]/100;
tau_bnds = [0 0.3 1.3 3.6 9.4 23 60 100000];
tau_levs = [0.15, 0.8, 2.45, 6.5, 16.2, 41.5, 100];




cf_plev = zeros([length(dat_load.dat_box_tau_plev{itau+1,iplev+1}) size(dat_load.dat_box_tau_plev,2)]);
cf_3d = NaN*ones([length(dat_load.dat_box_tau_plev{itau+1,iplev+1}) size(dat_load.dat_box_tau_plev,1) size(dat_load.dat_box_tau_plev,2)]);
for itau=0:6    
    for iplev=0:6                       
        dat_trends_tau_plev(itau+1,iplev+1) = dat_load.trends_tau_plev{itau+1,iplev+1}.coeffs(2) / 100; %convert from % to cloud fraction
        uncer_tau_plev(itau+1,iplev+1) = dat_load.trends_tau_plev{itau+1,iplev+1}.uncer_max / 100; %convert from % to cloud fraction
        % dat_load.dat_box_tau_plev{itau+1,iplev+1};
        cf_plev(:,iplev+1) = cf_plev(:,iplev+1) + dat_load.dat_box_tau_plev{itau+1,iplev+1}';
        cf_3d(:,itau+1,iplev+1) = dat_load.dat_box_tau_plev{itau+1,iplev+1}';
    end   
end

cf_time_mean_tau_plev = meanNoNan(cf_3d,1);
cf_time_mean_plev = meanNoNan(cf_plev,1);


clear tau_labs
for i=1:length(tau_bnds)
    tau_labs(i) = {num2str(tau_bnds(i))};
end


clear p_labs
for i=1:length(plev_bnds)
    p_labs(i) = {num2str(plev_bnds(i))};
end

%% Plot the cloud fraction trend

qpcolor(dat_trends_tau_plev'); %swap so that pressure is on the y-axis
set(gcf,'color','w');
set(gca,'xticklabels',tau_labs);
set(gca,'yticklabels',p_labs);
caxis([-13 13]*1e-5);
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
xlabel('Optical depth','fontsize',16);
ylabel('Pressure (hPa)','fontsize',16);
titlename = 'Cloud fraction trend 1985-2014 (yr^{-1})';
title(titlename,'fontsize',16);
set(gca,'fontsize',14);

savename=[savedir_date titlename ' ' region_str ' ' land_ocean];
clear opts
%opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% Plot uncertainty
qpcolor(uncer_tau_plev'); %swap so that pressure is on the y-axis
set(gcf,'color','w');
set(gca,'xticklabels',tau_labs);
set(gca,'yticklabels',p_labs);
caxis([0 5.5]*1e-5);
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
xlabel('Optical depth','fontsize',16);
ylabel('Pressure (hPa)','fontsize',16);
titlename='Uncertainty in cloud fraction trend 1985-2014 (yr^{-1})';
title(titlename,'fontsize',16);
set(gca,'fontsize',14);

savename=[savedir_date titlename ' ' land_ocean];
clear opts
%opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% Plot percent uncertainty
dat = abs(100*uncer_tau_plev'./dat_trends_tau_plev');
qpcolor(dat); %swap so that pressure is on the y-axis
set(gcf,'color','w');
set(gca,'xticklabels',tau_labs);
set(gca,'yticklabels',p_labs);
caxis([0 100]);
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
xlabel('Optical depth','fontsize',16);
ylabel('Pressure (hPa)','fontsize',16);
titlename='Percent uncertainty in cloud fraction trend 1985-2014 (%)';
title(titlename,'fontsize',16);
set(gca,'fontsize',14);

savename=[savedir_date titlename ' ' land_ocean];
clear opts
%opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% Plot overall CF occurrence
qpcolor(cf_time_mean_tau_plev'/100); %swap so that pressure is on the y-axis
set(gcf,'color','w');
set(gca,'xticklabels',tau_labs);
set(gca,'yticklabels',p_labs);
caxis([0 7]/100);
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
xlabel('Optical depth','fontsize',16);
ylabel('Pressure (hPa)','fontsize',16);
titlename = 'Time Average Cloud Fraction 1985-2014';
title(titlename,'fontsize',16);
set(gca,'fontsize',14);

savename=[savedir_date titlename ' ' land_ocean];
clear opts
%opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

%% Plot profile of time mean CF across all tau
figure; set(gcf,'color','w');
h=plot(cf_time_mean_plev/100,plevs,'bo-','linewidth',3,'markerfacecolor','b');

set(gca,'ydir','reverse');
set(gca,'fontsize',14);
grid on
xlabel('Cloud Fraction','fontsize',16);
ylabel('Pressure (hPa)','fontsize',16);
titlename = 'Cloud Fraction for all \tau';
title(titlename);

savename=[savedir_date titlename ' ' land_ocean];
clear opts
%opts.iplot_png=1;
opts.iplot_eps=1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)


