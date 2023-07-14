%% UKESM eval
%Need to run 
%ACSIS_Robson_paper_load_data  - sets var_ukesm in there
%before this to load in the data. Also sets savedir_date.

% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

%% Choose region
% Box for regional averages
%box_region_DRIVER = '1'; %US small box
%box_region_DRIVER = '3'; %SW of Spain ('Europe outflow')
%box_region_DRIVER = '4'; %All NA region
%box_region_DRIVER = '5'; %UK region
%box_region_DRIVER = '6'; %West of UK region
%box_region_DRIVER = '7'; %20 deg further west of there
%box_region_DRIVER = '8'; %middle of Atlantic at US box latitude
%box_region_DRIVER = '9'; %20 deg further west of there

box_region = box_region_DRIVER;

% ---
ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.
% ---

iscreen_land=1;
tail_test=2; %1 or 2-tailed test for signficance
iplot_CERES=0;

%%



tit_str_clean='';
i_plot_all_boxes=1;

% Load CF data using :-
%read_calipso_monthly_IPSL_2007_2017

% var_ukesm = 'Nd_cf_weighted_UKESM';
% var_ukesm = 'calipso_low_cloud_amount';

yr_start_trend=1850;        
%yr_start_trend=2003;
yr_end_trend=2014;

iscreen_sig=1; %Whether to screen for signficance
marker_size=1; %For non-signficant points
%marker_size=1; %For non-signficant points
iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid

p_conf = 95; % Confidence limit (%) for the trend significance
nthresh_days = 3;
%nthresh_days = 0;




UKESM_Nd_case = 'to ztop';
UKESM_Nd_case = 'to 3.2km';



icoarse_grain=0;

time_round='';
time_format_str='';
icontour_DRIVER=0;
isave_plot=0;
iplot_wind_arrows=0;

cont_col_str_DRIVER='k';

% This is the region to be plotted as a map (whole map region, not the box
% region).
region_choice = 'Southern NA';
%region_choice = 'Northern NA';
region_choice = 'Rosenfeld VOCALS';
region_choice = 'Rosenfeld ALL';
%region_choice = 'VOCALS CPT';
%region_choice = 'VOCALS coastal';
region_choice = 'Northern Hemisphere Yawen';
region_choice = 'Northern Hemisphere Yawen, further north';
region_choice = 'North Atlantic further north';
%region_choice = 'global';

irestrict_domain_DRIVER=1;
switch region_choice
    case 'global'
        irestrict_domain_DRIVER=0;
end

% Run the script to pick the region
[LAT_val_DRIVER2, LON_val_DRIVER2, region_shortname] = UM_ACSIS_choose_region(region_choice);

LAT_val_DRIVER_override = LAT_val_DRIVER2;
LON_val_DRIVER_override = LON_val_DRIVER2;

% var_UM = 'Nd_cf_weighted_UKESM_ztop';
% um_case='UKESM/r1i1p1f2_u-bc179/'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %
% dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
% dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);


%% Do the box mean, etc. calculations


%Run script to choose the data for the selected season and then to average
%over the box region.
ACSIS_Robson_paper_box_means

% %For the shorter timeseries (calculated SW)
% [dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std,trend_dat_box,...
%     trend_dat_box_ens,trend_ens_dat_box, dat_annual_box_ukesm_ens] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
%     yr_start_trend_box,yr_end_trend_box,p_conf,ibox);
% 
% %For the full ukesm (actual SW) timeseries
% %[dat_annual_box_ukesm_actual,dat_annual_box_ukesm_actual_ens_std,trend_dat_box_actual,...
% %    trend_dat_box_ens_actual,trend_ens_dat_box_actual, dat_annual_box_ukesm_ens_actual] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm_save,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
% %    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);

%% Obs data - calc regional means, etc.

inc_CERES = 0;
if inc_CERES==1
% CERES
obs_str='CERES';
years_obs = years_obs_CERES;

if iscreen_land==1
    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(SW_TOA_ceres,gcm_Plat2D_CERES,gcm_Plon2D_CERES);
end

[obs_monthly_box,obs_annual_box,obs_annual_map] =...
    ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES,season);
end

% Also do for DEEP-C SW dataset for SW
obs_str2='DEEP-C';
years_obs2 = [1985:2014];

if iscreen_land==1
    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_deepc.sw_up_toa,gcm_Plat2D_UM,gcm_Plon2D_UM);
end

[obs_monthly_box2,obs_annual_box2,obs_annual_map2] ... 
    = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_val,LON_val,years_obs2,gcm_area_UM,season);


%% Calculate the OBS CERES trend for the boxed region
for it_trend=1:1 %length(yr_start_trend_box_obs)
    yr_start=years_obs(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
    switch season
        case 'DJF'
            yr_start = yr_start + 1; %start a year later for now - although could do DJF of Y2000 for CERES
    end
    yr_start_trend_used_box_obs=yr_start; yr_end_trend_used_box_obs=yr_end;       
    
    %trend_dat_box_obs = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box,years_obs,yr_start,yr_end,p_conf,it_trend,ibox);    
    [trend_dat_box_obs{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs, obs_annual_box , p_conf);
end

%% Calculate the OBS DEEP-C trend for the boxed region
for it_trend=1:1 %length(yr_start_trend_box_obs)
    yr_start=years_obs2(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);   
    switch season
        case 'DJF'
            yr_start = yr_start + 1; %start a year later because data starts in Jan 1985
    end
    yr_start_trend_used_box_obs2=yr_start; yr_end_trend_used_box_obs2=yr_end;        
    
    %trend_dat_box_obs2 = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box2,years_obs2,yr_start,yr_end,p_conf,it_trend,ibox);
    [trend_dat_box_obs2{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs2, obs_annual_box2 , p_conf);
end   

%% plot annual mean timeseries of calculated model SW vs actual
% Could put this UKESM code in a seprate script to use in modular form
%Full record is 1850 to 2014. We want 1984 to end.
%istart=1984-1850+1;
istart=yr_start_plot-1850+1;

titlenam = [season ' mean ' var_str  ' for region ' box_region];
tit_str_clean = [season ', ' box_region_str];
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
increase_font_size_map_figures
%set(gcf,'position',[3         297        1256         322]);
set(gcf,'position',[3         297        1256         422]);
%set(gca,'position',[0.1300    0.1500    0.7750    0.7150]);
set(gca,'position',[0.1300    0.100    0.7750    0.7150]);
clear leg_str
ileg=1;

y_ens_upper = dat_annual_box_ukesm(istart:end) + dat_annual_box_ukesm_ens_std(istart:end)*2;
y_ens_lower = dat_annual_box_ukesm(istart:end) - dat_annual_box_ukesm_ens_std(istart:end)*2;
x_ens = dat_ukesm.years_ukesm_1d(istart:end)';
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
%h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
%p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
%Did the above to avoid weird legend behaviour - didn't quite work
hold on
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'bo-');

fscale=1e2; fscale_str=['\times10^{-2}'];
istart_trend=yr_start_trend_box-1850+1;
if iplot_individual_ens==1
   
    ens_str = ['iens=' num2str(iens_plot)];
    %h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_ens_actual(istart:end,iens_plot),'b-');
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
    %trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens_actual{ibox,1,iens_plot}.coeffs(2),'%.3f')];
    %trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),'%.3f')];
    trend_str = [ens_str ', trend=(' num2str(fscale*trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_ens{ibox,1,iens_plot}.uncer_max,'%.2f') ') ' fscale_str];
        
    
    model_dat_trend_period = dat_annual_box_ukesm_ens(istart_trend:end,iens_plot);
else
    clear corr_timser
    iens_plot=1;
    
    ens_str = 'ensemble mean';
    %h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_actual(istart:end),'b-');
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm(istart:end),'b-');    
    %trend_str=['trend=' num2str(trend_dat_box_actual{ibox,1}.coeffs(2),'%.3f')];
    %trend_str=['trend=' num2str(trend_dat_box{ibox,1}.coeffs(2),'%.3f')];
    trend_str=[ens_str ', trend=(' num2str(fscale*trend_dat_box{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];    
    
    model_dat_trend_period = dat_annual_box_ukesm(istart_trend:end);
end
%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}=['UKESM1, ' trend_str]; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on

if iplot_CERES==1

%CERES
if length(years_obs)>1
    h2=plot(years_obs,obs_annual_box,'r-');
else
    h2=plot(years_obs,obs_annual_box,'rs');
end
trend_str=['trend=(' num2str(fscale*trend_dat_box_obs{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
leg_str{ileg}=[obs_str ', ' trend_str]; ileg=ileg+1;
set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');

end

%Deep-C
if length(years_obs)>1
    h2=plot(years_obs2,obs_annual_box2,'k-');
else
    h2=plot(years_obs2,obs_annual_box2,'ks');
end
trend_str=['trend=(' num2str(fscale*trend_dat_box_obs2{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs2{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
leg_str{ileg}=[obs_str2 ', ' trend_str]; ileg=ileg+1;

set(h2,'linewidth',4);
set(h2,'markerfacecolor','k')


loc='SouthWest';
loc='NorthEast';
loc='NorthWest';

%set(gca,'xlim',[1980 2014]);
%set(gca,'xlim',[1844 2018]);
set(gca,'xlim',[yr_start_plot-2 2018]);

%Correlation coefficient over time of the model vs obs for the trend
%period.
switch season
    case 'DJF'        
        obs_dat = obs_annual_box2(2:end);
    otherwise
        obs_dat = obs_annual_box2;
end

corr_timser(iens_plot) = corr(model_dat_trend_period(:),obs_dat(:));

%% Plot the linear trend lines

if iplot_lin==1
   conf_str_tit='';
   
   %UKESM SW
   for it_trend=1:length(yr_start_trend_box)
       if iplot_trend(it_trend)==1
           if iplot_individual_ens==1              
               %x_trend = trend_dat_box_ens_actual{ibox,it_trend,iens_plot}.x;
               %y_trend = trend_dat_box_ens_actual{ibox,it_trend,iens_plot}.ylin;
               x_trend = trend_dat_box_ens{ibox,it_trend,iens_plot}.x;
               y_trend = trend_dat_box_ens{ibox,it_trend,iens_plot}.ylin;
           else
               %x_trend = trend_dat_box_actual{ibox,it_trend}.x;
               %y_trend = trend_dat_box_actual{ibox,it_trend}.ylin;
               x_trend = trend_dat_box{ibox,it_trend}.x;
               y_trend = trend_dat_box{ibox,it_trend}.ylin;
           end
           plot( x_trend , y_trend , 'b-','linewidth',3,'HandleVisibility','Off');
       end
      yr_start_trend_box_str = num2str(yr_start_trend_box(it_trend));
      yr_end_trend_box_str = num2str(yr_end_trend_box(it_trend));
      switch tail_test
          case 1
              if iplot_individual_ens==1
                  %conf_str = num2str(100*(1 - trend_dat_box_ens_actual{ibox,it_trend,iens_plot}.T1),'%1.1f');
                  conf_str = num2str(100*(1 - trend_dat_box_ens{ibox,it_trend,iens_plot}.T1_max),'%1.1f');
              else
                  %conf_str = num2str(100*(1 - trend_dat_box_actual{ibox,it_trend}.T1),'%1.1f');
                  conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T1_max),'%1.1f');
              end
          case 2     
              if iplot_individual_ens==1
                  %conf_str = num2str(100*(1 - trend_dat_box_ens_actual{ibox,it_trend,ens_plot}.T2),'%1.1f');
                  conf_str = num2str(100*(1 - trend_dat_box_ens{ibox,it_trend,iens_plot}.T2_max),'%1.1f');
              else
                  %conf_str = num2str(100*(1 - trend_dat_box_actual{ibox,it_trend}.T2),'%1.1f');
                  conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T2_max),'%1.1f');
              end
      end            
      
      conf_str_tit = [conf_str_tit '; p_{UKESM1}=' conf_str '%'];
      %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
   end
   
   
   if iplot_CERES==1
       
   %plot CERES trend
   it_trend=1; %just one trend for obs since timeseries is short
   plot( trend_dat_box_obs{ibox,it_trend}.x , trend_dat_box_obs{ibox,it_trend}.ylin , 'r-','linewidth',3,'HandleVisibility','Off');   
   switch tail_test
       case 1
           conf_str = num2str(100*(1 - trend_dat_box_obs{ibox,it_trend}.T1_max),'%1.1f');
       case 2
           conf_str = num2str(100*(1 - trend_dat_box_obs{ibox,it_trend}.T2_max),'%1.1f');
   end
   conf_str_tit = [conf_str_tit '; p_{Ceres}=' conf_str '%'];
   
   end
   
   
   %Deep-C trend
   it_trend=1; %just one trend for obs since timeseries is short
   plot( trend_dat_box_obs2{ibox,it_trend}.x , trend_dat_box_obs2{ibox,it_trend}.ylin , 'k-','linewidth',3,'HandleVisibility','Off');   
   switch tail_test
       case 1
           conf_str = num2str(100*(1 - trend_dat_box_obs2{ibox,it_trend}.T1_max),'%1.1f');
       case 2
           conf_str = num2str(100*(1 - trend_dat_box_obs2{ibox,it_trend}.T2_max),'%1.1f');
   end
   conf_str_tit = [conf_str_tit '; p_{Deep-C}=' conf_str '%'];
   
   tit_str_clean = [tit_str_clean conf_str_tit];
   

end



% -- plot the error bar for the PI --
%errorbar(1980,me_t_PI,std_t_PI*2,'CapSize',2,'ro','markerfacecolor','r');
herr=errorbarYY('vert',yr_PI_bar,me_t_PI,std_t_PI*2,gca,'b','o',2,0.01);
set(herr,'HandleVisibility','Off');
%errorbarYY('vert',1960,me_t_PI-std_t_PI,std_t_PI,gca,'r','o',2,0.01);

hpi_mean = plot(yr_PI_bar,me_t_PI,'ko','markerfacecolor','b','markersize',10);
set(hpi_mean,'HandleVisibility','Off');

        
switch var_ukesm
    case 'Nd_cf_weighted_UKESM'
        
        
        set(gca,'ylim',ylims_main_Nd);
        %title(titlenam,'position',[1925 200]);
        %title(tit_str_clean,'position',[1945 200]);
        if iplot_legend==1
            legend(leg_str,'location',loc);
        end
        
    case 'calipso_low_cloud_amount'
        
        
        set(gca,'ylim',ylims_main_cf);
        if iplot_legend_cf==1
            legend(leg_str,'location',loc);
        end
        
    case 'SW_up_TOA'
        
        set(gca,'ylim',ylims_main_sw);
        
        if iplot_legend_sw==1
            legend(leg_str,'location',loc);
        end
        
end




ylims=get(gca,'ylim'); 
%title(tit_str_clean,'position',[1995 ylims(2)]);
tit_str_clean = [tit_str_clean ', r_{corr}=' num2str(corr_timser(iens_plot),'%.2f')];
tit2 = {['Region=' box_region],'', tit_str_clean};
title(tit_str_clean,'position',[yr_title,ylims(2)]);
title(tit2,'position',[yr_title,ylims(2)]);

xlabel('Year');
ylabel([var_str ' (' units_str ')']);

%switch box_region_DRIVER
%    case '1'
%        title('');
%        title(titlenam,'position',[1925 200]);
%    otherwise
%        title(titlenam,'position',[1925 200]);
%end


grid on


savename=[savedir_date titlenam ' ' ens_str];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%Calculate the bias between model and obs for each year.

model_interp = interp1(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,years_obs);
prc_bias_annual = 100*(model_interp./obs_annual_box - 1);
[stats_regions{ibox}]=calc_stats(1e-10,obs_annual_box,model_interp);
stats_regions{ibox}.box_region_str = box_region_str;

%% plot INSET annual mean timeseries of model vs obs for 2003 to 2014
iplot_inset=0;

if iplot_inset==1
titlenam = ['Annual mean ' var_str  ' for region ' box_region];
tit_str_clean = ['Region ' box_region_str];
%figure

%ax2 = axes('Position',[0.17 0.7 0.25 0.25]);
%ax2 = axes('Position',[0.1748    0.7683    0.2500    0.2500]);

switch var_ukesm
    case 'Nd_cf_weighted_UKESM'
        ax2 = axes('Position',inset_ax_pos);    
        xlims_inset = [2002 2015];
    case 'calipso_low_cloud_amount'
        ax2 = axes('Position',inset_ax_pos_cf); 
        xlims_inset = [2006 2015];
    case 'SW_up_TOA'  
        ax2 = axes('Position',inset_ax_pos_sw); 
        xlims_inset = [2000 2015];
end

box on
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%increase_font_size_map_figures
set(gca,'fontsize',16);
%set(gcf,'position',[3         297        1256         322]);
clear leg_str
ileg=1;
%h1=errorbar(dat.years_ukesm_1d,dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std*2,'bo-');
y_ens_upper = dat_annual_box_ukesm + dat_annual_box_ukesm_ens_std*2;
y_ens_lower = dat_annual_box_ukesm - dat_annual_box_ukesm_ens_std*2;
x_ens = dat_ukesm.years_ukesm_1d';
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
hold on
h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,'bo-');
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'b-');

%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}='UKESM'; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on


h2=plot(years_obs,obs_annual_box,'rs-');
%leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
%leg_str{ileg}=['CALIPSO'];
leg_str{ileg}=obs_str;
set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');
loc='SouthWest';
%loc='NorthEast';

set(gca,'xlim',xlims_inset);
%set(gca,'xlim',[1850 2018]);

%set(gca,'ylim',[100 175]);
%set(gca,'ylim',[45 210]);
%set(gca,'ylim',[60 170]);
switch var_ukesm
    case 'Nd_cf_weighted_UKESM'        
        set(gca,'ylim',ylims_inset);
    case 'calipso_low_cloud_amount'
        set(gca,'ylim',ylims_inset_cf);
        %set(gca,'ylim',[0.3 0.45]);
    case 'SW_up_TOA'    
        set(gca,'ylim',ylims_inset_sw);
end



%legend(leg_str,'location',loc);
%xlabel('Year');
%ylabel([var_str ' ' units_str]);
title(titlenam);
grid on


savename=[savedir_date 'inset ' titlenam ' ' ens_str];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
title('');
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%title(titlenam);
    
end


%% Do box and whisker plot
if iplot_individual_ens~=1
    ACSIS_Robson_paper_plot_ens_trend_Box_Whiskier_generic
end
