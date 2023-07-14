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
ACSIS_Robson_paper_choose_regional_box %run script - also chooses ylims, etc.
% ---

iscreen_land=1;
tail_test=2; %1 or 2-tailed test for signficance

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


%var_str = 'Total cloud fraction'; %set in
%ACSIS_Robson_paper_choose_clims_etc, which is run from ACSIS_Robson_paper_load_data.m

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


%% Do the calculations


%Run script to choose the data for the selected season and then to average
%over the box region.
ACSIS_Robson_paper_box_means




%bint; %The 95% confidence intervals of the coeff values (y-int and slope)
%stats (4 element vector); %Stats contains R^2, the F-statistic and it's p-value,
        %and an estimate of its error variance.
% The p-value for the F-statistic looks like it is the same as the T2 value (2-tailed T-test t value)        
        
%% Obs data - calc regional means, etc.

    
    % CALIPSO obs


                obs_str_cal='CALIPSO';
                                
                %Load Calipso data - currently 2007-2017
                % Load CF data using :-
                %read_calipso_monthly_IPSL_2007_2017
                years_obs_cal = years_requested; %from read_calipso_monthly_IPSL_2007_2017
                
                dat_obs = 0.01*cltcalipso_monthly_AVERAGE;
                if iscreen_land==1
                    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);
                end
                
                [obs_monthly_box_cal,obs_annual_box_cal,obs_annual_map_cal] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,LAT_val,LON_val,years_obs_cal,gcm_area_CALIPSO_monthly,season);
                
          %CCI dataset
                obs_str_cci='ESA CCI';
                                
                years_obs_cci = unique(cci_dat.Y_out);
                
                if ~exist('cfc_ocean')
                    dat_obs = cci_dat.cfc; %total CF (day and night)
                    if iscreen_land==1
                        [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs,gcm_Plat2D_CCI,gcm_Plon2D_CCI);
                    end
                    cfc_ocean = dat_tmp;
                else
                    dat_tmp = cfc_ocean;
                end
                
                [obs_monthly_box_cci,obs_annual_box_cci,obs_annual_map_cci] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CCI,gcm_Plon2D_CCI,LAT_val,LON_val,years_obs_cci,gcm_area_CCI,season);
                
                
                
                
        


%% Calculate the OBS CALIPSO trend for the boxed region
for it_trend=1:1 %length(yr_start_trend_box_obs)
    yr_start=years_obs_cal(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
    yr_start_trend_used_box_obs=yr_start; yr_end_trend_used_box_obs=yr_end;    
    
    %trend_dat_box_obs_cal = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box_cal,years_obs_cal,yr_start,yr_end,p_conf,it_trend,ibox);
    
    [trend_dat_box_obs_cal{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs_cal, obs_annual_box_cal , p_conf);
            
            
    
end

%% Calculate the OBS CCI trend for the boxed region
for it_trend=1:1 %length(yr_start_trend_box_obs)
    yr_start=years_obs_cci(1); 
    yr_start = 1985;
    %yr_start=years_obs_cal(1);
    yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);   
    yr_start_trend_used_box_obs2=yr_start; yr_end_trend_used_box_obs2=yr_end;    
    
    %trend_dat_box_obs_cci= ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box_cci,years_obs_cci,yr_start,yr_end,p_conf,it_trend,ibox);
   [trend_dat_box_obs_cci{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs_cci, obs_annual_box_cci , p_conf);
end   

%% plot annual mean timeseries of calculated model SW vs actual
%Full record is 1850 to 2014. We want 1984 to end.
istart=yr_start_trend_box-1850+1;
istart=yr_start_plot-1850+1;
%istart=1850-1850+1;

titlenam = ['Annual mean ' var_str  ' for region ' box_region];
tit_str_clean = [box_region_str];
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
increase_font_size_map_figures
set(gcf,'position',[3         297        1256         322]);
set(gca,'position',[0.1300    0.1500    0.7750    0.7150]);
clear leg_str
ileg=1;
%h1=errorbar(dat.years_ukesm_1d,dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std*2,'bo-');
y_ens_upper = dat_annual_box_ukesm(istart:end) + dat_annual_box_ukesm_ens_std(istart:end)*2;
y_ens_lower = dat_annual_box_ukesm(istart:end) - dat_annual_box_ukesm_ens_std(istart:end)*2;
x_ens = dat_ukesm.years_ukesm_1d(istart:end)';
inan = find(isnan(y_ens_upper)==1);
%Deal with NaNs since they stop patch from working
y_ens_upper(inan)=[];
y_ens_lower(inan)=[];
x_ens(inan)=[];
%
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
hold on
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'bo-');

fscale=1e4; fscale_str=['\times10^{-4}'];
if iplot_individual_ens==1
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
    trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),'%.2e') '\pm' num2str(trend_dat_box_ens{ibox,1,iens_plot}.uncer_max)];
else
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm(istart:end),'b-');
    %trend_str=['trend=' num2str(trend_dat_box{ibox,1}.coeffs(2),'%.2e') '\pm' num2str(trend_dat_box{ibox,1}.uncer2)];
    trend_str=['trend=(' num2str(fscale*trend_dat_box{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
end
%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}=['UKESM1, ' trend_str]; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on

%CALIPSO
if length(years_obs)>1
    h2=plot(years_obs_cal,obs_annual_box_cal,'r-');
else
    h2=plot(years_obs_cal,obs_annual_box_cal,'rs');
end
trend_str=['trend=(' num2str(fscale*trend_dat_box_obs_cal{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs_cal{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
%leg_str{ileg}=[obs_str_cal ', trend=' num2str(trend_dat_box_obs_cal{ibox,1}.coeffs(2),'%.2e')]; ileg=ileg+1;
leg_str{ileg}=[obs_str_cal ', ' trend_str]; ileg=ileg+1;

set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');

%CCI
if length(years_obs)>1
    h2=plot(years_obs_cci,obs_annual_box_cci,'k-');
else
    h2=plot(years_obs_cci,obs_annual_box_cci,'ks');
end
trend_str=['trend=(' num2str(fscale*trend_dat_box_obs_cci{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs_cci{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
%leg_str{ileg}=[obs_str_cci ', trend=' num2str(trend_dat_box_obs_cci{ibox,1}.coeffs(2),'%.2e')]; ileg=ileg+1;
leg_str{ileg}=[obs_str_cci ', ' trend_str];
set(h2,'linewidth',4);
set(h2,'markerfacecolor','k')


loc='SouthWest';
loc='NorthEast';
loc='NorthWest';

%set(gca,'xlim',[1980 2014]);
%set(gca,'xlim',[1844 2018]);
%set(gca,'xlim',[1850 2014]);
%set(gca,'xlim',[yr_start_trend_box 2014]);
set(gca,'xlim',[yr_start_plot-2 2016]);

%% Plot the linear trend lines

if iplot_lin==1
   conf_str_tit='';
   
   %UKESM SW
   for it_trend=1:length(yr_start_trend_box)
       if iplot_trend(it_trend)==1
           if iplot_individual_ens==1              
               x_trend = trend_dat_box_ens{ibox,it_trend,iens_plot}.x;
               y_trend = trend_dat_box_ens{ibox,it_trend,iens_plot}.ylin;
           else
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
                  conf_str = num2str(100*(1 - trend_dat_box_ens{ibox,it_trend,iens_plot}.T1),'%1.1f');
              else
                  conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T1),'%1.1f');
              end
          case 2     
              if iplot_individual_ens==1
                  conf_str = num2str(100*(1 - trend_dat_box_ens{ibox,it_trend,iens_plot}.T2),'%1.1f');
              else
                  conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T2),'%1.1f');
              end
      end            
      
      conf_str_tit = [conf_str_tit '; p_{UKESM1}=' conf_str '%'];
      %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
   end
   
   %plot CALIPSO trend
   it_trend=1; %just one trend for obs since timeseries is short
   plot( trend_dat_box_obs_cal{ibox,it_trend}.x , trend_dat_box_obs_cal{ibox,it_trend}.ylin , 'r-','linewidth',3,'HandleVisibility','Off');   
   switch tail_test
       case 1
           conf_str = num2str(100*(1 - trend_dat_box_obs_cal{ibox,it_trend}.T1),'%1.1f');
       case 2
           conf_str = num2str(100*(1 - trend_dat_box_obs_cal{ibox,it_trend}.T2),'%1.1f');
   end
   conf_str_tit = [conf_str_tit '; p_{Calipso}=' conf_str '%'];
   
   
   
   %CCI trend
   it_trend=1; %just one trend for obs since timeseries is short
   plot( trend_dat_box_obs_cci{ibox,it_trend}.x , trend_dat_box_obs_cci{ibox,it_trend}.ylin , 'k-','linewidth',3,'HandleVisibility','Off');   
   switch tail_test
       case 1
           conf_str = num2str(100*(1 - trend_dat_box_obs_cci{ibox,it_trend}.T1),'%1.1f');
       case 2
           conf_str = num2str(100*(1 - trend_dat_box_obs_cci{ibox,it_trend}.T2),'%1.1f');
   end
   conf_str_tit = [conf_str_tit '; p_{CCI}=' conf_str '%'];
   
   tit_str_clean = [tit_str_clean conf_str_tit];
   

end



% -- plot the error bar for the PI --
%errorbar(1980,me_t_PI,std_t_PI*2,'CapSize',2,'ro','markerfacecolor','r');
herr=errorbarYY('vert',1847.5,me_t_PI,std_t_PI*2,gca,'k','o',2,0.01);
set(herr,'HandleVisibility','Off');
%errorbarYY('vert',1960,me_t_PI-std_t_PI,std_t_PI,gca,'r','o',2,0.01);



switch var_ukesm
    case 'Nd_cf_weighted_UKESM'
        set(gca,'ylim',ylims_main_Nd);
        %title(titlenam,'position',[1925 200]);
        %title(tit_str_clean,'position',[1945 200]);
        if iplot_legend==1
            legend(leg_str,'location',loc);
        end
    case 'calipso_low_cloud_amount'
        plot(1847.5,me_t_PI,'ko','markerfacecolor','k','markersize',10);
        set(gca,'ylim',ylims_main_cf);        
        if iplot_legend_cf==1
            legend(leg_str,'location',loc);
        end
        
    case 'SW_up_TOA'
        plot(1847.5,me_t_PI,'ko','markerfacecolor','k','markersize',10);
        set(gca,'ylim',ylims_main_sw);
        
        if iplot_legend_sw==1
            legend(leg_str,'location',loc);
        end
    case 'calipso_total_cloud_amount'
        plot(1847.5,me_t_PI,'ko','markerfacecolor','k','markersize',10);
        set(gca,'ylim',ylims_main_totcf);
        %if iplot_legend_cf==1
            legend(leg_str,'location',loc);
        %end
end

ylims=get(gca,'ylim'); 
%title(tit_str_clean,'position',[1995 ylims(2)]);
title(tit_str_clean,'position',[yr_title,ylims(2)]);

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


savename=[savedir_date titlenam];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

%saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%Calculate the bias between model and obs for each year.

model_interp = interp1(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,years_obs);
prc_bias_annual = 100*(model_interp./obs_annual_box_cci - 1);
[stats_regions{ibox}]=calc_stats(1e-10,obs_annual_box_cci,model_interp);
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


savename=[savedir_date 'inset ' titlenam];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
title('');
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%title(titlenam);
    
end


