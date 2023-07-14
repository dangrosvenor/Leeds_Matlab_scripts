%% UKESM eval
%Need to run 
%ACSIS_Robson_paper_load_data  - sets var_ukesm in there
%before this to load in the data. Also sets savedir_date.

% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

%% Choose region
% Box for regional averages
%box_region_DRIVER = '1'; %US small box
%box_region_DRIVER = '3'; %SW of Spain
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
tail_test=1; %1 or 2-tailed test for signficance

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


%% calculate some timeseries in the regional box for the model
[dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std,trend_dat_box,...
    trend_dat_box_ens,trend_ens_dat_box, dat_annual_box_ukesm_ens,me_t_PI,N_t_PI,std_t_PI] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);
    

%% Obs data - calc regional means, etc.

switch var_ukesm
    
     case 'SW_up_TOA';
        obs_str='CERES';        
        years_obs = [2001:2018];  
        years_obs = [2001:2015];
        
        if iscreen_land==1
            [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(SW_TOA_ceres,gcm_Plat2D_CERES,gcm_Plon2D_CERES);
        end
        
       % [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES);
        
        [obs_monthly_box,obs_annual_box,obs_annual_map] =...
    ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES,season);

                        
     case 'DEEPC_fluxes';
        obs_str='DEEPC';        
        years_obs = [1985:2015];            
        
        if iscreen_land==1
            [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_deepc.sw_up_toa,gcm_Plat2D_UM,gcm_Plon2D_UM);
        end
        
        [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_deepc,gcm_Plon2D_deepc,LAT_val,LON_val,years_obs,gcm_area_deepc);
        
    
    case {'calipso_low_cloud_amount','calipso_total_cloud_amount'};
        cf_obs = 'CALIPSO';
        cf_obs = 'ESA CCI';
        switch cf_obs
            case 'CALIPSO'
                obs_str='CALIPSO';
                                
                %Load Calipso data - currently 2007-2017
                % Load CF data using :-
                %read_calipso_monthly_IPSL_2007_2017
                years_obs = years_requested;
                
                dat_obs = 0.01*cllcalipso_monthly_AVERAGE;
                if iscreen_land==1
                    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly);
                end
                
                [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,LAT_val,LON_val,years_obs,gcm_area_CALIPSO_monthly);
                
            case 'ESA CCI'
                obs_str='ESA CCI';
                                
                years_obs = unique(cci_dat.Y_out);
                
                dat_obs = cci_dat.cfc; %total CF (day and night)
                if iscreen_land==1
                    [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_obs,gcm_Plat2D_CCI,gcm_Plon2D_CCI);
                end
                
                [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CCI,gcm_Plon2D_CCI,LAT_val,LON_val,years_obs,gcm_area_CCI);
                
                
                
                
        end
        
    case 'Nd_cf_weighted_UKESM';
        obs_str='MODIS';
        
        
        % MODIS data
        % Using the data given to Jane - screened for sea-ice etc.
        cf_screen_str = 'CF>80';
        cf_screen_str = 'CF>0';
        
        str_2137='21';
        str_2137='37';
        
        res_str='1deg';
        res_str='1km';
        
        
        switch cf_screen_str
            case 'CF>80'
                file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';
                dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km_screened_for_seaice__2week_max';
                dataset_str = 'SZA_LT_65_CF_GT_80_CTH_LT_3.2km';
                
            case 'CF>0'
                file_dir='/home/disk/eos1/d.grosvenor/mock_L3/CF_0.0_meanCTT_173_meanCTH_3.2km_SZA_65/';
                dataset_str = 'SZA_LT_65_CF_GT_0_CTH_LT_3.2km';
        end
        
        
        
        
        switch str_2137
            case '21'
                str_label_2137='2.1 um';
            case '37'
                str_label_2137='3.7 um';
        end
        
        years_MODIS2=[2003:2014];        
        clear mon_me_MODIS2 mon_me_MODIS2_Ndatap
        for iy=1:length(years_MODIS2)
            year_str = num2str(years_MODIS2(iy));
            filename = [file_dir 'Nd_monthly_' str_2137 '_' res_str '_' year_str '_' dataset_str '.mat.nc'];
            nc=netcdf(filename);
            if iy==1
                lat=nc{'lat'}(:);
                lon=nc{'lon'}(:);
                [gcm_Plon2D_AMSRE,gcm_Plat2D_AMSRE]=meshgrid(lon,lat);
            end
            
            mon_me_MODIS2{iy} = nc{['Nd_' res_str '_mean']}(:);
            mon_me_MODIS2_Ndatap{iy} = nc{['Nd_' res_str '_Ndatap']}(:);
            inan = find(mon_me_MODIS2_Ndatap{iy} < nthresh_days);
            mon_me_MODIS2{iy}(inan)=NaN;
            %mon_me_filter{iy}(inan)=NaN;
        end
        
        if iscreen_land==1
        	land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
            lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
            lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean                      
        end
        
        [mon_me_region_MODIS2] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val,iscreen_land,lmask_MODIS,1,gcm_area_AMSRE);
        [mon_me_region_MODIS2_Ndatap] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2_Ndatap,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val,iscreen_land,lmask_MODIS,1,gcm_area_AMSRE);
        %[mon_me_filter_SO] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_filter,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,thresh_LAT,thresh_LON);
        
        for iy=1:length(years_MODIS2)
            Nd_MODIS(iy,:) = mon_me_region_MODIS2{iy}(:); %order is [year month]
        end
        
        
        
        
        clear Nd_annual_box_MODIS Nd_annual_MODIS obs_annual_box
        for iy=1:size(Nd_MODIS,1)
            %Nd_annual_box_MODIS(iy) = meanNoNan(Nd_MODIS(iy,:),2);
            obs_annual_box(iy) = meanNoNan(Nd_MODIS(iy,:),2);
            Nd_annual_MODIS(iy,:,:) = meanNoNan(mon_me_MODIS2{iy},3);
        end
        
        years_obs = years_MODIS2;
        
end

%% Calculate the OBS trend for the boxed region
for it_trend=1:1 %length(yr_start_trend_box_obs)
    if ioverride_obs_trend_yr_start==1
        yr_start = yr_trend_obs_start;
    else
        yr_start=years_obs(1);
    end
    yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
    yr_start_trend_used_box_obs=yr_start; yr_end_trend_used_box_obs=yr_end;    
    
    %trend_dat_box_obs = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box,years_obs,yr_start,yr_end,p_conf,it_trend,ibox);
     [trend_dat_box_obs{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs, obs_annual_box , p_conf);
end

%% plot annual mean timeseries of model vs obs for chosen period
istart=max([1 yr_start_plot-1850+1]);

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
y_ens_upper = dat_annual_box_ukesm + dat_annual_box_ukesm_ens_std*2;
y_ens_lower = dat_annual_box_ukesm - dat_annual_box_ukesm_ens_std*2;
x_ens = dat_ukesm.years_ukesm_1d';
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
hold on

if iplot_individual_ens==1
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
    %trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),sigfig_str)];
    trend_str = [''];
else
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm(istart:end),'b-');
    %trend_str=[', trend=' num2str(trend_dat_box{ibox,1}.coeffs(2),sigfig_str)];
    trend_str=[''];
end

%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
%leg_str{ileg}='UKESM1'; ileg=ileg+1;
leg_str{ileg}=['UKESM1' trend_str]; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on


%h2=plot(years_obs,obs_annual_box,'rs-');
if length(years_obs)>1
    h2=plot(years_obs,obs_annual_box,'r-');
else
    h2=plot(years_obs,obs_annual_box,'rs');
end
%leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
%leg_str{ileg}=['CALIPSO'];
%leg_str{ileg}=obs_str; ileg=ileg+1;
trend_str = [', trend=' num2str(trend_dat_box_obs{ibox,1}.coeffs(2),sigfig_str)];
trend_str='';
leg_str{ileg}=[obs_str trend_str]; ileg=ileg+1;
set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');


%set(gca,'xlim',[1980 2018]);
set(gca,'xlim',[yr_start_plot 2018]);

%% Plot the linear trends
if iplot_lin==1
   conf_str_tit='';
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
      
      conf_str_tit = [conf_str_tit '; p_' num2str(it_trend) '=' conf_str '%'];
      %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
   end
   
   %plot obs trend
   it_trend=1; %just one trend for obs since timeseries is short
   plot( trend_dat_box_obs{ibox,it_trend}.x , trend_dat_box_obs{ibox,it_trend}.ylin , 'r-','linewidth',3,'HandleVisibility','Off');   
   switch tail_test
       case 1
           conf_str = num2str(100*(1 - trend_dat_box_obs{ibox,it_trend}.T1),'%1.1f');
       case 2
           conf_str = num2str(100*(1 - trend_dat_box_obs{ibox,it_trend}.T2),'%1.1f');
   end
   conf_str_tit = [conf_str_tit '; p_{obs}=' conf_str '%'];
   
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
            legend(leg_str,'location',loc_DRIVER);
        end
    case 'calipso_low_cloud_amount'
        plot(1847.5,me_t_PI,'ko','markerfacecolor','k','markersize',10);
        set(gca,'ylim',ylims_main_cf);        
        if iplot_legend_cf==1
            legend(leg_str,'location',loc_DRIVER);
        end
        
    case 'SW_up_TOA'
        plot(1847.5,me_t_PI,'ko','markerfacecolor','k','markersize',10);
        set(gca,'ylim',ylims_main_sw);
        
        if iplot_legend_sw==1
            legend(leg_str,'location',loc_DRIVER);
        end
        
    case 'calipso_total_cloud_amount'
        set(gca,'ylim',ylims_main_totcf);
        legend(leg_str,'location',loc_DRIVER);
        
    otherwise
        legend(leg_str,'location',loc_DRIVER);
        
end

ylims=get(gca,'ylim'); 
title(tit_str_clean,'position',[yr_title ylims(2)]);

xlabel('Year');
ylabel([var_str ' ' units_str]);

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
prc_bias_annual = 100*(model_interp./obs_annual_box - 1);
[stats_regions{ibox}]=calc_stats(1e-10,obs_annual_box,model_interp);
stats_regions{ibox}.box_region_str = box_region_str;

%% plot INSET annual mean timeseries of model vs obs for 2003 to 2014
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
loc='NorthEast';

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
%title(titlenam);
title('');
grid on

end

savename=[savedir_date 'inset ' titlenam];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%title(titlenam);
    



