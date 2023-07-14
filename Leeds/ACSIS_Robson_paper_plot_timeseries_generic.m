%% UKESM eval
% Generally this is run from :-
% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic.m
% Need to run 
% ACSIS_Robson_paper_load_data  - sets var_ukesm in there
% before this to load in the data. Also sets savedir_date.

%ACSIS_Robson_paper_choose_regional_box2 is run to set the ylims for each
%variable, plus some other stuff.

% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

try

    if ~exist('ioverride_vals')
        ioverride_vals=0;
    end
    
%% Choose region
box_region = box_region_DRIVER;

iscreen_land=1;

%These are used in ACSIS_Robson_paper_CALC_timeseries2.m (called from ACSIS_Robson_paper_box_means.m) for the model data
% and in SW_Robson_obs_land_mask_FUNC2 (called from ACSIS_Robson_paper_plot_timeseries_obs_FUNC.m) for the obs
land_ocean = 'land+ocean';
land_ocean = 'ocean only';
%land_ocean = 'land only'; 
land_ocean = 'ocean only, no sea-ice';
%land_ocean = 'ocean, sea-ice regions only'

tail_test=2; %1 or 2-tailed test for signficance
%iplot_obs01=0;
inc_CERES = 0;
%inc_obs_01 = 0; %whether to include the other set of obs
%inc_obs_02 = 1; %whether to include the main set of obs
iload_amip_now=0;
 

if ioverride_vals==0;
    %iadd_amip = 1;
    %iload_amip=1;
    iadd_nudged = 0;    
end


end_str='';

%Set these to NaN to start with since they get saved later, but are only
%created if using obs.
trend_dat_box_obs1 = {NaN}; %need to be cells, though
trend_dat_box_obs2 = {NaN};
trend_dat_box_ens = {NaN};
model_dat_trend_period = {NaN};


%%
% ---
ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.
% ---



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

no_PI=1;




UKESM_Nd_case = 'to ztop';
UKESM_Nd_case = 'to 3.2km';



icoarse_grain=0;

time_round='';
time_format_str='';
icontour_DRIVER=0;
isave_plot=0;
iplot_wind_arrows=0;

cont_col_str_DRIVER='k';

switch var_str
    case 'Total cloud fraction'
        var_str_short = 'totCF'; %For saving figure, etc.   
    otherwise
        var_str_short = var_str;
end

switch var_ukesm
    case 'SWTOA Calc' %Override here since var_str is the same for the normal SW and the calculated one.
        %var_str_short is the one used for the savename for the figure.
        var_str_short = 'SW TOA up calc'; %For saving figure, etc.   
        
    case 'SO2_low_anthropogenic_emissions'
        iscreen_land=0;
        inc_obs_02=0;
end

%% Do the box mean and trend calculations

i_amip=0;
i_nudged=0;


%Run script to choose the data for the selected season and then to average
%over the box region.
yr_start_trend_box2 = yr_start_trend_box; yr_end_trend_box2 = yr_end_trend_box;
ACSIS_Robson_paper_box_means
%Returns itrend_start and -trend_end too

% %For the shorter timeseries (calculated SW)
% [dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std,trend_dat_box,...
%     trend_dat_box_ens,trend_ens_dat_box, dat_annual_box_ukesm_ens] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
%     yr_start_trend_box,yr_end_trend_box,p_conf,ibox);
% 
% %For the full ukesm (actual SW) timeseries
% %[dat_annual_box_ukesm_actual,dat_annual_box_ukesm_actual_ens_std,trend_dat_box_actual,...
% %    trend_dat_box_ens_actual,trend_ens_dat_box_actual, dat_annual_box_ukesm_ens_actual] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm_save,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
% %    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);

%% Do the AMIP box mean, etc. calculations if requested

if iadd_amip==1
    i_amip=1;    
    
    if iload_amip==1
        iload_amip_now=1;
        %Use this script to load the different .mat files, etc.:-
        ACSIS_Robson_paper_load_data_generic_multi_model
        iload_amip_now=0;
    end
    
    
%     %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_all_rsut.mat';
%     load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_all_' var_ukesm '.mat'];    
%     
%     load_file_PI = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_PI_PI_control_SW_up_TOA.mat';
%     
%     dat_PI_amip = load(load_file_PI,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     dat_ukesm_amip = load(load_file,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     
%     dat_ukesm_amip_DJF=load(load_file,'dat_annual_ens_DJF','years_ukesm_1d','dat_annual_DJF','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     dat_ukesm_amip_DJF.dat_annual_ens = dat_ukesm_amip_DJF.dat_annual_ens_DJF;
%     dat_ukesm_amip_DJF.dat_annual = dat_ukesm_amip_DJF.dat_annual_DJF;
%     
%     dat_ukesm_amip_MAM=load(load_file,'dat_annual_ens_MAM','years_ukesm_1d','dat_annual_MAM','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     dat_ukesm_amip_MAM.dat_annual_ens = dat_ukesm_amip_MAM.dat_annual_ens_MAM;
%     dat_ukesm_amip_MAM.dat_annual = dat_ukesm_amip_MAM.dat_annual_MAM;
%     
%     dat_ukesm_amip_JJA=load(load_file,'dat_annual_ens_JJA','years_ukesm_1d','dat_annual_JJA','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     dat_ukesm_amip_JJA.dat_annual_ens = dat_ukesm_amip_JJA.dat_annual_ens_JJA;
%     dat_ukesm_amip_JJA.dat_annual = dat_ukesm_amip_JJA.dat_annual_JJA;
%     
%     dat_ukesm_amip_SON=load(load_file,'dat_annual_ens_SON','years_ukesm_1d','dat_annual_SON','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
%     dat_ukesm_amip_SON.dat_annual_ens = dat_ukesm_amip_SON.dat_annual_ens_SON;
%     dat_ukesm_amip_SON.dat_annual = dat_ukesm_amip_SON.dat_annual_SON;
%     
%     %Run script to choose the data for the selected season and then to average
%     %over the box region.    
%     yr_start_trend_box2 = yr_start_trend_box_AMIP; yr_end_trend_box2 = yr_end_trend_box_AMIP;
%     end_str = '_amip';
%     ACSIS_Robson_paper_box_means
    
end

if iadd_nudged==1
    i_nudged=1;
    
    var_ukesm_nudged = '';
    switch var_ukesm
        case 'rsut'
            var_ukesm_nudged = 'SW_up_TOA';
        case 'clt'
            var_ukesm_nudged = 'tot_cloud_amount_in_rad';
    end
        
    
    %load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_AMIP_all_rsut.mat';
    %load_file = ['/home/disk/eos15/d.grosvenor/UM/ACSIS_Nd_trends/ensemble_timeseries_ukesm_nudged_CONTROL_all_' var_ukesm_nudged '.mat'];  %temp and wind nudging  
    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_Nudged_u-by844_all_' var_ukesm_nudged '.mat'];  %wind only nudging  
    load_file = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/EGSF_ensemble_timeseries_Nudged_u-bz785_all_' var_ukesm_nudged '.mat'];  %wind only nudging      
    load_file_PI = load_file;
    
    %load_file_PI = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_PI_PI_control_SW_up_TOA.mat';
    
    dat_PI_nudged = load(load_file_PI,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    dat_ukesm_nudged = load(load_file,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    
    dat_ukesm_nudged_DJF=load(load_file,'dat_annual_ens_DJF','years_ukesm_1d','dat_annual_DJF','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    dat_ukesm_nudged_DJF.dat_annual_ens = dat_ukesm_nudged_DJF.dat_annual_ens_DJF;
    dat_ukesm_nudged_DJF.dat_annual = dat_ukesm_nudged_DJF.dat_annual_DJF;
    
    dat_ukesm_nudged_MAM=load(load_file,'dat_annual_ens_MAM','years_ukesm_1d','dat_annual_MAM','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    dat_ukesm_nudged_MAM.dat_annual_ens = dat_ukesm_nudged_MAM.dat_annual_ens_MAM;
    dat_ukesm_nudged_MAM.dat_annual = dat_ukesm_nudged_MAM.dat_annual_MAM;
    
    dat_ukesm_nudged_JJA=load(load_file,'dat_annual_ens_JJA','years_ukesm_1d','dat_annual_JJA','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    dat_ukesm_nudged_JJA.dat_annual_ens = dat_ukesm_nudged_JJA.dat_annual_ens_JJA;
    dat_ukesm_nudged_JJA.dat_annual = dat_ukesm_nudged_JJA.dat_annual_JJA;
    
    dat_ukesm_nudged_SON=load(load_file,'dat_annual_ens_SON','years_ukesm_1d','dat_annual_SON','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    dat_ukesm_nudged_SON.dat_annual_ens = dat_ukesm_nudged_SON.dat_annual_ens_SON;
    dat_ukesm_nudged_SON.dat_annual = dat_ukesm_nudged_SON.dat_annual_SON;
    
    %Run script to choose the data for the selected season and then to average
    %over the box region.    
    yr_start_trend_box2 = yr_start_trend_box_nudged; yr_end_trend_box2 = yr_end_trend_box_nudged;
    end_str = '_nudged';
    ACSIS_Robson_paper_box_means
    
end


if iadd_DAMIP==1
    i_DAMIP=1;            
    
   % if iload_DAMIP==1
        %Use this script to load the different .mat files, etc.:-
        ACSIS_Robson_paper_load_data_generic
   % end
    
    inc_obs_01=0;
    inc_obs_02=0;
    
end

if iadd_AerChemMIP==1
    i_AerChemMIP=1;            
    
    if iload_AerChemMIP==1
        %Use this script to load the different .mat files, etc.:-
        ACSIS_Robson_paper_load_data_AerChemMIP
    end
    
    inc_obs_01=0;
    inc_obs_02=0;
    
end

if i_CMIP6_multimodel==1 | iadd_HADGEM==1           
    
    if iload_multi_model==1 | iload_HADGEM==1
        %Use this script to load the different .mat files, etc.:-
        ACSIS_Robson_paper_load_data_generic_multi_model
    end
    if iadd_HADGEM==0
        inc_obs_01=0;
        %inc_obs_02=0;
    end
    
    
    no_PI=1;
    i_CMIP6_multimodel2 = 1;
    
end

switch var_ukesm
    case 'SO2_low_anthropogenic_emissions'
        fscale=1e13; fscale_str=['\times10^{-13}'];
        
    case {'SW_up_TOA'}
        fscale=1e2; fscale_str=['\times10^{-2}'];
        
    case {'rsut'}
        fscale=1e2; fscale_str=['\times10^{-2}'];
        
    case {'calipso_total_cloud_amount','clt'}
        fscale=1e4; fscale_str=['\times10^{-4}'];
        
    case {'Nd_cf_weighted_UKESM','scldncl','Nd_clw_weighted_ESGF'}
        %fscale=1e4; fscale_str=['\times10^{-4}'];
        fscale=1; fscale_str=[''];            
        
    case {'dust_od550'}
        fscale=1e5; fscale_str=['\times10^{-5}'];
        
    case {'ts'}
        fscale=1e2; fscale_str=['\times10^{-2}'];      
        
end



%% Obs data - calc regional means, etc.

if inc_obs_02==1

    % Run FUNC here    
    obs_var_end_str='2'; %string to go on the end of the variables created, e.g. [obs_annual_box obs_var_str] which = obs_var_str2
    % Also does eval(['obs_str_in = obs_str' obs_var_end_str ';']); in the
    % script below - only requried for variables where can choose between
    % obs datasets.
    ACSIS_Robson_paper_plot_timeseries_obs_FUNC; %not really a function in fact        

%     if inc_CERES==1
%         % CERES
%         obs_str01='CERES';
%         years_obs = years_obs_CERES;
%         
%         %     if iscreen_land==1
%         %         [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(SW_TOA_ceres,gcm_Plat2D_CERES,gcm_Plon2D_CERES);
%         %     end
%         
%         clear opts; opts.screen_type_optional = land_ocean;
%         [dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC2(SW_TOA_ceres,gcm_Plat2D_CERES,gcm_Plon2D_CERES,opts);
%         
%         [obs_monthly_box,obs_annual_box,obs_annual_map] =...
%             ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES,season);
%     end


end


if inc_obs_01==1
    % Run FUNC here    
    obs_var_end_str='1'; %string to go on the end of the variables created, e.g. [obs_annual_box obs_var_str] which = obs_var_str2
    % Also does eval(['obs_str_in = obs_str' obs_var_end_str ';']); in the
    % script below - only requried for variables where can choose between
    % obs datasets.
    ACSIS_Robson_paper_plot_timeseries_obs_FUNC; %not really a function in fact        

end


%% Calculate the main trend for the boxed region
if inc_obs_02==1
    for it_trend=1:length(yr_start_trend_obs2)
        
        switch season
            case 'DJF'
                yr_start_trend_obs2 = yr_start_trend_obs2 + 1; %start a year later because data starts in Jan 1985
        end
        %yr_start_trend_used_box_obs2=yr_start_trend_obs2(it_trend); yr_end_trend_used_box_obs2=yr_end_trend_obs2(it_trend);
        
        %trend_dat_box_obs2 = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box2,years_obs2,yr_start,yr_end,p_conf,it_trend,ibox);
        [trend_dat_box_obs2{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_obs2, yr_end_trend_obs2, ...
            it_trend, years_obs2, obs_annual_box2 , p_conf);
    end
end

%% Calculate the 2nds OBS trend for the boxed region
if inc_obs_01==1
%     for it_trend=1:1 %length(yr_start_trend_box_obs)
%         yr_start=years_obs(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
%         switch season
%             case 'DJF'
%                 yr_start = yr_start + 1; %start a year later for now - although could do DJF of Y2000 for CERES
%         end
%         yr_start_trend_used_box_obs=yr_start; yr_end_trend_used_box_obs=yr_end;
%         
%         %trend_dat_box_obs = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box,years_obs,yr_start,yr_end,p_conf,it_trend,ibox);
%         [trend_dat_box_obs{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
%             it_trend, years_obs, obs_annual_box , p_conf);
%     end
    
     for it_trend=1:1 %length(yr_start_trend_box_obs)
        
        switch season
            case 'DJF'
                yr_start_trend_obs1 = yr_start_trend_obs1 + 1; %start a year later because data starts in Jan 1985
        end
        yr_start_trend_used_box_obs1=yr_start_trend_obs1; yr_end_trend_used_box_obs1=yr_end_trend_obs1;
        
        %trend_dat_box_obs1 = ACSIS_SW_paper_obs_trend_FUNC(obs_annual_box1,years_obs1,yr_start,yr_end,p_conf,it_trend,ibox);
        [trend_dat_box_obs1{ibox,it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_obs1, yr_end_trend_obs1, ...
            it_trend, years_obs1, obs_annual_box1 , p_conf);
    end
        
end

%% plot annual mean timeseries of calculated model SW vs obs
% Could put this UKESM code in a seprate script to use in modular form
%Full record is 1850 to 2014. We want 1984 to end.
%istart=1984-1850+1;
istart=yr_start_plot-dat_ukesm.years_ukesm_1d(1)+1;

titlenam = [season ' mean ' var_str_short  ' for region ' box_region];
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

if iadd_DAMIP==1 | iadd_AerChemMIP
        %istart_nudged=yr_start_plot-dat_ukesm_nudged.years_ukesm_1d(1)+1;        
    istart_DAMIP=1;
    
    %istart_trend_DAMIP=yr_start_trend_box - dat_ukesm_DAMIP.years_ukesm_1d(1)+1;
    
    %baseline = me_t_PI;
    
    baseline_type = 'individual';
    %baseline_type = 'zero';
    switch baseline_type
        case 'common'
            baseline = mean([dat_annual_box_ukesm_hist_aer(istart_DAMIP:istart_DAMIP+9) dat_annual_box_ukesm_hist_GHG(istart_DAMIP:istart_DAMIP+9)]);
        case 'individual'
            baseline = mean([dat_annual_box_ukesm(istart_DAMIP:istart_DAMIP+9)]); %use the HADGEM one for the HAdGEM line
        case 'zero'
            baseline = 0;
    end
    %baseline = 0;
    ylab_pre_str='\Delta';
else
    baseline = 0;
    ylab_pre_str='';
end

y_ens_upper = dat_annual_box_ukesm(istart:end) - baseline + dat_annual_box_ukesm_ens_std(istart:end)*2;
y_ens_lower = dat_annual_box_ukesm(istart:end) - baseline - dat_annual_box_ukesm_ens_std(istart:end)*2;
x_ens = dat_ukesm.years_ukesm_1d(istart:end)';
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
if iadd_DAMIP==0
    %p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
    p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none');
    %This line stops the line being listed in the legend :-
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

        
%h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
%p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
%Did the above to avoid weird legend behaviour - didn't quite work
hold on
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'bo-');


istart_trend=yr_start_trend_box-dat_ukesm.years_ukesm_1d(1)+1;
if iplot_individual_ens==1
   
    ens_str = ['iens=' num2str(iens_plot)];
    %h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_ens_actual(istart:end,iens_plot),'b-');
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm_ens(istart:end,iens_plot),'b-');
    %trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens_actual{ibox,1,iens_plot}.coeffs(2),'%.3f')];
    %trend_str = [',iens=' num2str(iens_plot) ', trend=' num2str(trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),'%.3f')];
    
    if iadd_trend_str==1
        trend_str = [ens_str ', trend=(' num2str(fscale*trend_dat_box_ens{ibox,1,iens_plot}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_ens{ibox,1,iens_plot}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str = [ens_str];
    end
    
    model_dat_trend_period{1} = dat_annual_box_ukesm_ens(itrend_start:itrend_end,iens_plot)';
else
    clear corr_timser
    iens_plot=1;
    
    switch var_ukesm
        case 'SWTOA Calc'
            ens_str = '';
        otherwise
            %ens_str = 'ensemble mean';
            ens_str='';
    end
    if iadd_DAMIP==1
        ens_str='';
    end
    
    %h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_actual(istart:end),'b-');
    h1=plot(dat_ukesm.years_ukesm_1d(istart:end),dat_annual_box_ukesm(istart:end)- baseline,'b-');    
    %trend_str=['trend=' num2str(trend_dat_box_actual{ibox,1}.coeffs(2),'%.3f')];
    %trend_str=['trend=' num2str(trend_dat_box{ibox,1}.coeffs(2),'%.3f')];
    if iadd_trend_str==1
        trend_str=[ens_str ', trend=(' num2str(fscale*trend_dat_box{ibox,itrend_lab}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box{ibox,itrend_lab}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str=[ens_str];
    end
    
    clear model_dat_trend_period
    for it_trend=1:length(itrend_start)
        model_dat_trend_period{it_trend} = dat_annual_box_ukesm(itrend_start(it_trend):itrend_end(it_trend))- baseline;
    end
end
%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}=[remove_character(model_str,'_',' ') trend_str]; ileg=ileg+1;
if ipad_legend==1
    %Pad the first entry to avoid box cutting into the text when have
    %superscripts
    leg_str{1}(2,:)=leg_str{1}(1,:); leg_str{1}(1,:)=' ';
end
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on

%% Plot AMIP
if iadd_amip==1
    %istart_amip=yr_start_plot-dat_ukesm_amip.years_ukesm_1d(1)+1;
    istart_amip=1;
    
    istart_trend_amip=yr_start_trend_box_AMIP - dat_ukesm_amip.years_ukesm_1d(1)+1;
    
    h_amip=plot(dat_ukesm_amip.years_ukesm_1d(istart_amip:end),dat_annual_box_ukesm_amip(istart_amip:end),'r-');
    if iadd_trend_str==1
        trend_str_amip=[', trend=(' num2str(fscale*trend_dat_box_amip{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_amip{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str_amip='';
    end
    
    %model_dat_trend_period = dat_annual_box_ukesm_amip(istart_trend_amip:end);
    
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['UKESM1-AMIP' trend_str_amip]; ileg=ileg+1;
    set(h_amip,'linewidth',4);
    set(h_amip,'markerfacecolor','r');
    
end

%% Plot nudged
if iadd_nudged==1
    %istart_nudged=yr_start_plot-dat_ukesm_nudged.years_ukesm_1d(1)+1;
    istart_nudged=1;
    
    istart_trend_nudged=yr_start_trend_box_nudged - dat_ukesm_nudged.years_ukesm_1d(1)+1;
    
    h_nudged=plot(dat_ukesm_nudged.years_ukesm_1d(istart_nudged:end),dat_annual_box_ukesm_nudged(istart_nudged:end),'g-');
    if iadd_trend_str==1
        trend_str_nudged=[', trend=(' num2str(fscale*trend_dat_box_nudged{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_nudged{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str_nudged='';
    end
    
    %model_dat_trend_period = dat_annual_box_ukesm_nudged(istart_trend_nudged:end);
    
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['UKESM1-NUDGED' trend_str_nudged]; ileg=ileg+1;
    set(h_nudged,'linewidth',4);
    set(h_nudged,'markerfacecolor','g');
    
end

%% Plot DAMIP
if iadd_DAMIP==1
    
    if iadd_trend_str==1
        trend_str_DAMIP=[', trend=(' num2str(fscale*trend_dat_box_DAMIP{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_DAMIP{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str_DAMIP='';
    end
    
     %Sum of perturbations from baseline    
     switch baseline_type
         case 'common'
             baseline_hist_aer = baseline;
             baseline_hist_GHG = baseline;
             baseline_hist_nat = baseline;
         case 'individual'
             baseline_hist_aer = mean([dat_annual_box_ukesm_hist_aer(istart_DAMIP:istart_DAMIP+9)]); %use the hist-aer one for the hist-aer line, etc.             
             baseline_hist_GHG = mean([dat_annual_box_ukesm_hist_GHG(istart_DAMIP:istart_DAMIP+9)]); %use the hist-aer one for the hist-aer line, etc.
             baseline_hist_nat = mean([dat_annual_box_ukesm_hist_nat(istart_DAMIP:istart_DAMIP+9)]); %use the hist-aer one for the hist-aer line, etc.
         case 'zero'
             baseline_hist_aer = 0;
             baseline_hist_GHG = 0;
             baseline_hist_nat = 0;
     end
    
    perts_hist_aer = dat_annual_box_ukesm_hist_aer(istart_DAMIP:end) - baseline_hist_aer;
    perts_hist_GHG = dat_annual_box_ukesm_hist_GHG(istart_DAMIP:end) - baseline_hist_GHG;
    perts_hist_nat = dat_annual_box_ukesm_hist_nat(istart_DAMIP:end) - baseline_hist_nat;
    
    h_hist_aer=plot(dat_ukesm_hist_aer.years_ukesm_1d(istart_DAMIP:end),perts_hist_aer,'k-'); 
    col_str=[0.6 0.6 0.6];  
    set(h_hist_aer,'color',col_str);
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['DAMIP-hist-aer' trend_str_DAMIP]; ileg=ileg+1;
    set(h_hist_aer,'linewidth',4);
    set(h_hist_aer,'markerfacecolor',col_str);
    
    h_hist_GHG=plot(dat_ukesm_hist_GHG.years_ukesm_1d(istart_DAMIP:end),perts_hist_GHG,'g-');        
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['DAMIP-hist-GHG' trend_str_DAMIP]; ileg=ileg+1;
    set(h_hist_GHG,'linewidth',4);
    set(h_hist_GHG,'markerfacecolor','g');
    
    h_hist_nat=plot(dat_ukesm_hist_nat.years_ukesm_1d(istart_DAMIP:end),perts_hist_nat,'r-');        
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['DAMIP-hist-nat' trend_str_DAMIP]; ileg=ileg+1;
    set(h_hist_nat,'linewidth',4);
    set(h_hist_nat,'markerfacecolor','r');
    
    %Test using the same method as the DAMIP lines
%     h_hist_all=plot(dat_ukesm_hist_ALL.years_ukesm_1d(istart_DAMIP:end),dat_annual_box_ukesm_hist_ALL(istart_DAMIP:end),'r-');        
%     %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
%     leg_str{ileg}=['DAMIP-ALL' trend_str_DAMIP]; ileg=ileg+1;
%     set(h_hist_all,'linewidth',4);
%     set(h_hist_all,'markerfacecolor','r');


%Calculate the sum of the individual DAMIP emissions lines for the
%"linear-sum" line.  
    perts_tot = perts_hist_aer + perts_hist_GHG + perts_hist_nat;
    dat_annual_box_ukesm_hist_linear = baseline+perts_tot;
    h_hist_linear=plot(dat_ukesm_hist_nat.years_ukesm_1d(istart_DAMIP:end),dat_annual_box_ukesm_hist_linear - baseline,'m-');        
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['DAMIP-linear sum' trend_str_DAMIP]; ileg=ileg+1;
    set(h_hist_linear,'linewidth',4);
    set(h_hist_linear,'markerfacecolor','m');
        
%Calculate  trends for linear sum line.        
    %The data is already averaged for the box region, so don't need to do
    %that again. Just calculate the trend using dat_annual_box_ukesm_hist_linear
    
    for it_trend=1:length(yr_start_trend_box)
        [trend_dat_box_hist_linear{ibox,it_trend},itrend_start(it_trend),itrend_end(it_trend)] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm_hist_linear, p_conf);
    end
    
% --- Proxy aerosol contribution ---    
%Calculate the difference between the HADGEM line and the DAMIP-HistGHG
%line as another estimate of the aerosol effect. I.e., are assuming that
%the HADGEM3 line is essentially GHGs+aerosols and so subtracting the GHG
%line would give the aerosol effect.

%dat_annual_box_ukesm(istart_DAMIP:end)
    perts_aerosol_proxy =  dat_annual_box_ukesm(istart:end)- baseline - perts_hist_GHG ;
    dat_annual_box_ukesm_hist_aer2 = baseline+perts_aerosol_proxy;
    col_str='k';
    %col_str=[0.6 0.6 0.6];
    h_hist_aer2=plot(dat_ukesm_hist_nat.years_ukesm_1d(istart_DAMIP:end),dat_annual_box_ukesm_hist_aer2 - baseline,[col_str '-']);        
    set(h_hist_aer2,'color',col_str);
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}=['DAMIP aerosol proxy' trend_str_DAMIP]; ileg=ileg+1;
    set(h_hist_aer2,'linewidth',4);
    set(h_hist_aer2,'markerfacecolor',col_str);
        
%Calculate  trends for the aerosol proxy line
    %The data is already averaged for the box region, so don't need to do
    %that again. Just calculate the trend using dat_annual_box_ukesm_hist_linear
    
    for it_trend=1:length(yr_start_trend_box)
        [trend_dat_box_hist_aer2{ibox,it_trend},itrend_start(it_trend),itrend_end(it_trend)] = ...
            ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm_hist_aer2, p_conf);
    end  

    
    %model_dat_trend_period = dat_annual_box_ukesm_DAMIP(istart_trend_DAMIP:end);
    

    
end

%% Plot AerChemMIP
if iadd_AerChemMIP==1
    
    istart_AerChemMIP_piaer = 1;
    
    if iadd_trend_str==1
        trend_str_AerChemMIP_piaer=[', trend=(' num2str(fscale*trend_dat_box_AerChemMIP_piaer{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_AerChemMIP_piaer{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str_AerChemMIP_piaer='';
    end
    
    %baseline_type = 'zero'; %set earlier
    %baseline_type = 'individual';
     %Sum of perturbations from baseline    
     switch baseline_type
         case 'common'
             baseline_AerChemMIP_piaer = baseline;
         case 'individual'
             baseline_AerChemMIP_piaer = mean([dat_annual_box_ukesm_AerChemMIP_piaer(istart_AerChemMIP_piaer:istart_AerChemMIP_piaer+9)]); %use the hist-aer one for the hist-aer line, etc.                          
         case 'zero'
             baseline_AerChemMIP_piaer = 0;             
     end
    
    perts_AerChemMIP_piaer = dat_annual_box_ukesm_AerChemMIP_piaer(istart_AerChemMIP_piaer:end) - baseline_AerChemMIP_piaer;
    
    h_AerChemMIP_piaer = plot(dat_ukesm_AerChemMIP_piaer.years_ukesm_1d(istart_AerChemMIP_piaer:end),perts_AerChemMIP_piaer,'g-');        
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    %leg_str{ileg}=['AerChemMIP-hist-PIaer' trend_str_AerChemMIP_piaer]; ileg=ileg+1;
    leg_str{ileg}=['AerChemMIP-GHG-only-proxy' trend_str_AerChemMIP_piaer]; ileg=ileg+1;
    set(h_AerChemMIP_piaer,'linewidth',4);
    set(h_AerChemMIP_piaer,'markerfacecolor','g'); 
    
%Calculate  trends for the GHG proxy line (PI aerosol)
    %The data is already averaged for the box region, so don't need to do
    %that again. Just calculate the trend using dat_annual_box_ukesm_hist_linear
    
    for it_trend=1:length(yr_start_trend_box)
        [trend_dat_box_AerChemMIP_piaer{ibox,it_trend},itrend_start(it_trend),itrend_end(it_trend)] = ...
            ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm_AerChemMIP_piaer, p_conf);
    end     
    
    
    % --- Proxy aerosol contribution ---    
%Calculate the difference between the UKESM line and the DAMIP-piAer
%line as another estimate of the aerosol effect. I.e., are assuming that
%the UKESM line is essentially GHGs+aerosols and so pi-Aer is just GHGs. So,
%subtracting the GHG line would give the aerosol effect.

%dat_annual_box_ukesm(istart_DAMIP:end)
    perts_aerosol_proxy =  dat_annual_box_ukesm(istart:end)- baseline - perts_AerChemMIP_piaer;
    dat_annual_box_ukesm_AerChemMIP_aero = baseline+perts_aerosol_proxy;
    %col_str='c';
    col_str=[0.6 0.6 0.6];
    h_AerChemMIP_aero=plot(dat_ukesm_AerChemMIP_piaer.years_ukesm_1d(istart_AerChemMIP_piaer:end),dat_annual_box_ukesm_AerChemMIP_aero...
        - baseline,'color',col_str,'linestyle','-');            
    %leg_str{ileg}=['AerChemMIP aerosol']; ileg=ileg+1;
    leg_str{ileg}=['AerChemMIP-aerosol-only-proxy']; ileg=ileg+1;    
    set(h_AerChemMIP_aero,'linewidth',4);
    set(h_AerChemMIP_aero,'markerfacecolor',col_str);
        
%Calculate  trends for the aerosol proxy line
    %The data is already averaged for the box region, so don't need to do
    %that again. Just calculate the trend using dat_annual_box_ukesm_hist_linear
    
    for it_trend=1:length(yr_start_trend_box)
        [trend_dat_box_AerChemMIP_aero{ibox,it_trend},itrend_start(it_trend),itrend_end(it_trend)] = ...
            ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend_box, yr_end_trend_box, ...
            it_trend, dat_ukesm.years_ukesm_1d, dat_annual_box_ukesm_AerChemMIP_aero, p_conf);
    end  
    
end


%% Plot CMIP6 models
if i_CMIP6_multimodel==1 | iadd_HADGEM==1
    %To check individual CMIP6 lines look at :-
    %[dat_annual_box_ukesm_' expt_str2 ]
    %E.g. dat_annual_box_ukesm_CESM2, etc. Years/time are in
    % dat_ukesm_CESM2.years_ukesm_1d
    ACSIS_Robson_paper_plot_multi_CMIP6_model
end

%% Plot obs
if inc_obs_01==1
    
    %col_01 = [0.35 0.625 0.1];
    col_01 = [0.5 0.5 0.5];
        
        %
        if length(years_obs)>1
            h2=plot(years_obs1,obs_annual_box1,'-','color',col_01);
        else
            h2=plot(years_obs1,obs_annual_box1,'s','markerfacecolor',col_01);
        end
        if iadd_trend_str==1
            trend_str=[', trend=(' num2str(fscale*trend_dat_box_obs1{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs1{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
        else
            trend_str='';
        end
        leg_str{ileg}=[obs_str1 trend_str]; ileg=ileg+1;
        set(h2,'linewidth',4);
        set(h2,'markerfacecolor',col_01);
        
end
    
if inc_obs_02==1       
    
    %Deep-C
    if length(years_obs2)>1
        h2=plot(years_obs2,obs_annual_box2,'k-');
    else
        h2=plot(years_obs2,obs_annual_box2,'ks');
    end
    if iadd_trend_str==1
        trend_str=[', trend=(' num2str(fscale*trend_dat_box_obs2{ibox,1}.coeffs(2),'%.2f') '\pm' num2str(fscale*trend_dat_box_obs2{ibox,1}.uncer_max,'%.2f') ') ' fscale_str];
    else
        trend_str='';
    end
    leg_str{ileg}=[obs_str2 trend_str]; ileg=ileg+1;        
    
    set(h2,'linewidth',4);
    set(h2,'markerfacecolor','k')            
    
    %Correlation coefficient over time of the model vs obs for the trend
    %period.
    for it_trend=1:length(yr_start_trend_obs2)
        yr_start_trend_used_box_obs2=yr_start_trend_obs2(it_trend); yr_end_trend_used_box_obs2=yr_end_trend_obs2(it_trend);        
        istart_trend_obs=find(years_obs2==yr_start_trend_used_box_obs2);
        iend_trend_obs=find(years_obs2==yr_end_trend_used_box_obs2);
        switch season
            case 'DJF'
                obs_dat = obs_annual_box2(istart_trend_obs+1:iend_trend_obs);
            otherwise
                obs_dat = obs_annual_box2(istart_trend_obs:iend_trend_obs);
        end
        
        
        
        if length(model_dat_trend_period{it_trend}) == length(obs_dat)
            corr_timser{it_trend}(iens_plot) = corr(model_dat_trend_period{it_trend}(:),obs_dat(:));
        else
            corr_timser{it_trend}=99;
        end
        
    end
    
    if length(dat_annual_box_ukesm(:)) == length(obs_annual_box2(:))
        corr_full_timser = corr(dat_annual_box_ukesm(:),obs_annual_box2(:));
    else
       corr_full_timser = 99;
    end
    
end


%These are set in ACSIS_Robson_paper_choose_regional_box2.m now
%loc='SouthWest';
%loc='NorthEast';
%loc='NorthWest';

%set(gca,'xlim',[1980 2014]);
%set(gca,'xlim',[1844 2018]);
%set(gca,'xlim',[yr_start_plot-2 yr_end_plot]);
set(gca,'xlim',plot_year_lims); %plot_year_lims set in RUN_multi

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
            %plot( x_trend , y_trend , 'b-','linewidth',3,'HandleVisibility','Off');
            hlin = plot( x_trend , y_trend , 'b-','linewidth',3);
            set(get(get(hlin,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            
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
    
    if iadd_HADGEM==1               
        for it_trend=1:length(eval(['trend_dat_box_' had_str]))
            if iplot_trend(it_trend)==1
                x_trend = eval(['trend_dat_box_' had_str '{ibox,it_trend}.x;']);
                y_trend = eval(['trend_dat_box_' had_str '{ibox,it_trend}.ylin;']);
                hlin_HadGEM3_GC31_LL = plot( x_trend , y_trend , 'color',cols{1},'linewidth',3);
                set(get(get(hlin_HadGEM3_GC31_LL,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end
    
    if iadd_amip==1
        for it_trend = itrend_amip
            x_trend = trend_dat_box_amip{ibox,it_trend}.x;
            y_trend = trend_dat_box_amip{ibox,it_trend}.ylin;
            hlin_amip = plot( x_trend , y_trend , 'r-','linewidth',3);
            set(get(get(hlin_amip,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    if iadd_nudged==1
        for it_trend=1:length(yr_start_trend_box_nudged)
            x_trend = trend_dat_box_nudged{ibox,it_trend}.x;
            y_trend = trend_dat_box_nudged{ibox,it_trend}.ylin;
            hlin_nudged = plot( x_trend , y_trend , 'g-','linewidth',3);
            set(get(get(hlin_nudged,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
    
    
    if inc_obs_01==1
        
        %plot CERES trend
        it_trend=1; %just one trend for obs since timeseries is short
        %plot( trend_dat_box_obs{ibox,it_trend}.x , trend_dat_box_obs{ibox,it_trend}.ylin , 'r-','linewidth',3,'HandleVisibility','Off');
        hobs = plot( trend_dat_box_obs1{ibox,it_trend}.x , trend_dat_box_obs1{ibox,it_trend}.ylin , 'm','color',col_01,'linewidth',3);
        %This line stops the line being listed in the legend :-
        set(get(get(hobs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        switch tail_test
            case 1
                conf_str = num2str(100*(1 - trend_dat_box_obs1{ibox,it_trend}.T1_max),'%1.1f');
            case 2
                conf_str = num2str(100*(1 - trend_dat_box_obs1{ibox,it_trend}.T2_max),'%1.1f');
        end
        conf_str_tit = [conf_str_tit '; p_{' obs_str1 '}=' conf_str '%'];
        
    end
    
    if inc_obs_02==1
        
        for it_trend=1:length(trend_dat_box_obs2)   
            %if iplot_trend(it_trend)==1
                %plot( trend_dat_box_obs2{ibox,it_trend}.x , trend_dat_box_obs2{ibox,it_trend}.ylin , 'k-','linewidth',3,'HandleVisibility','Off');
                hobs2 = plot( trend_dat_box_obs2{ibox,it_trend}.x , trend_dat_box_obs2{ibox,it_trend}.ylin , 'k-','linewidth',3);
                %This line stops the line being listed in the legend :-
                set(get(get(hobs2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                switch tail_test
                    case 1
                        conf_str = num2str(100*(1 - trend_dat_box_obs2{ibox,it_trend}.T1_max),'%1.1f');
                    case 2
                        conf_str = num2str(100*(1 - trend_dat_box_obs2{ibox,it_trend}.T2_max),'%1.1f');
                end
                conf_str_tit = [conf_str_tit '; p_{' obs_str2 '}=' conf_str '%'];
                tit_str_clean = [tit_str_clean conf_str_tit];
            %end
        end
        
    end
    
    
end


if no_PI==0 %set in ACSIS_Robson_paper_load_data.m
    % -- plot the error bar for the PI --
    %errorbar(1980,me_t_PI,std_t_PI*2,'CapSize',2,'ro','markerfacecolor','r');
    herr=errorbarYY('vert',yr_PI_bar,me_t_PI,std_t_PI*2,gca,'b','o',2,0.01);
    %This line stops the line being listed in the legend :-
    %set(get(get(herr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        %above didn't work for hrr
    set(herr,'HandleVisibility','Off'); 
    
    %errorbarYY('vert',1960,me_t_PI-std_t_PI,std_t_PI,gca,'r','o',2,0.01);
    
    hpi_mean = plot(yr_PI_bar,me_t_PI,'ko','markerfacecolor','b','markersize',10);
    %set(hpi_mean,'HandleVisibility','Off');
    %This line stops the line being listed in the legend :-
    set(get(get(hpi_mean,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

ylims_timser=get(gca,'ylim');        
switch var_ukesm
    case 'SO2_low_anthropogenic_emissions'
       ylims_timser =ylims_main_SO2;    
        %if iplot_legend==1
        %    legend(leg_str,'location',loc);
        %end
        iplot_legend = 0;
        
    case {'Nd_cf_weighted_UKESM','scldncl','Nd_clw_weighted_ESGF','Nd_clw_weighted_ESGF_no_dz'}
        ylims_timser = ylims_main_Nd;
        %title(titlenam,'position',[1925 200]);
        %title(tit_str_clean,'position',[1945 200]);
        %if iplot_legend==1
        %    legend(leg_str,'location',loc);
        %end
        
    case 'reffclwtop'
        ylims_timser = ylims_main_reff;
        
    case 'calipso_low_cloud_amount'
        
        
        ylims_timser =ylims_main_cf;
        %if iplot_legend_cf==1
        %    legend(leg_str,'location',loc);
        %end
        
    case {'SW_up_TOA','rsut','SWTOA Calc'}
        
        ylims_timser =ylims_main_sw;
        
%         if iplot_legend_sw==1
%             legend(leg_str,'location',loc);
%         end
        
    case {'SWTOA Calc'}        
        ylims_timser = ylims_main_swtoa_calc;
        loc='NorthEast';

    case {'calipso_total_cloud_amount','clt'}        
        ylims_timser =ylims_main_totcf;
        %if iplot_legend_cf==1
        %legend(leg_str,'location',loc);
        %end
        
    case {'ts'} %surface temperature
        
        ylims_timser = ylims_main_ts;
        
%         if iplot_legend_sw==1
%             legend(leg_str,'location',loc);
%         end    
        
     case {'dust_od550'} %surface temperature
        
        %set(gca,'ylim',ylims_main);
        
        %if iplot_legend_sw==1
            %legend(leg_str,'location',loc);
        %end    
        
    case {'lwp'}
        ylims_timser =ylims_main_lwp;
        
    case {'lwpic'}
        ylims_timser = ylims_main_lwpic;
    case {'rsutcs'}        
        ylims_timser = ylims_main_rsutcs;
    case {'od550tot'}
       ylims_timser =  ylims_main_od550tot;
            
end

% if iadd_DAMIP==1
%     ylims_timser = ylims_timser - baseline;
% end
set(gca,'ylim',ylims_timser);



if iplot_legend==1
    %h_leg = legend(leg_str,'location',loc,'fontsize',16);
    max_per_col = 2;
    ncols_leg = ceil(length(leg_str)/max_per_col);
    if ncols_leg>1
        h_leg = columnlegend(ncols_leg,leg_str,'location',loc,'fontsize',16);
    else
        h_leg = legend(leg_str,'location',loc,'fontsize',16);
    end
end




ylims=get(gca,'ylim'); 
%title(tit_str_clean,'position',[1995 ylims(2)]);
if inc_obs_02==1
    tit_str_clean = [tit_str_clean ', r_{corr}=' num2str(corr_full_timser(iens_plot),'%.2f')];
end

if ino_tit==1
    tit_str_clean='';
end

tit2 = {['Region=' box_region ' ' land_ocean],'', tit_str_clean};
title(tit_str_clean,'position',[yr_title,ylims(2)]);
title(tit2,'position',[yr_title,ylims(2)]);


switch var_str
    case 'Aerosol Optical Depth at 550nm'
        ylab_str = '\tau_{a}';
    otherwise
        ylab_str = var_str;
end
xlabel('Year');
if length(units_str)==0
    ylabel([ylab_pre_str ylab_str]);
else
    ylabel([ylab_pre_str ylab_str ' (' units_str ')']);
end

%switch box_region_DRIVER
%    case '1'
%        title('');
%        title(titlenam,'position',[1925 200]);
%    otherwise
%        title(titlenam,'position',[1925 200]);
%end


grid on

%return
%error('Ending early to avoid saves, etc.'); %end early to save time on saving, etc.
%Or could use return to pass back to invoking script - but want to clear the overrides, etc.

switch obs_str
    case {'none',''}
        save_obs_str='_no_obs';
    otherwise
        save_obs_str='';
end

%% Save main fig. and the trend, etc. data
if iadd_DAMIP==1
    damip_str='DAMIP';
    model_str_save='';
else
    damip_str='';
    model_str_save=model_str;
end
savename=[savedir_date titlenam ' ' ens_str ' ' land_ocean '_' model_str_save '_' damip_str period_str save_obs_str];
%savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
opts.iplot_png = 0;
opts.iplot_jpg = 0;
opts.isavefig = 1;
savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)


%Do another plot without xlabels and with the top ytick and title removed to allow
%them to stack together better in Latex
fig_file = [savename_out '.fig']
open(fig_file);
ytick_labs = get(gca,'yticklabel');
%ytick_labs(end,:)=' ';
%set(gca,'yticklabel',ytick_labs);
title('');
if iremove_xlabs==1
    xlabel(''); %remove xlabel
    set(gca,'xticklabel','');
end
if iplot_legend==1    
     if length(leg_str)>3
        h_leg = columnlegend(ncols_leg,leg_str,'location',loc,'fontsize',16);
    else
        h_leg = legend(leg_str,'location',loc,'fontsize',16);
    end
end

savename_out_nolabs = saveas_ps_fig_emf(gcf,[savename '_nolabs'],'',0,0,0,'',[],0,opts)

%Save trends, etc. - when doing this for the no-obs table (SW trend calcs,
%etc.) set obs_str='none' before running the RUN_multi script
save_trend_dat_name = [savename_out '.mat'];
save_trend_dat_name2 = remove_character(save_trend_dat_name,' ','_');
save_trend_dat_name2 = remove_character(save_trend_dat_name2,'%','_');
save_trend_dat_name2 = remove_character(save_trend_dat_name2,':','_');


years_ukesm_1d = dat_ukesm.years_ukesm_1d;
save(save_trend_dat_name2,'trend_dat_box','trend_dat_box_obs1','trend_dat_box_obs2','trend_dat_box_ens','yr_start_trend_box2'...
    ,'yr_end_trend_box2','years_ukesm_1d','dat_annual_box_ukesm','dat_annual_box_ukesm_ens'...
    ,'-V7.3');
if iadd_DAMIP==1
   save(save_trend_dat_name2,'trend_dat_box_hist_GHG','trend_dat_box_hist_aer','trend_dat_box_hist_nat','trend_dat_box_hist_linear'...
       ,'trend_dat_box_hist_aer2'...
       ,'dat_annual_box_ukesm_hist_aer','dat_annual_box_ukesm_hist_GHG','dat_annual_box_ukesm_hist_nat','dat_annual_box_ukesm_hist_linear' ...
       ,'dat_annual_box_ukesm_hist_aer2'...
    ,'-V7.3','-APPEND'); 
end
if iadd_AerChemMIP==1
   save(save_trend_dat_name2,'trend_dat_box_AerChemMIP_aero','trend_dat_box_AerChemMIP_piaer'...
       ,'dat_annual_box_ukesm_AerChemMIP_aero','dat_annual_box_ukesm_AerChemMIP_piaer'...
       ,'-V7.3','-APPEND'); 
end

%Calculate the bias between model and obs for each year.


%model_interp = interp1(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,years_obs2);
if (inc_obs_01==1 | inc_obs_02==1) & (length(model_dat_trend_period) == length(obs_dat))
    prc_bias_annual = 100*(model_dat_trend_period./ obs_dat- 1);
    [stats_regions{ibox}]=calc_stats(1e-10,obs_dat,model_dat_trend_period);
    stats_regions{ibox}.box_region_str = box_region_str;
end


%% plot INSET annual mean timeseries of model vs obs for 2003 to 2014
iplot_inset=0;

if iplot_inset==1
    titlenam = ['Annual mean ' var_str  ' for region ' box_region];
    tit_str_clean = ['Region ' box_region_str];
    %figure
    
    %ax2 = axes('Position',[0.17 0.7 0.25 0.25]);
    %ax2 = axes('Position',[0.1748    0.7683    0.2500    0.2500]);
    
    switch var_ukesm
        case {'Nd_cf_weighted_UKESM','scldncl'}
            ax2 = axes('Position',inset_ax_pos); %inset_ax_pos set in ACSIS_Robson_paper_choose_regional_box2.m
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
    %patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
    p = patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
    %This line stops the line being listed in the legend :-
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
    hold on
    h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,'bo-');
    %h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'b-');
    
    %leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
    leg_str{ileg}='UKESM'; ileg=ileg+1;
    set(h1,'linewidth',4);
    set(h1,'markerfacecolor','b');
    hold on
    
    
    %h2=plot(years_obs,obs_annual_box,'rs-');
    h2=plot(years_obs,obs_annual_box2,'ks-');
    
    %leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
    %leg_str{ileg}=['CALIPSO'];
    leg_str{ileg}=obs_str;
    set(h2,'linewidth',4);
    set(h2,'markerfacecolor','k');
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
        case {'SW_up_TOA','rsut'}
            set(gca,'ylim',ylims_inset_sw);
        case {'scldncl'}
            set(gca,'ylim',[100 145]);
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
    [savename_out] = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
    %title(titlenam);
    
end


%% Do box and whisker plot
iplot_box_whisker = 1;
if iplot_individual_ens~=1 & iplot_box_whisker==1
    ACSIS_Robson_paper_plot_ens_trend_Box_Whiskier_generic
    [savename_out_box_whisk] = saveas_ps_fig_emf(gcf,[savename '_box_whisk'],'',0,0,0,'',[],0,opts);
end

%clear the override and also clear if there is an error in the code
clear ioverride_vals

catch timser_error
    clear ioverride_vals
    rethrow(timser_error);
end



