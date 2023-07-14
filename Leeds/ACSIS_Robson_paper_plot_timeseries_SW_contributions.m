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
ACSIS_Robson_paper_choose_regional_box_SW_contributions %run script - also chooses ylims, etc.
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


%% Do the calculations

%For the shorter timeseries (calculated SW - all varying)
[dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std,trend_dat_box] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);

%CF constant
[dat_annual_box_ukesm_cf,dat_annual_box_ukesm_ens_std_cf,trend_dat_box_cf] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm_sw_sens,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);

%dat_ukesm_sw_sens seems to be set in ACSIS_Robson_paper_timeseries_etc_FUNC.m

%Nd constant
%[dat_annual_box_ukesm_Nd,dat_annual_box_ukesm_ens_std_Nd,trend_dat_box_Nd] = ACSIS_Robson_paper_CALC_timeseries(dat_ukesm_Nd,dat_PI,LAT_val,LON_val,iscreen_land,gcm_area_UM,...
%    yr_start_trend_box,yr_end_trend_box,p_conf,ibox);

%% Obs data - calc regional means, etc.

switch var_ukesm
    
     case 'SW_up_TOA';
         switch SW_dataset
             case 'CERES'                 
                 obs_str='CERES';
                 years_obs = [2001:2018];
                 
                 dat_tmp = SW_TOA_ceres;
                 if iscreen_land==1
                     land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
                     lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
                     lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean
                     lmask_CERES = griddata(land_mask_MODIS.gcm_Plat2D_AMSRE,land_mask_MODIS.gcm_Plon2D_AMSRE,lmask_MODIS,gcm_Plat2D_CERES,gcm_Plon2D_CERES,'nearest');
                     for it=1:size(SW_TOA_ceres,1)
                         dat_tmp(it,:,:) = squeeze(dat_tmp(it,:,:)).*lmask_CERES;
                     end
                     
                     
                 end
                 
                 [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES);
                 
             case 'DEEP-C';
                 obs_str='DEEP-C';
                 years_obs = [1985:2014];                 
                 dat_tmp = dat_deepc.sw_up_toa;
                 
                 if iscreen_land==1
                     land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
                     lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
                     lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean
                     lmask_deepc = griddata(land_mask_MODIS.gcm_Plat2D_AMSRE,land_mask_MODIS.gcm_Plon2D_AMSRE,lmask_MODIS,gcm_Plat2D_UM,gcm_Plon2D_UM,'nearest');
                     for it=1:size(dat_deepc.sw_up_toa,1)
                         dat_tmp(it,:,:) = squeeze(dat_tmp(it,:,:)).*lmask_deepc;
                     end                                          
                 end
                 
                 [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_val,LON_val,years_obs,gcm_area_UM);
         end
    
    case 'calipso_low_cloud_amount';
        obs_str='CALIPSO';
        
        
        %Load Calipso data - currently 2007-2017
        % Load CF data using :-
        %read_calipso_monthly_IPSL_2007_2017
        years_obs = years_requested;
        
        dat_tmp = 0.01*cllcalipso_monthly_AVERAGE;
        if iscreen_land==1
        	land_mask_MODIS=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
            lmask_MODIS = flipdim(land_mask_MODIS.amsre_land_mask,1);
            lmask_MODIS = lmask_MODIS + 1; %Make it ones where have ocean    
            lmask_CAL = griddata(land_mask_MODIS.gcm_Plat2D_AMSRE,land_mask_MODIS.gcm_Plon2D_AMSRE,lmask_MODIS,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,'nearest');
            for it=1:size(cllcalipso_monthly_AVERAGE,1)
                dat_tmp(it,:,:) = squeeze(dat_tmp(it,:,:)).*lmask_CAL;
            end

            
        end
        
        [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,LAT_val,LON_val,years_obs,gcm_area_CALIPSO_monthly);
        
%         ilat = find(gcm_Plat2D_CALIPSO_monthly(:,1)>LAT_val(1) & gcm_Plat2D_CALIPSO_monthly(:,1)<LAT_val(2));
%         ilon = find(gcm_Plon2D_CALIPSO_monthly(1,:)>LON_val(1) & gcm_Plon2D_CALIPSO_monthly(1,:)<LON_val(2));
%         %     clear dat_time_mean
%         %     for iy=1:size(Nd_ukesm,1)
%         %         dat_time_mean(iy,:) = meanNoNan(meanNoNan(dat.dat_ukesm(iy,:,ilat,ilon),4),2);
%         %     end
%         
%         
%         clear obs_annual_box
%         for it=1:size(cllcalipso_monthly_AVERAGE,1)
%             obs_monthly_box(it) = meanNoNan(meanNoNan(0.01*cllcalipso_monthly_AVERAGE(it,ilat,ilon),3),2);
%             %dat_tmp = meanNoNan(meanNoNan(dat.dat_annual_ens(:,it,ilat,ilon),4),2);
%             %dat_annual_box_ukesm_ens_std(it) = std(dat_tmp); %std dev across the ensemble
%         end
%         
%         istart=1; iend=12;
%         clear obs_annual_box obs_annual_map
%         for it=1:length(years_obs)
%             obs_annual_box(it) = meanNoNan(obs_monthly_box(istart:iend),2);
%             obs_annual_map(it,:,:) = meanNoNan(0.01*cllcalipso_monthly_AVERAGE(istart:iend,:,:),1);
%             istart=istart+12;
%             iend=iend+12;
%         end
        

        
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
    yr_start=years_obs(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
    yr_start_trend_used_box_obs=yr_start; yr_end_trend_used_box_obs=yr_end;    
    istart=find(years_obs==yr_start);
    iend=find(years_obs==yr_end);    
    x = [yr_start:yr_end]';
    y = obs_annual_box(istart:iend);

    %test whether the significance is affected by subtracting the mean -
    % -- makes no difference --
    %me_lin = meanNoNan(obs_annual_box(istart:iend),2);
    %y = obs_annual_box(istart:iend) - me_lin;
    
if length(y)>2
    
   [coeffs,t_trend] = trend_linear_fit_nd_data(x,y,2); 
   
   ylin=coeffs(1)+coeffs(2).*x; %The straight line for the linear trend
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    %t_thresh = tinv(p_conf/100,n_dof); %find the t value needed for 95% confidence using a one-tailed t-test
        %N.B. - here is the function for a 2-tailed test :-
            % E.g. 
            % t=4; v=10;
            % tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));   % 2-tailed t-distribution function
            % tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2; %1-tailed t-distribution function
            % OR just :-  
            % tail2P = 2*tcdf(-abs(t),v);
            % tail1P = tcdf(-abs(t),v);
            % Test :-            
            % T2 = [1-tdist2T(t,v)  tail2P]; %give the same answer
            % T1 = [1-tdist1T(t,v)  tail1P]; %give the same answer
    %itrend_not_sig=find(abs(t_trend)<t_thresh); 
    
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend),n_dof);
    T1 = T2/2;    
    itrend_not_sig = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig = find((1-T1)<=p_conf/100); %1-tailed

    trend_dat_box_obs{ibox,it_trend}.coeffs = coeffs;
    trend_dat_box_obs{ibox,it_trend}.t_trend = t_trend;
    trend_dat_box_obs{ibox,it_trend}.ylin = ylin;
    trend_dat_box_obs{ibox,it_trend}.x = x;
    trend_dat_box_obs{ibox,it_trend}.itrend_not_sig = itrend_not_sig;
    trend_dat_box_obs{ibox,it_trend}.T1 = T1;
    trend_dat_box_obs{ibox,it_trend}.T2 = T2;
    
end
    
end

%%
titlenam = ['Annual mean ' var_str  ' for region ' box_region];
tit_str_clean = [box_region_str];
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
increase_font_size_map_figures
set(gcf,'position',[3         297        1256         322]);
set(gca,'position',[0.1300    0.1500    0.7750    0.7150]);


%% Plot the linear trend lines

if iplot_lin==1
   conf_str_tit='';
   
   %Model calculated SW
   for it_trend=1:length(yr_start_trend_box)
       if iplot_trend(it_trend)==1
           plot( trend_dat_box{ibox,it_trend}.x , trend_dat_box{ibox,it_trend}.ylin , 'b-','linewidth',3,'HandleVisibility','Off');
           hold on
       end
      yr_start_trend_box_str = num2str(yr_start_trend_box(it_trend));
      yr_end_trend_box_str = num2str(yr_end_trend_box(it_trend));
      switch tail_test
          case 1
              conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T1),'%1.1f');     
          case 2      
              conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T2),'%1.1f');
      end            
      
      conf_str_tit = [conf_str_tit '; p_' num2str(it_trend) '=' conf_str '%'];
      %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
   end
   
   %Model calculated SW w/ constant CF
   for it_trend=1:length(yr_start_trend_box)
       if iplot_trend(it_trend)==1
           plot( trend_dat_box_cf{ibox,it_trend}.x , trend_dat_box_cf{ibox,it_trend}.ylin , 'g-','linewidth',3,'HandleVisibility','Off');
       end
      yr_start_trend_box_str = num2str(yr_start_trend_box(it_trend));
      yr_end_trend_box_str = num2str(yr_end_trend_box(it_trend));
      switch tail_test
          case 1
              conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T1),'%1.1f');     
          case 2      
              conf_str = num2str(100*(1 - trend_dat_box{ibox,it_trend}.T2),'%1.1f');
      end            
      
      conf_str_tit = [conf_str_tit '; p_' num2str(it_trend) '=' conf_str '%'];
      %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
   end
   
%    %Model actual SW
%    for it_trend=1:length(yr_start_trend_box)
%        if iplot_trend(it_trend)==1
%            plot( trend_dat_box_actual{ibox,it_trend}.x , trend_dat_box_actual{ibox,it_trend}.ylin , 'k-','linewidth',3,'HandleVisibility','Off');
%        end
%       yr_start_trend_box_str = num2str(yr_start_trend_box(it_trend));
%       yr_end_trend_box_str = num2str(yr_end_trend_box(it_trend));
%       switch tail_test
%           case 1
%               conf_str = num2str(100*(1 - trend_dat_box_actual{ibox,it_trend}.T1),'%1.1f');     
%           case 2      
%               conf_str = num2str(100*(1 - trend_dat_box_actual{ibox,it_trend}.T2),'%1.1f');
%       end            
%       
%       conf_str_tit = [conf_str_tit '; p_' num2str(it_trend) '=' conf_str '%'];
%       %leg_str{ileg}=['Trend ' yr_start_trend_box_str ' to ' yr_end_trend_box_str '; ' conf_str '%']; ileg=ileg+1;
%    end
%    
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


%% plot annual mean timeseries of calculated model SW vs actual


clear leg_str
ileg=1;

%Full record is 1850 to 2014. We want 1984 to end.
% istart=1984-1850+1;
% h3=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_actual(istart:end),'k-');
% leg_str{ileg}='Actual'; ileg=ileg+1;
% set(h3,'linewidth',4);
% set(h3,'markerfacecolor','k');
% hold on

%h1=errorbar(dat.years_ukesm_1d,dat_annual_box_ukesm,dat_annual_box_ukesm_ens_std*2,'bo-');
y_ens_upper = dat_annual_box_ukesm + dat_annual_box_ukesm_ens_std*2;
y_ens_lower = dat_annual_box_ukesm - dat_annual_box_ukesm_ens_std*2;
x_ens = dat_ukesm.years_ukesm_1d';
y_ens_patch = [y_ens_upper fliplr(y_ens_lower) y_ens_upper(1)];
x_ens_patch = [x_ens fliplr(x_ens) x_ens(1)];  
p=patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
hold on
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'bo-');
h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,'b-');

%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}=['Calculated x ' num2str(f_sw_calc) ', trend=' num2str(trend_dat_box{ibox,it_trend}.coeffs(2),'%.3f') ]; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on

%Full record is 1850 to 2014. We want 1984 to end.
istart=1984-1850+1;
h3=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm_cf,'g-');
leg_str{ileg}=[sw_str ', trend=' num2str(trend_dat_box_cf{ibox,it_trend}.coeffs(2),'%.3f')]; ileg=ileg+1;
set(h3,'linewidth',4);
set(h3,'markerfacecolor','g');




%h2=plot(years_obs,obs_annual_box,'rs-');
if length(years_obs)>1
    h2=plot(years_obs,obs_annual_box,'r-');
else
    h2=plot(years_obs,obs_annual_box,'rs');
end
%leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
%leg_str{ileg}=['CALIPSO'];
leg_str{ileg}=[obs_str ', trend=' num2str(trend_dat_box_obs{ibox,it_trend}.coeffs(2),'%.3f')]; ileg=ileg+1;
set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');
loc='SouthWest';
%loc='NorthEast';

set(gca,'xlim',[1980 2014]);
%set(gca,'xlim',[1844 2018]);



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
        
end

ylims=get(gca,'ylim'); 
title(tit_str_clean,'position',[1965 ylims(2)]);

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

model_interp = interp1(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,years_obs)
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


