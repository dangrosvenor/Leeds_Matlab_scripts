%% UKESM eval
savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
eval(['!mkdir ' savedir_date]); 

%% Choose region
% Box for regional averages
box_region_DRIVER = '1'; %US small box
box_region_DRIVER = '3'; %SW of Spain
%box_region_DRIVER = '4'; %All NA region
%box_region_DRIVER = '5'; %UK region
%box_region_DRIVER = '6'; %West of UK region
%box_region_DRIVER = '7'; %20 deg further west of there
%box_region_DRIVER = '8'; %middle of Atlantic at US box latitude
box_region_DRIVER = '9'; %
box_region_DRIVER = '10'; %

box_region = box_region_DRIVER;
ACSIS_Robson_paper_choose_regional_box %run script

iscreen_land=0;

SW_dataset = 'CERES';
SW_dataset = 'DEEP-C';

%%



tit_str_clean='';
i_plot_all_boxes=1;

% Load CF data using :-
%read_calipso_monthly_IPSL_2007_2017

% var_ukesm = 'Nd_cf_weighted_UKESM';
% var_ukesm = 'calipso_low_cloud_amount';
% var_ukesm = 'SW_up_TOA';

yr_start_trend=1850;        
%yr_start_trend=1970;
%yr_start_trend=1985;
%yr_start_trend=1982;
%yr_start_trend=2003;
%yr_start_trend=2007;

%yr_start_trend=2001;
%yr_start_trend=2007;

yr_end_trend=2014;
%yr_end_trend=1970;

iscreen_sig=0; %Whether to screen for signficance
marker_size=1; %For non-signficant points
%marker_size=1; %For non-signficant points
iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;

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


%% Run the script to load the 9 member ensemble

%ACSIS_dat_trends_load_ensemble(var_UM); %switch to loading the .mat file
% where var_UM is the name of the variable in the UM files.


switch var_ukesm
    case 'Nd_cf_weighted_UKESM_ztop'
        %load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Nd_trends_ukesm.mat');
    case 'Nd_cf_weighted_UKESM'        
        load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
        load_file_PI = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_PI_PI_control_Nd_cf_weighted_UKESM.mat';
        dat_PI = load(load_file_PI,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');
    case 'calipso_low_cloud_amount'       
        load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat'; 
case 'calipso_total_cloud_amount'       
    load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_total_cloud_amount.mat';
end


dat_ukesm=load(load_file,'dat_annual_ens','years_ukesm_1d','dat_annual','gcm_Plon2D_UM','gcm_Plat2D_UM','gcm_Plon2D_edges_UM','gcm_Plat2D_edges_UM');


ioverride_LAT_plots=1;

gcm_Plon2D_UM = dat_ukesm.gcm_Plon2D_UM;
gcm_Plat2D_UM = dat_ukesm.gcm_Plat2D_UM;
gcm_Plon2D_edges_UM = dat_ukesm.gcm_Plon2D_edges_UM;
gcm_Plat2D_edges_UM = dat_ukesm.gcm_Plat2D_edges_UM;


%% Trend analysis (linear least squares fit) - model
    yr_start=yr_start_trend; yr_end=yr_end_trend;
    yr_start_trend_used=yr_start_trend; yr_end_trend_used=yr_end_trend;    
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end);    
    x = [yr_start:yr_end]';
%     y = Nd_annual(istart:iend,100,100);
%     [Nd_trend,bint,residuals,rint,stats] = regress(y,x);
%     
%     [coeffs,bint,residuals,rint,stats] = trend_linear_fit_nd_data(Nd_annual(istart:iend,:,:),1);
    
   [coeffs,t_trend] = trend_linear_fit_nd_data(x,dat_ukesm.dat_annual(istart:iend,:,:),1); 
   
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
    %itrend_not_sig = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    itrend_not_sig = find((1-T1)<=p_conf/100); %1-tailed
   
%% Map of linear MODEL trend  
    yr_start=yr_start_trend_used; yr_end=yr_end_trend_used;
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end);    
    dat_modis = squeeze(coeffs(2,:,:)); 
    if iscreen_sig==1
        switch region_choice
            case 'global'
                %dat_modis(itrend_not_sig)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end  
    
ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on 
% var_ukesm   
   
    var_UM = ['UKESM ' var_str ' trend of ensemble mean between y' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and y' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
    tit_str_clean = ['UKESM ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];    
    subtitle_str = tit_str_clean;
    add_str = [' ' units_str_trend];
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims);   
    %caxis([-2 2]);
    xlabel(hc,units_str_trend); %label the colour bar
   
    
    % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;    
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end    
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_UM(itrend_not_sig),gcm_Plat2D_UM(itrend_not_sig),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end
    
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
%% Map of MODEL trend using last 5 years minus first 5 as suggested by Ken - Spatial pattern looks very similar, but with slightly smaller trend
    nyears = 5;
    yr_start=yr_start_trend_used; yr_end=yr_end_trend_used;
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end); 
    start_mean = meanNoNan(dat_ukesm.dat_annual(istart:istart+nyears-1,:,:),1);
    year_start_mean = meanNoNan(dat_ukesm.years_ukesm_1d(istart:istart+nyears-1),1);
    year_end_mean = meanNoNan(dat_ukesm.years_ukesm_1d(iend-nyears+1:iend),1);    
    end_mean = meanNoNan(dat_ukesm.dat_annual(iend-nyears+1:iend,:,:),1);
    dat_modis = ( end_mean - start_mean ) / ( year_end_mean - year_start_mean );
    trend_model_5years = dat_modis;
    if iscreen_sig==1
        switch region_choice
            case 'global'
                %dat_modis(itrend_not_sig)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end  
    
ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on 
% var_ukesm   
   
    var_UM = ['UKESM ' var_str ' end minus start trend of ensemble mean between y' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and y' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
    tit_str_clean = ['UKESM ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];    
    %run plotting script
    figure
    ioverride_proj_type=0;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=1;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-1 1]);   
    %caxis(clims);
    xlabel(hc,units_str_trend); %label the colour bar
   
    
    % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;    
    if i_plot_all_boxes==1
        ACSIS_Robson_paper_plot_all_boxes
    else
        plot_box_on_map
    end    
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_UM(itrend_not_sig),gcm_Plat2D_UM(itrend_not_sig),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end
    
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    

    

%% Map of MODEL time means
    switch var_ukesm
        case 'Nd_cf_weighted_UKESM'
            clims = [0 250];
            cbar_label = 'cm^{-3}'; %label the colour bar
            yr_start=2003; yr_end=2014;
        case 'calipso_low_cloud_amount'
            clims = [0 0.8];            
            cbar_label = ''; %label the colour bar
            yr_start=2007; yr_end=2014;
        case 'calipso_total_cloud_amount'
            clims = [0 1];            
            cbar_label = ''; %label the colour bar
            yr_start=2007; yr_end=2014;     
            subtitle_str = 'UKESM1 total cloud fraction';
        case 'SW_up_TOA'
            clims = [0 150];
            cbar_label = 'W m^{-2}'; %label the colour bar
            yr_start=2001; yr_end=2014;
    end
    
   
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end);    
    dat_modis = meanNoNan(dat_ukesm.dat_annual(istart:iend,:,:),1); 
    model_map = dat_modis;
     
    
%ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on 
% var_ukesm   
   
    %var_UM = ['UKESM ' var_ukesm ' time-mean between y' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and y' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
    %tit_str_clean = ['UKESM ' var_ukesm ' time-mean ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];        
    %subtitle_str = [subtitle_str ' ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];        
    subtitle_str = ['UKESM time-mean ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];
    add_str = [' ' remove_character(var_ukesm,'_',' ')];
    var_UM = [subtitle_str add_str];
    
    
    %run plotting script
    figure    
    set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
    %set(gcf,'position',[5 30 1252 590]);
    set(gcf,'position',[5 30 500 590]);
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    ioverride_ticks_DRIVER=1;
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));    
    caxis(clims);   
    %increase_font_size_map_figures;   %This creates a gap between map and
    %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?    
    xlabel(hc,cbar_label); %label the colour bar
    title(subtitle_str);
    
   
    
    % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;    
  
    switch var_ukesm
        case 'Nd_cf_weighted_UKESM'
            %if i_plot_all_boxes==1
            ACSIS_Robson_paper_plot_all_boxes
            %else
            %plot_box_on_map
            %end
        case 'calipso_low_cloud_amount'
            '';
    end
   
    
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    save_map_plot_data
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
    
%% calculate some timeseries in the regional box
    ilat = find(gcm_Plat2D_UM(:,1)>LAT_val(1) & gcm_Plat2D_UM(:,1)<LAT_val(2));
    ilon = find(gcm_Plon2D_UM(1,:)>LON_val(1) & gcm_Plon2D_UM(1,:)<LON_val(2));
%     clear dat_time_mean
%     for iy=1:size(Nd_ukesm,1)
%         dat_time_mean(iy,:) = meanNoNan(meanNoNan(dat.dat_ukesm(iy,:,ilat,ilon),4),2);        
%     end

    if iscreen_land==1
        load_type = 'merged netCDF';
        var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
        um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
        dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
        dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
              
        
        land_mask_mean = squeeze(ones(size(dat_ukesm.dat_annual(1,:,:))));
        land_mask_mean(dat_global.dat==1)=NaN;
                
        land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
        land_mask_ens = permute(land_mask_ens,[3 1 2]);
        
    else
        land_mask_mean = squeeze(ones(size(dat_ukesm.dat_annual(1,:,:))));
        land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
        land_mask_ens = permute(land_mask_ens,[3 1 2]);
                     
    end
    
  
    clear dat_annual_box_ukesm dat_annual_box_ukesm_ens
    for it=1:size(dat_ukesm.dat_annual,1)
        dat_mean = squeeze(dat_ukesm.dat_annual(it,:,:)) .* land_mask_mean;        
        dat_mean_ens = squeeze(dat_ukesm.dat_annual_ens(:,it,:,:)) .* land_mask_ens;
        dat_annual_box_ukesm(it) = meanNoNan(meanNoNan(dat_mean(ilat,ilon),2),1);        
        dat_tmp = meanNoNan(meanNoNan(dat_mean_ens(:,ilat,ilon),3),2);
        dat_annual_box_ukesm_ens_std(it) = std(dat_tmp); %std dev across the ensemble
    end
    
    land_mask_mean_all_t = repmat(land_mask_mean,[1 1 size(dat_PI.dat_annual,1)]);
    land_mask_mean_all_t = permute(land_mask_mean_all_t,[3 1 2]);
    dat_PI_box = meanNoNan(meanNoNan(dat_PI.dat_annual(:,ilat,ilon) .* land_mask_mean_all_t(:,ilat,ilon),3),2);
    [me_t_PI,N_t_PI,std_t_PI] = meanNoNan(dat_PI_box,1);

%% Obs data - load data and calc regional means, etc.
 % -- Think should use ACSIS_Robson_paper_load_data rather than this script
 % for loading in obs.
switch var_ukesm
    case 'SW_up_TOA'
%         switch SW_dataset
%             case 'CERES'
                
                %ceres_file = '/home/disk/eos15/d.grosvenor/eos8/CERES/ACSIS/CERES_EBAF-TOA_Ed4.0_Subset_200903-201004.nc';
                ceres_file = '/home/disk/eos15/d.grosvenor/eos8/CERES/ACSIS/CERES_EBAF_Ed4.1_Subset_200003-201905.nc';
                nc_ceres = netcdf(ceres_file);
                lon_ceres = nc_ceres{'lon'}(:);
                lat_ceres = nc_ceres{'lat'}(:);
                time_ceres = nc_ceres{'time'}(:); %days since 2000-03-01
                    %Runs from 15th March 2001 to 15th May 2009
                SW_TOA_ceres_in = nc_ceres{'toa_sw_all_mon'}(:);
                %solar_mon:long_name = "Incoming Solar Flux, Monthly Means" ;
                %SW_TOA_in_ceres = nc_ceres{'solar_mon'}(:);
                %cf_ceres = nc_ceres{'cldarea_total_daynight_mon'}(:);
                %tau_ceres = nc_ceres{'cldtau_total_day_mon'}(:);
                
                time_matlab_CERES = datenum('01-Mar-2000') + time_ceres;
                %[out, time_out_CERES, time_inds_CERES, dtime_match] = get_time_range_of_array(array_in,time_matlab_CERES,time_choice,dim);
                
                i180=find(lon_ceres>180);
                lon_ceres(i180) = lon_ceres(i180)-360;
                [gcm_Plon2D_CERES, gcm_Plat2D_CERES] = meshgrid(lon_ceres,lat_ceres);
                
                %Ceres data here runs from 15th March 2001 to 15th May 2009
                %- so select 2002 to 2018
                years_obs=[2002:2018];
                [Y,M,D] = datevec(time_matlab_CERES);
                istart = find(Y==years_obs(1));
                iend = find(Y==years_obs(end));
                SW_TOA_ceres = SW_TOA_ceres_in(istart(1):iend(end),:,:);
                [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(SW_TOA_ceres,gcm_Plat2D_CERES,gcm_Plon2D_CERES,LAT_val,LON_val,years_obs,gcm_area_CERES);
                
            %case 'DEEP-C'                                
                dat_file = ['/home/disk/eos15/d.grosvenor/eos8/DEEPC_CERES/sw_up_TOA_calculated.mat'];
                dat_file = ['/home/disk/eos15/d.grosvenor/eos8/DEEPC_CERES/sw_up_TOA_calculated_non_coarse_grained.mat'];
                %dat_file = ['/home/disk/eos15/d.grosvenor/eos8/DEEPC_CERES/sw_up_TOA_calculated_non_coarse_grained_non_GEOD.mat'];
                
                % Deep-C data file created using SW_DEEPC_data_read_FUNC.
                dat_deepc = load(dat_file,'time_matlab','sw_up_toa');              
                                                                
                years_obs_DeepC=[1985:2014];
                [obs_monthly_box_DeepC,obs_annual_box_DeepC,obs_annual_map_DeepC] = ACSIS_Robson_paper_process_obs(dat_deepc.sw_up_toa,gcm_Plat2D_UM,gcm_Plon2D_UM,LAT_val,LON_val,years_obs_DeepC,gcm_area_UM);
                
        %end
        
        
    
    case 'calipso_low_cloud_amount';
        obs_str='CALIPSO';
        
        
        %Load Calipso data - currently 2007-2017
        % Load CF data using :-
        read_calipso_monthly_IPSL_2007_2017
        years_obs = years_requested;
        
        [gcm_area_CALIPSO_monthly] = calc_area_lat_lon2d(gcm_Plat2D_edges_CALIPSO_monthly,gcm_Plon2D_edges_CALIPSO_monthly);
        [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(0.01*cllcalipso_monthly_AVERAGE,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,LAT_val,LON_val,years_obs,gcm_area_CALIPSO_monthly);
        
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
%         


    case 'calipso_total_cloud_amount';
        obs_str='CALIPSO';
                
        %Load Calipso data - currently 2007-2017
        % Load CF data using :-
        read_calipso_monthly_IPSL_2007_2017        
        years_obs = years_requested;
        
        [gcm_area_CALIPSO_monthly] = calc_area_lat_lon2d(gcm_Plat2D_edges_CALIPSO_monthly,gcm_Plon2D_edges_CALIPSO_monthly);        
        [obs_monthly_box,obs_annual_box,obs_annual_map] = ACSIS_Robson_paper_process_obs(0.01*cltcalipso_monthly_AVERAGE,gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,LAT_val,LON_val,years_obs,gcm_area_CALIPSO_monthly);        

    case 'CCI_total_cloud'
        obs_str='ESA CCI';
        cci_dat = load('/home/disk/eos8/d.grosvenor/ESA_Cloud_CCI/AVHRR_PMv3_L3C_Monthly/ESA_Cloud_CCI_Monthly_Cloud_Fraction.mat');
        
        gcm_Plat2D_CCI = cci_dat.gcm_Plat2D_CCI;
        gcm_Plon2D_CCI = cci_dat.gcm_Plon2D_CCI;
        [gcm_Plat2D_edges_CCI,gcm_Plon2D_edges_CCI] = get_edges_lat_lon(gcm_Plat2D_CCI,gcm_Plon2D_CCI);
        [gcm_area_CCI] = calc_area_lat_lon2d(gcm_Plat2D_edges_CCI,gcm_Plon2D_edges_CCI);        

        %For the map plots :-
        [obs_monthly_box_cci_TEMP,obs_annual_box_cci_TEMP,obs_annual_map_cci] = ACSIS_Robson_paper_process_obs(dat_tmp,gcm_Plat2D_CCI,gcm_Plon2D_CCI,LAT_val,LON_val,years_obs_cci,gcm_area_CCI);
        
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
                %dataset_str = 'SZA_LT_65_CF_GT_0_CTH_LT_3.2km_screened_for_seaice__2week_max';
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
        
        [mon_me_region_MODIS2] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val,iscreen_land,lmask_MODIS,0,1);
        [mon_me_region_MODIS2_Ndatap] = monthly_means_restrict_lat_lon(years_MODIS2,mon_me_MODIS2_Ndatap,gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,LAT_val,LON_val,iscreen_land,lmask_MODIS,0,1);
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

%% Trend analysis (linear least squares fit) - MODIS
    yr_start_MODIS=2003; yr_end_MODIS=2014;
    %yr_start_MODIS=2007; yr_end_MODIS=2014;
    istart_MODIS=find(years_MODIS2==yr_start_MODIS);
    iend_MODIS=find(years_MODIS2==yr_end_MODIS);                
    x = [yr_start_MODIS:yr_end_MODIS]';
    %The following takes a litte while...
   [coeffs_MODIS,t_trend_MODIS] = trend_linear_fit_nd_data(x,Nd_annual_MODIS(istart_MODIS:iend_MODIS,:,:),1); 
   
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
    %itrend_not_sig_MODIS = find(abs(t_trend_MODIS)<=t_thresh); %make NaN for now, but can put a dot on, etc.
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_MODIS),n_dof);
    T1 = T2/2;    
    itrend_not_sig_MODIS = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig_MODIS = find((1-T1)<=p_conf/100); %1-tailed

    
    %% Plot MODIS trend map
    dat_modis2 = squeeze(coeffs_MODIS(2,:,:));
%     if iscreen_land==1
%         	land_mask=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%             lmask = flipdim(land_mask.amsre_land_mask,1);
%             dat_modis2(isnan(lmask)==1)=NaN;            
%     end

    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' N_d trend between y' num2str(years_MODIS2(istart_MODIS)) ' and y' num2str(years_MODIS2(iend_MODIS)) '; (cm^{-3} yr^{-1})' add_str];
    
    tit_str_clean = ['MODIS N_d trend ' num2str(years_MODIS2(istart_MODIS)) ' to ' num2str(years_MODIS2(iend_MODIS))];
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);
    xlabel(hc,units_str_trend); %label the colour bar
    
   % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end 
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_AMSRE(itrend_not_sig_MODIS),gcm_Plat2D_AMSRE(itrend_not_sig_MODIS),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end

    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
%% Plot MODIS time mean Nd
    
    dat_modis2 = meanNoNan(Nd_annual_MODIS,1);
    %dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
    
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    obs_map = dat_modis;
    var_UM = ['MODIS time-mean ' str_label_2137 ' ' cf_screen_str ' Nd for ' num2str(years_MODIS2(1)) ' to ' num2str(years_MODIS2(end)) '; (cm^{-3} yr^{-1})'];
    tit_str_clean = ['MODIS N_d time-mean ' num2str(years_MODIS2(1)) ' to ' num2str(years_MODIS2(end))];
    subtitle_str =  tit_str_clean;
    %run plotting script
    figure
    set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
    %set(gcf,'position',[5 30 1252 590]);
    set(gcf,'position',[5 30 500 590]);
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 250]);
     %increase_font_size_map_figures;   %This creates a gap between map and
    %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?    
    xlabel(hc,'cm^{-3}'); %label the colour bar
    title(subtitle_str);
    
%     % Plot the box    
%     ioverride_box_colour=1;  
%     irotated_pole_box=0;    
%     itext_in_box=0;
%     imap=1;
%     col_str='k-';
%     box_lwidth = 3;
%     %i_plot_all_boxes=1;
%     %if i_plot_all_boxes==1
%         ACSIS_Robson_paper_plot_all_boxes
%     %else
%         %plot_box_on_map
%     %end 
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
%% Plot model bias 
    
    dat_modis = 100*(model_map - obs_map) ./ obs_map;
    %var_UM = ['UKESM prc Nd bias for ' num2str(years_MODIS2(1)) ' to ' num2str(years_MODIS2(end))];
    %tit_str_clean = ['UKESM bias ' num2str(years_MODIS2(1)) ' to ' num2str(years_MODIS2(end))];    
    subtitle_str = ['UKESM N_d bias ' num2str(years_MODIS2(1)) ' to ' num2str(years_MODIS2(end))];    
    add_str = [' ' remove_character(var_ukesm,'_',' ')];
    var_UM = [remove_character(var_ukesm,'_',' ') subtitle_str];
    
    %run plotting script
    figure
    set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
    %set(gcf,'position',[5 30 1252 590]);
    set(gcf,'position',[5 30 500 590]);
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);
     %increase_font_size_map_figures;   %This creates a gap between map and
    %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?    
    %xlabel(hc,'cm^{-3}'); %label the colour bar
    xlabel(hc,'%'); %label the colour bar
    title(subtitle_str);
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);    
    
%% Plot CALIPSO time mean CF
    yr_start_CALIPSO=2007; yr_end_CALIPSO=2014;
    istart_CALIPSO=find(years_obs==yr_start_CALIPSO);
    iend_CALIPSO=find(years_obs==yr_end_CALIPSO);    
    
    dat_modis2 = meanNoNan(obs_annual_map(istart_CALIPSO:iend_CALIPSO,:,:),1);
    %dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
    
    dat_modis =  griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    obs_map = dat_modis;
    %var_UM = ['CALIPSO time-mean ' var_ukesm ' for ' num2str(yr_start_CALIPSO) ' to ' num2str(yr_end_CALIPSO) ';'];
        %The above used for the filename, so keep this, or base it on
        %subtitle_str and add_str
    %tit_str_clean = ['CALIPSO time-mean ' num2str(yr_start_CALIPSO) ' to ' num2str(yr_end_CALIPSO)];
    %tit_str_clean = var_UM;
    add_str = [' ' remove_character(var_ukesm,'_',' ')];
    subtitle_str = ['CALIPSO time-mean ' num2str(yr_start_CALIPSO) ' to ' num2str(yr_end_CALIPSO)]
    var_UM = [subtitle_str add_str];
    
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    switch var_ukesm
        case 'calipso_low_cloud_amount'
            caxis([0 0.8]);
        case 'calipso_total_cloud_amount'
            caxis([0 1]);
    end
    %xlabel(hc,'cm^{-3}'); %label the colour bar
    
%     % Plot the box    
%     ioverride_box_colour=1;  
%     irotated_pole_box=0;    
%     itext_in_box=0;
%     imap=1;
%     col_str='k-';
%     box_lwidth = 3;
%     %i_plot_all_boxes=1;
%     %if i_plot_all_boxes==1
%         ACSIS_Robson_paper_plot_all_boxes
%     %else
%         %plot_box_on_map
%     %end 
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    save_map_plot_data
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
 
    
%% Plot model CF bias 
    
    dat_modis = 100*(model_map - obs_map) ./ obs_map;
    var_UM = ['UKESM prc CF bias for ' num2str(yr_start_CALIPSO) ' to ' num2str(yr_end_CALIPSO) ';'];
    %tit_str_clean = ['UKESM bias ' num2str(yr_start_CALIPSO) ' to ' num2str(yr_end_CALIPSO)];
    subtitle_str = ['Model bias']
    add_str = [' ' remove_character(var_ukesm,'_',' ')];
    
    
    
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-50 50]);  
    xlabel(hc,'%'); %label the colour bar
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    save_map_plot_data
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);    
        

%% Plot CERES SW TOA time mean
    yr_start_CERES=2001; yr_end_CERES=2014;    
    istart_CERES=find(years_obs_CERES==yr_start_CERES);
    iend_CERES=find(years_obs_CERES==yr_end_CERES);    
    
    dat_modis2 = meanNoNan(obs_annual_map(istart_CERES:iend_CERES,:,:),1);
    %dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
    
    dat_modis =  griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    obs_map = dat_modis;
    var_UM = ['CERES time-mean for ' num2str(yr_start_CERES) ' to ' num2str(yr_end_CERES) ';'];
    tit_str_clean = ['CERES time-mean ' num2str(yr_start_CERES) ' to ' num2str(yr_end_CERES)];
    subtitle_str = tit_str_clean;
    
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 150]);  
    xlabel(hc,'W m^{-2}'); %label the colour bar
    
%     % Plot the box    
%     ioverride_box_colour=1;  
%     irotated_pole_box=0;    
%     itext_in_box=0;
%     imap=1;
%     col_str='k-';
%     box_lwidth = 3;
%     %i_plot_all_boxes=1;
%     %if i_plot_all_boxes==1
%         ACSIS_Robson_paper_plot_all_boxes
%     %else
%         %plot_box_on_map
%     %end 
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
 
    
%% Plot model SW bias 
    %dat_modis = model_map - obs_map;
    dat_modis = 100*(model_map - obs_map) ./ obs_map;
    var_UM = ['UKESM prc SWTOA bias for ' num2str(yr_start_CERES) ' to ' num2str(yr_end_CERES) ';'];
    tit_str_clean = ['UKESM bias ' num2str(yr_start_CERES) ' to ' num2str(yr_end_CERES)];
    
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-40 40]);  
    xlabel(hc,'%'); %label the colour bar
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);    
    
    
    
    
 %% Plot MODIS Nd for 2006
    istart=find(years_MODIS2==2006);
    Nd_annual_MODIS(istart,:,:);
    dat_modis2 = squeeze(Nd_annual_MODIS(istart,:,:));
    %dat_modis2(itrend_not_sig_MODIS)=NaN; %make NaN for now, but can put a dot on, etc.
    
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd for y' num2str(years_MODIS2(istart)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([0 300]);



%% plot MODIS difference map 2003 to 2014
%2003 to 2014    
    yr_start=2003; yr_end=2014;
    istart=find(years_MODIS2==yr_start);
    iend=find(years_MODIS2==yr_end);
    dat = squeeze (Nd_annual_MODIS(iend,:,:) - Nd_annual_MODIS(istart,:,:) ) / (iend-istart);   
    dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat,gcm_Plat2D_UM,gcm_Plon2D_UM);
    var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd trend between y' num2str(years_MODIS2(istart)) ' and y' num2str(years_MODIS2(iend)) '; (cm^{-3} yr^{-1})'];
    %run plotting script
    figure
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-6 6]);   
    %caxis([-2 2]);
    
        
% Plot the box    
    ioverride_box_colour=1;    
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        ACSIS_Robson_paper_plot_all_boxes
    else
        plot_box_on_map
    end 
    
    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

%% Trend analysis (linear least squares fit) - CALIPSO low cloud fraction
    yr_start_CALIPSO=2007; yr_end_CALIPSO=2014;
    istart_CALIPSO=find(years_obs==yr_start_CALIPSO);
    iend_CALIPSO=find(years_obs==yr_end_CALIPSO);                
    x = [yr_start_CALIPSO:yr_end_CALIPSO]';
    %The following takes a litte while...
   [coeffs_CALIPSO,t_trend_CALIPSO] = trend_linear_fit_nd_data(x,obs_annual_map(istart_CALIPSO:iend_CALIPSO,:,:),1); 
   
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
    %itrend_not_sig_MODIS = find(abs(t_trend_MODIS)<=t_thresh); %make NaN for now, but can put a dot on, etc.
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_CALIPSO),n_dof);
    T1 = T2/2;    
    itrend_not_sig_CALIPSO = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    %itrend_not_sig_CALIPSO = find((1-T1)<=p_conf/100); %1-tailed    
    
    
    
%% Plot CALIPSO trend map
    dat_modis2 = squeeze(coeffs_CALIPSO(2,:,:));
%     if iscreen_land==1
%         	land_mask=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%             lmask = flipdim(land_mask.amsre_land_mask,1);
%             dat_modis2(isnan(lmask)==1)=NaN;            
%     end

    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_CALIPSO)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    dat_modis =  griddata(gcm_Plat2D_CALIPSO_monthly,gcm_Plon2D_CALIPSO_monthly,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    
    
    var_UM = ['CALIPSO cf trend between y' num2str(years_obs(istart_CALIPSO)) ' and y' num2str(years_obs(iend_CALIPSO)) '; (cm^{-3} yr^{-1})' add_str];
    
    tit_str_clean = ['CALIPSO f_c trend ' num2str(years_obs(istart_CALIPSO)) ' to ' num2str(years_obs(iend_CALIPSO))];
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-0.01 0.01]);
    xlabel(hc,units_str_trend); %label the colour bar
    xlabel(hc,'yr^{-1}'); %label the colour bar
    
   % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end 
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
                otherwisexs
                m_plot(gcm_Plon2D_CALIPSO_monthly(itrend_not_sig_CALIPSO),gcm_Plat2D_CALIPSO_monthly(itrend_not_sig_CALIPSO),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end

    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
    
    
%% Trend analysis (linear least squares fit) - CERES SWTOA
    %yr_start_CERES=2001; yr_end_CERES=2014;
    yr_start_CERES=2007; yr_end_CERES=2014;
    istart_CERES=find(years_obs==yr_start_CERES);
    iend_CERES=find(years_obs==yr_end_CERES);                
    x = [yr_start_CERES:yr_end_CERES]';
    %The following takes a litte while...
   [coeffs_CERES,t_trend_CERES] = trend_linear_fit_nd_data(x,obs_annual_map(istart_CERES:iend_CERES,:,:),1); 
   
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
    %itrend_not_sig_MODIS = find(abs(t_trend_MODIS)<=t_thresh); %make NaN for now, but can put a dot on, etc.
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_CERES),n_dof);
    T1 = T2/2;    
    
    
    switch tail_test
        case 1
            itrend_not_sig_CERES = find((1-T1)<=p_conf/100); %1-tailed
        case 2
            itrend_not_sig_CERES = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    end
    
    
    
    
%% Plot CERES trend map
    dat_modis2 = squeeze(coeffs_CERES(2,:,:));
%     if iscreen_land==1
%         	land_mask=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%             lmask = flipdim(land_mask.amsre_land_mask,1);
%             dat_modis2(isnan(lmask)==1)=NaN;            
%     end

    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_CERES)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    dat_modis =  griddata(gcm_Plat2D_CERES,gcm_Plon2D_CERES,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    
    
    var_UM = ['CERES SWTOA trend between y' num2str(years_obs(istart_CERES)) ' and y' num2str(years_obs(iend_CERES)) '; (cm^{-3} yr^{-1})' add_str];
    
    tit_str_clean = ['CERES SWTOA trend ' num2str(years_obs(istart_CERES)) ' to ' num2str(years_obs(iend_CERES))];
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis([-1 1]);
    xlabel(hc,units_str_trend); %label the colour bar
    %xlabel(hc,'W m^{-2} yr^{-1}'); %label the colour bar
    
   % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end 
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_CERES(itrend_not_sig_CERES),gcm_Plat2D_CERES(itrend_not_sig_CERES),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end

    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
%% Trend analysis (linear least squares fit) - Deep-C SWTOA
    %yr_start_DeepC=2001; yr_end_DeepC=2014;
    yr_start_DeepC=1985; yr_end_DeepC=2014;
    istart_DeepC=find(years_obs_DeepC==yr_start_DeepC);
    iend_DeepC=find(years_obs_DeepC==yr_end_DeepC);                
    x = [yr_start_DeepC:yr_end_DeepC]';
    %The following takes a litte while...
   [coeffs_DeepC,t_trend_DeepC] = trend_linear_fit_nd_data(x,obs_annual_map_DeepC(istart_DeepC:iend_DeepC,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
    T2 = 2*tcdf(-abs(t_trend_DeepC),n_dof);
    T1 = T2/2;    
    
    
    switch tail_test
        case 1
            itrend_not_sig_DeepC = find((1-T1)<=p_conf/100); %1-tailed
        case 2
            itrend_not_sig_DeepC = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    end
    
    
    
    
%% Plot DeepC trend map
    dat_modis2 = squeeze(coeffs_DeepC(2,:,:));
%     if iscreen_land==1
%         	land_mask=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%             lmask = flipdim(land_mask.amsre_land_mask,1);
%             dat_modis2(isnan(lmask)==1)=NaN;            
%     end

    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_DeepC)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    %dat_modis =  griddata(gcm_Plat2D_DeepC,gcm_Plon2D_DeepC,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    dat_modis = dat_modis2;
    
    var_UM = ['DeepC SWTOA trend between y' num2str(years_obs_DeepC(istart_DeepC)) ' and y' num2str(years_obs_DeepC(iend_DeepC)) '; (cm^{-3} yr^{-1})' add_str];    
    tit_str_clean = ['DeepC SWTOA trend ' num2str(years_obs_DeepC(istart_DeepC)) ' to ' num2str(years_obs_DeepC(iend_DeepC))];
    subtitle_str = tit_str_clean;

    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims);
    xlabel(hc,units_str_trend); %label the colour bar
    %xlabel(hc,'W m^{-2} yr^{-1}'); %label the colour bar
    
   % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end 
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                m_plot(gcm_Plon2D_DeepC(itrend_not_sig_DeepC),gcm_Plat2D_DeepC(itrend_not_sig_DeepC),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end

    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
       
 
%% Trend analysis (linear least squares fit) - CCI Cloud Fraction
    %yr_start_cci=2001; yr_end_cci=2014;
    yr_start_cci=1982; yr_end_cci=2014;
    yr_start_cci=2007; yr_end_cci=2014;
    yr_start_cci=1985; yr_end_cci=2014;
    istart_cci=find(years_obs_cci==yr_start_cci);
    iend_cci=find(years_obs==yr_end_cci);                
    x = [yr_start_cci:yr_end_cci]';
    %The following takes a litte while...
   [coeffs_cci,t_trend_cci] = trend_linear_fit_nd_data(x,obs_annual_map_cci(istart_cci:iend_cci,:,:),1); 
   
   %t-test threshold for 95% confidence
    n_dof = length(x)-2; %number of degrees of freedom
   
    % Significance of our t values :-
    T2 = 2*tcdf(-abs(t_trend_cci),n_dof);
    T1 = T2/2;    
    
    
    switch tail_test
        case 1
            itrend_not_sig_cci = find((1-T1)<=p_conf/100); %1-tailed
        case 2
            itrend_not_sig_cci = find((1-T2)<=p_conf/100); %2-tailed test - is more appropriate I think since trend can be positive or neg
    end
    
    
    
    
%% Plot CCI trend map
    ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on 
% var_ukesm   

    dat_modis2 = squeeze(coeffs_cci(2,:,:));
%     if iscreen_land==1
%         	land_mask=load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%             lmask = flipdim(land_mask.amsre_land_mask,1);
%             dat_modis2(isnan(lmask)==1)=NaN;            
%     end

    if iscreen_sig==1
        switch region_choice
            case 'global'
                dat_modis2(itrend_not_sig_cci)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end    
    
    dat_modis =  griddata(gcm_Plat2D_CCI,gcm_Plon2D_CCI,dat_modis2,gcm_Plat2D_UM,gcm_Plon2D_UM);
    
    
    var_UM = ['CCI total cloud fraction trend between y' num2str(years_obs(istart_cci)) ' and y' num2str(years_obs(iend_cci)) '; (cm^{-3} yr^{-1})' add_str];        
    tit_str_clean = ['CCI total cloud fraction trend ' num2str(years_obs(istart_cci)) ' to ' num2str(years_obs(iend_cci))];
    subtitle_str = tit_str_clean;
    %run plotting script
    figure
    ioverride_proj_type=1;
    proj_type_DRIVER='ortho';
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims);
    xlabel(hc,units_str_trend); %label the colour bar
    %xlabel(hc,'W m^{-2} yr^{-1}'); %label the colour bar
    
   % Plot the box    
    ioverride_box_colour=1;  
    irotated_pole_box=0;    
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if i_plot_all_boxes==1
        %ACSIS_Robson_paper_plot_all_boxes
    else
        %plot_box_on_map
    end 
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                %m_plot(gcm_Plon2D_CCI(itrend_not_sig_cci),gcm_Plat2D_CCI(itrend_not_sig_cci),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end

    
    savename=[savedir_date var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    
    
       
    
    
%% Linear trend calc for the annual means over the box region.
yr_start=2003; yr_end=2014;
istart=find(years_ukesm_1d==yr_start);
iend=find(years_ukesm==yr_end);
y = Nd_annual_box_ukesm(istart:iend);
x1 = [yr_start:yr_end];
x = [ones(size(y)); x1];
[Nd_trend,bint,residuals,rint,stats] = regress(y',x');
ylin=Nd_trend(1)+Nd_trend(2).*x1; %this is the linear trend fit line

%% plot annual mean timeseries of model vs obs for 1850 to 2014
titlenam = ['Annual mean ' var_str  ' for region ' box_region_str];
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
patch(x_ens_patch,y_ens_patch,[0.5 0.5 1],'linestyle','none','HandleVisibility','Off');
hold on
%h1=plot(dat.years_ukesm_1d,dat_annual_box_ukesm,'bo-');
h1=plot(dat_ukesm.years_ukesm_1d,dat_annual_box_ukesm,'b-');

%errorbar(1980,me_t_PI,std_t_PI*2,'CapSize',2,'ro','markerfacecolor','r');
herr=errorbarYY('vert',1847.5,me_t_PI,std_t_PI*2,gca,'k','o',2,0.01);
set(herr,'HandleVisibility','Off');
%errorbarYY('vert',1960,me_t_PI-std_t_PI,std_t_PI,gca,'r','o',2,0.01);

%leg_str{ileg}='\pm 2 std. dev.'; ileg=ileg+1;
leg_str{ileg}='UKESM'; ileg=ileg+1;
set(h1,'linewidth',4);
set(h1,'markerfacecolor','b');
hold on


%h2=plot(years_obs,obs_annual_box,'rs-');
h2=plot(years_obs,obs_annual_box,'r-');
%leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
%leg_str{ileg}=['CALIPSO'];
leg_str{ileg}=obs_str;
set(h2,'linewidth',4);
set(h2,'markerfacecolor','r');
loc='SouthWest';
loc='NorthEast';

%set(gca,'xlim',[1980 2018]);
set(gca,'xlim',[1844 2018]);

set(gca,'ylim',[45 230]);

legend(leg_str,'location',loc);
xlabel('Year');
ylabel([var_str ' ' units_str]);

%switch box_region_DRIVER
%    case '1'
%        title('');
%        title(titlenam,'position',[1925 200]);
%    otherwise
%        title(titlenam,'position',[1925 200]);
%end
%title(titlenam,'position',[1925 200]);
title(tit_str_clean,'position',[1945 200]);

grid on


savename=[savedir_date titlenam];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;

%saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

    

%% plot INSET annual mean timeseries of model vs obs for 2003 to 2014
titlenam = ['Annual mean ' var_str  ' for region ' box_region_str];
tit_str_clean = ['Region ' box_region_str];
%figure

%ax2 = axes('Position',[0.17 0.7 0.25 0.25]);
%ax2 = axes('Position',[0.1748    0.7683    0.2500    0.2500]);
ax2 = axes('Position',inset_ax_pos);

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

set(gca,'xlim',[2002 2015]);
%set(gca,'xlim',[1850 2018]);

%set(gca,'ylim',[100 175]);
%set(gca,'ylim',[45 210]);
%set(gca,'ylim',[60 170]);
set(gca,'ylim',ylims_inset);


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
    



%% plot annual mean timeseries with MODIS
titlenam = 'Annual mean N_d';
figure
clear leg_str
ileg=1;
h1=errorbar(years_ukesm_1d,Nd_annual_box_ukesm,Nd_annual_box_ukesm_ens_std*2,'bo-');
leg_str{ileg}='UKESM'; ileg=ileg+1;
set(h1,'linewidth',3);
set(h1,'markerfacecolor','b');
hold on


h2=plot(years_MODIS2,Nd_annual_box_MODIS,'rs-');
leg_str{ileg}=['MODIS ' str_label_2137 ' ' cf_screen_str]; ileg=ileg+1;
set(h2,'linewidth',3);
set(h2,'markerfacecolor','r');
loc='SouthWest';
loc='NorthEast';

set(gca,'xlim',[1999 2015]);
set(gca,'XMinorTick','off');
fontsize_figure(gcf,gca,18);

legend(leg_str,'location',loc);
xlabel('Year');
ylabel('N_d (cm^{-3})');
title(titlenam);
grid on

savename=[savedir_date titlenam ' ' str_label_2137 ' ' cf_screen_str];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);



%% Find a correction for the MODIS data to adjust for the fact that we are using y2000 emissions, but comparing to y2009.
% Will apply the correction to tbe obs - so need to take the obs back to y2000
% So, if the trend is positive then need to reduce the obs and so dNd will
% be negative
MODIS_trends = squeeze(coeffs_MODIS(2,:,:));
MODIS_trends(itrend_not_sig_MODIS)=0; %make the non-significant trends zero
ny_corr = 2009-2000; %number of years that are correcting for.
dNd_corr = -MODIS_trends*ny_corr;
dNd_corr_UM_grid = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dNd_corr,gcm_Plat2D_UM,gcm_Plon2D_UM); %ADD theses to MODIS values

save_file_dNd_trend = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/dNd_trend_correction.mat'];
save(save_file_dNd_trend,'dNd_corr','dNd_corr_UM_grid','-V7.3');
%ref year for percentage diff plot
yr=2009;
%yr=2010;
iyr=find(years_MODIS2==yr);
Nd_ref_year = squeeze (Nd_annual_MODIS(iyr,:,:));

prc_diff = 100*dNd_corr./Nd_ref_year;


dat_modis = dNd_corr_UM_grid;
%dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,prc_diff,gcm_Plat2D_UM,gcm_Plon2D_UM);

figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([-60 60]);
    
%%
dat_modis = griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,Nd_ref_year,gcm_Plat2D_UM,gcm_Plon2D_UM);
var_UM = ['MODIS ' str_label_2137 ' ' cf_screen_str ' Nd for y' num2str(years_MODIS2(iyr)) ' (cm^{-3})'];
figure
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 300]);


%% Make the 3 panel plot of Total CF evaluation vs CALIPSO
UM_Robson_SUBPLOT_commands_Fig_tot_CF_eval2

%% Map of variability of SWTOA for all Januarys

% Load data
%dat_all_ukesm = load(load_file,'dat_ens'); %size [9        1980         144         192]
%1980 is all 12 months of all 165 years.

dat_all_ukesm = load(load_file_PI,'dat_ens'); %size [1        9000         144         192]
%9000 is all 12 months for 750 years (PI control)

%% Extract and plot
% Pick one ensemble and all Januarys.
iens=1;
imonth=1;
dat_Jan  = squeeze( dat_all_ukesm.dat_ens(iens,imonth:12:end,:,:) );
[me,N,std] = meanNoNan(dat_Jan,1);


%% std dev
dat_modis = std;
subtitle_str = ['Interannual Std. Dev. of SWTOA for Januarys'];
add_str = [' ' remove_character(var_ukesm,'_',' ')];
var_UM = [subtitle_str add_str];
cbar_label=('W m^{-2}');

%run plotting script
figure
ioverride_proj_type=0;
proj_type_DRIVER='ortho';
ioverride_ticks_DRIVER=1;
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 30]);
xlabel(hc,cbar_label); %label the colour bar


savename=[savedir_date var_UM];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
save_map_plot_data
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);


%% Relative std dev
dat_modis = 100*std./me;
subtitle_str = ['Relative Interannual Std. Dev. (%) of SWTOA for Januarys'];
add_str = [' ' remove_character(var_ukesm,'_',' ')];
var_UM = [subtitle_str add_str];
cbar_label=('%');

%run plotting script
figure
ioverride_proj_type=0;
proj_type_DRIVER='ortho';
ioverride_ticks_DRIVER=1;
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 40]);
xlabel(hc,cbar_label); %label the colour bar


savename=[savedir_date var_UM];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
save_map_plot_data
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);



%% mean
dat_modis = me;
subtitle_str = ['Mean SWTOA for Januarys'];
add_str = [' ' remove_character(var_ukesm,'_',' ')];
var_UM = [subtitle_str add_str];
cbar_label=('W m^{-2}');

%run plotting script
figure
ioverride_proj_type=0;
proj_type_DRIVER='ortho';
ioverride_ticks_DRIVER=1;
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
%caxis(clims);
xlabel(hc,cbar_label); %label the colour bar


savename=[savedir_date var_UM];
clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
save_map_plot_data
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);



%----------------------------------------------------------------------
%----- SW calculation/recreation using model monthly means
%----------------------------------------------------------------------
%% 1) Load ensemble monthly mean data

%SW down TOA
 load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_SW_down_TOA_partial_SW_down_TOA.mat';
 load(load_file,'dat_ens_mean'); %size [380         144         192];
 %380 is 8 months of 1983 and all 12 months for 1984-2014
 %So ignore the 1st 8 months to make it simpler
 SW_down_TOA_dat_ens_mean = dat_ens_mean(9:end,:,:);  %1984-2014

 %Pick out just 1984-2014 from others (all months for 1850-2014)
 %(1984-1850)*12 = 1608 
 years_sw_calc = [1984:2014];
 t_inds = 1609:1980;
 
% Nd
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_Nd_cf_weighted_UKESM.mat';  
load(load_file,'dat_ens_mean'); %size [1980         144         192];
%1980 is all 12 months of all 165 years.
Nd_dat_ens_mean = dat_ens_mean(t_inds,:,:);

%total CF
%load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat'; 
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_total_cloud_amount.mat';
load(load_file,'dat_ens_mean'); %size [1980         144         192];
%1980 is all 12 months of all 165 years.
totcf_dat_ens_mean = dat_ens_mean(t_inds,:,:);

%low CF
%load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat'; 
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat';
load(load_file,'dat_ens_mean'); %size [1980         144         192];
%1980 is all 12 months of all 165 years.
lowcf_dat_ens_mean = dat_ens_mean(t_inds,:,:);

%SW TOA up (to compare against the calculation)
%load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_calipso_low_cloud_amount.mat'; 
load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_SW_up_TOA.mat';
load(load_file,'dat_ens_mean'); %size [1980         144         192];
%1980 is all 12 months of all 165 years.
SW_up_TOA_dat_ens_mean = dat_ens_mean(t_inds,:,:);

%LWP (kg/m2)
 load_file = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_timeseries_ukesm_all_LWP_2-391.mat';
 load(load_file,'dat_ens_mean'); %size [1980         144         192];
 %1980 is all 12 months of all 165 years.
 lwp_dat_ens_mean = dat_ens_mean(t_inds,:,:);

%% 2) Do SW TOA calcs from monthly means - likely inaccuracies due to using monthly mean SWTOA down, etc.
tr_for_SW_bias = 'constant value';
%tr_for_SW_bias = 'time mean from model';

tr_atmos_constant_val = 0.6;
%tr_atmos_constant_val = 0.75;

cf_min = 0.01;
i_Liu = 1;

% Also need the clear-sky albedo - can estimate from cloud free regions?
A_clear = 0.15; %Guess for now - value for CF quite sensitive to this (ranging up to 30%), so should prob
  %check what this is from the model
A_clear = 0.12; %value estimated from the 12 noon snapshot from VOCALS reginal run - N.B. - this uses the estimated transmission since
 %it is the surface albedo assuming no atmosphere above
%A_clear = 0.0;  
A_clear = 0.09;

fac_tr = 1/1.1;
 fac_tr = 1;
 
 fac_SWin=1.7;
 fac_SWin=1;
 
 fac_Nd = 1;
 %fac_Nd = 2;%1.3;
 
 fac_cf = 1;
 %fac_cf = 1.2;
 
 fac_SW_calc=1;
 %fac_SW_calc=1.4;       
  
 
 switch tr_for_SW_bias
     case 'constant value'
        tr_atmos_mean = tr_atmos_constant_val; %if constant value
     case 'time mean from model'
         SW_thresh=200;
         transmission_atmos = SW_clean_clear_PD_ALL ./ SW_in_PD_ALL;
         transmission_atmos(SW_in_PD_ALL<SW_thresh)=0.5;
         tr_atmos_mean = meanNoNan(transmission_atmos,3);
 end
 
 
 
%% Do actual calc
time_mean=0;
time_mean=1;

%cf = lowcf_dat_ens_mean;
cf = totcf_dat_ens_mean;
cf2 = cf;
cfmin_lwp = cf_min;
%cfmin_lwp = 0.05;
cf2(cf<cfmin_lwp)=cfmin_lwp;

lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%lwp = lwp_dat_ens_mean

nd = Nd_dat_ens_mean;

switch tr_for_SW_bias
    case 'constant value'
        trans = tr_atmos_constant_val; %if constant value
    case 'time mean from model'
        trans = tr_atmos_mean;
end


str_type_calc = 'time-mean';

%sw = meanNoNan(SW_in,3);
sw = fac_SWin * SW_down_TOA_dat_ens_mean .* trans;

[Ac_calc_model,tau_calc_model,A_calc_model,SWTOA_calc_model] = calc_SW(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu); 
SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual
for iy=1:length(years_sw_calc)
   istart=(iy-1)*12+1;
   SWTOA_calc_model_annual(iy,:,:) = meanNoNan(SWTOA_calc_model(istart:istart+11,:,:),1);       
end



%% Do sensitivity test SW calcs

%% Keep Cf constant

istart=1984-1850+1; %starting index to take average over to use as constant value

%Keep cf constant first

%cf = lowcf_dat_ens_mean;
cf = totcf_dat_ens_mean;
cf2 = cf;

mean_cf = meanNoNan(cf(istart:end,:,:),1);
cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
cf_const = permute(cf_const,[3 1 2]);
cf = cf_const;


cfmin_lwp = cf_min;
%cfmin_lwp = 0.05;
cf2(cf2<cfmin_lwp)=cfmin_lwp;

lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%lwp = lwp_dat_ens_mean

nd = Nd_dat_ens_mean;


[Ac_calc_model_cf,tau_calc_model_cf,A_calc_model_cf,SWTOA_calc_model_cf] = calc_SW(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu); 
SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual_cf
for iy=1:length(years_sw_calc)
   istart=(iy-1)*12+1;
   SWTOA_calc_model_annual_cf(iy,:,:) = meanNoNan(SWTOA_calc_model_cf(istart:istart+11,:,:),1);       
end


%% Keep Nd constant

%cf = lowcf_dat_ens_mean;
cf = totcf_dat_ens_mean;

cf2 = cf;
cfmin_lwp = cf_min;
%cfmin_lwp = 0.05;
cf2(cf<cfmin_lwp)=cfmin_lwp;

lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%lwp = lwp_dat_ens_mean

nd = Nd_dat_ens_mean;
mean_val = meanNoNan(nd(istart:end,:,:),1);
val_const = repmat(mean_val,[1 1 size(nd,1)]);
val_const = permute(val_const,[3 1 2]);
nd = val_const;


[Ac_calc_model_Nd,tau_calc_model_Nd,A_calc_model_Nd,SWTOA_calc_model_Nd] = calc_SW(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu); 
SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual_Nd
for iy=1:length(years_sw_calc)
   istart=(iy-1)*12+1;
   SWTOA_calc_model_annual_Nd(iy,:,:) = meanNoNan(SWTOA_calc_model_Nd(istart:istart+11,:,:),1);       
end


%% Keep LWP constant

%cf = lowcf_dat_ens_mean;
cf = totcf_dat_ens_mean;

cf2 = cf;
cfmin_lwp = cf_min;
%cfmin_lwp = 0.05;
cf2(cf2<cfmin_lwp)=cfmin_lwp;

lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%lwp = lwp_dat_ens_mean
mean_val = meanNoNan(lwp(istart:end,:,:),1);
val_const = repmat(mean_val,[1 1 size(lwp,1)]);
val_const = permute(val_const,[3 1 2]);
lwp = val_const;

nd = Nd_dat_ens_mean;

[Ac_calc_model_lwp,tau_calc_model_lwp,A_calc_model_lwp,SWTOA_calc_model_lwp] = calc_SW(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu); 
SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual_lwp
for iy=1:length(years_sw_calc)
   istart=(iy-1)*12+1;
   SWTOA_calc_model_annual_lwp(iy,:,:) = meanNoNan(SWTOA_calc_model_lwp(istart:istart+11,:,:),1);       
end

%% Keep cf and Nd constant

%cf = lowcf_dat_ens_mean;
cf = totcf_dat_ens_mean;
cf2 = cf;

mean_cf = meanNoNan(cf(istart:end,:,:),1);
cf_const = repmat(mean_cf,[1 1 size(cf,1)]);
cf_const = permute(cf_const,[3 1 2]);
cf = cf_const;


cfmin_lwp = cf_min;
%cfmin_lwp = 0.05;
cf2(cf2<cfmin_lwp)=cfmin_lwp;

lwp = lwp_dat_ens_mean ./ cf2; %convert to in-cloud LWP
%lwp = lwp_dat_ens_mean

nd = Nd_dat_ens_mean;
mean_val = meanNoNan(nd(istart:end,:,:),1);
val_const = repmat(mean_val,[1 1 size(nd,1)]);
val_const = permute(val_const,[3 1 2]);
nd = val_const;

[Ac_calc_model_Nd_cf,tau_calc_model_Nd_cf,A_calc_model_Nd_cf,SWTOA_calc_model_Nd_cf] = calc_SW(cf,lwp*1e3,nd,sw,A_clear,trans,cf_min,0,NaN,i_Liu); 
SWTOA_calc_model_mean = meanNoNan(SWTOA_calc_model,1);

clear SWTOA_calc_model_annual_Nd_cf
for iy=1:length(years_sw_calc)
   istart=(iy-1)*12+1;
   SWTOA_calc_model_annual_Nd_cf(iy,:,:) = meanNoNan(SWTOA_calc_model_Nd_cf(istart:istart+11,:,:),1);       
end




%% plot SW
var_UM = ['SW TOA model ' str_type_calc ' calc'];
dat_modis = SWTOA_calc_model_mean;
subtitle_str = var_UM;
add_str='';
tit_str_clean='Model values';
figure
ioverride_proj_type=1;
proj_type_DRIVER='ortho';
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);
xlabel(hc,'W m^{-2}'); %label the colour bar

%Plot the actual SW TOA up
subtitle_str = ['SW TOA model ACTUAL'];
dat_modis = meanNoNan(SW_up_TOA_dat_ens_mean,1);
figure
ioverride_proj_type=1;
proj_type_DRIVER='ortho';
irestrict_domain_DRIVER=0;
UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
caxis([0 100]);
xlabel(hc,'W m^{-2}'); %label the colour bar

%%
if isave_plot_SWTOA_calcs_draft==1
    savename=[savedir_date titlenam_driver];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    save_map_plot_data; %sets opts values for saving data
    out_file = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
    %        close(gcf);
end

% %save the map data for later plotting in a subplot plot
% dat_file = [out_file '.mat'];
% save(dat_file,'dat_modis','-V7.3');


%% Create dat_ukesm for use the the timeseries plotting script
f_sw_calc=1.65;

%dat_ukesm_save = dat_ukesm;
clear dat_ukesm
dat_ukesm.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm.years_ukesm_1d = years_sw_calc';

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_calcs


%% Set up for plotting timeseries of sensitivities (CF, Nd, etc. constant)
f_sw_calc=1.65;

%dat_ukesm_save = dat_ukesm;

%All vars varying
clear dat_ukesm
dat_ukesm.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm.years_ukesm_1d = years_sw_calc';

dat_ukesm.dat_annual = f_sw_calc * SWTOA_calc_model_annual;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm.dat_annual_ens=a;

%% Nd and CF constant
sw_str = 'CF and Nd constant';
clear dat_ukesm_Nd_cf
dat_ukesm_Nd_cf.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm_Nd_cf.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm_Nd_cf.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm_Nd_cf.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm_Nd_cf.years_ukesm_1d = years_sw_calc';

dat_ukesm_Nd_cf.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm_Nd_cf.dat_annual_ens=a;

ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions


%% Nd constant
sw_str = 'Nd constant';
clear dat_ukesm_Nd
dat_ukesm_Nd.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm_Nd.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm_Nd.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm_Nd.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm_Nd.years_ukesm_1d = years_sw_calc';

dat_ukesm_Nd.dat_annual = f_sw_calc * SWTOA_calc_model_annual_Nd;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual_Nd,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm_Nd.dat_annual_ens=a;

dat_ukesm_sw_sens = dat_ukesm_Nd;
ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions



%% CF constant
sw_str = 'CF constant';
clear dat_ukesm_cf
dat_ukesm_cf.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm_cf.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm_cf.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm_cf.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm_cf.years_ukesm_1d = years_sw_calc';

dat_ukesm_cf.dat_annual = f_sw_calc * SWTOA_calc_model_annual_cf;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual_cf,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm_cf.dat_annual_ens=a;

dat_ukesm_sw_sens = dat_ukesm_cf;
ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions


%% LWP constant
sw_str = 'LWP constant';
clear dat_ukesm_lwp
dat_ukesm_lwp.gcm_Plon2D_UM = dat_ukesem_save.gcm_Plon2D_UM;
dat_ukesm_lwp.gcm_Plat2D_UM = dat_ukesem_save.gcm_Plat2D_UM;
dat_ukesm_lwp.gcm_Plon2D_edges_UM = dat_ukesem_save.gcm_Plon2D_edges_UM;
dat_ukesm_lwp.gcm_Plat2D_edges_UM = dat_ukesem_save.gcm_Plat2D_edges_UM;
dat_ukesm_lwp.years_ukesm_1d = years_sw_calc';

dat_ukesm_lwp.dat_annual = f_sw_calc * SWTOA_calc_model_annual_lwp;
%Just repeat with dummy data for ensemble values for now - could actually
%calculate the individual ensembles.
a = repmat(f_sw_calc * SWTOA_calc_model_annual_lwp,[1 1 1 9]);
a=permute(a,[4 1 2 3]);
dat_ukesm_lwp.dat_annual_ens=a;

dat_ukesm_sw_sens = dat_ukesm_lwp;
ACSIS_Robson_paper_plot_timeseries_RUN_multi_SW_contributions