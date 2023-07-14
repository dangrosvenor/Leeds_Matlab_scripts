%% UKESM eval









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



switch trend_plot_type
    
    case 'Model least squares trend'

   
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
    
    
    case 'Model 5 years start and end'
    
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
    

    case 'MODIS trend'    
    
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
        ACSIS_Robson_paper_plot_all_boxes
    else
        plot_box_on_map
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
    
end
   