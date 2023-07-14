i_calc=1;
iregrid_obs=1; %regrid MODIS data 
iload_obs_regridded=0;

LAT_val_DRIVER_override = [-1e9 1e9]; LON_val_DRIVER_override = [-1e9 1e9];

isave_plot=0;

proj_type_DRIVER='ortho'; plot_region_str='NA'; %Spherical globe projection - "angle of view" is chosen in plot_global_maps at present 
proj_type_DRIVER='other'; plot_region_str='global'; %Full global map in miller projection

season='Annual';


iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;
iscreen_sig=0;

iplot_boxes=0;

gcm_Plat2D_GENERIC = gcm_Plat2D_UM;
gcm_Plon2D_GENERIC = gcm_Plon2D_UM;
[gcm_Plon2D_edges_GENERIC,gcm_Plat2D_edges_GENERIC] = get_edges_lat_lon(gcm_Plon2D_GENERIC,gcm_Plat2D_GENERIC);
daynum_timeseries3_GENERIC = 1;
gcm_time_UTC_GENERIC = 1;
igeneric_plot = 1;

iplot_individual_ens=0;

savedat_name = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/' obs_str '_trend_map_dat_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
savedat_name_regridded = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/' obs_str '_map_dat_regridded_to_UM.mat'];

switch obs_str %set in ACSIS_Robson_paper_load_data.m script
     case 'ESA_CCI'
         yr_start=1985; yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
         %yr_start=1985;
         years_obs_map = years_obs_CCI;
         obs_annual_map = obs_annual_map_CCI;
         clims=[-0.005 0.005];
         
         if iregrid_obs==1
             clear obs_annual_map
             for iy=1:size(obs_annual_map_CCI,1)
                 obs_annual_map(iy,:,:) =  griddata(gcm_Plat2D_CCI,gcm_Plon2D_CCI,squeeze(obs_annual_map_CCI(iy,:,:)),gcm_Plat2D_UM,gcm_Plon2D_UM);
             end
             save(savedat_name_regridded,'-V7.3','obs_annual_map');
         elseif iload_obs_regridded==1
             load(savedat_name_regridded,'obs_annual_map');; %load the regridded data    
         else
             gcm_Plat2D_GENERIC = gcm_Plat2D_CCI;
             gcm_Plon2D_GENERIC = gcm_Plon2D_CCI;
             [gcm_Plon2D_edges_GENERIC,gcm_Plat2D_edges_GENERIC] = get_edges_lat_lon(gcm_Plon2D_GENERIC,gcm_Plat2D_GENERIC);
             daynum_timeseries3_GENERIC = 1;
             gcm_time_UTC_GENERIC = 1;
             
         end
    
    case 'HadISST'
         yr_start=years_obs_ts(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
         yr_start=1985;
         years_obs_map = years_obs_ts;
         obs_annual_map = obs_annual_map_ts;
         clims=[-0.4 0.4];
         
    case 'Deep-C'
         yr_start=years_obs_DeepC(1); yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
         yr_start=1985;
         %yr_start=1995; yr_end=2014;
         %yr_start=1985; yr_end=2009;
         yr_start=1985; yr_end=2001;
         years_obs_map = years_obs_DeepC;
         obs_annual_map = obs_annual_map_DeepC;
         clims=[-0.4 0.4];

         
    case 'MODIS'
         yr_start=2003; yr_end=2014; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
         years_obs_map = years_obs_modis;
         
         %obs_annual_map =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,obs_annual_map_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
         igeneric_plot = 1;
         
         clims=[-6 6];         

         
         if iregrid_obs==1
             clear obs_annual_map
             for iy=1:size(obs_annual_map_modis,1)
                 obs_annual_map(iy,:,:) =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,squeeze(obs_annual_map_modis(iy,:,:)),gcm_Plat2D_UM,gcm_Plon2D_UM);
             end
             save(savedat_name_regridded,'-V7.3','obs_annual_map');
         elseif iload_obs_regridded==1
             load(savedat_name_regridded,'obs_annual_map'); %load the MODIS data
         else
             obs_annual_map = obs_annual_map_modis;
             
             gcm_Plat2D_GENERIC = gcm_Plat2D_AMSRE;
             gcm_Plon2D_GENERIC = gcm_Plon2D_AMSRE;
             [gcm_Plon2D_edges_GENERIC,gcm_Plat2D_edges_GENERIC] = get_edges_lat_lon(gcm_Plon2D_GENERIC,gcm_Plat2D_GENERIC);
             daynum_timeseries3_GENERIC = 1;
             gcm_time_UTC_GENERIC = 1;
         end
         
    case 'PATMOSx'   
        %Data only goes from 1983 to 2009; so choose the same period for
        %the model when comparing.
         yr_start=1983; yr_end=2009; %yr_end=years_obs(end);  %yr_end_trend_box_obs(it_trend);
         yr_start=1985; yr_end=2009;
         years_obs_map = years_obs_PATMOS;
         obs_annual_map = obs_annual_map_PATMOS/100; %Is in % anomaly, so convert to 0-1 form.
         clims=[-2.5 2.5]*1e-3;         
         
    otherwise
            error('Need to set obs_str here');
         
         
end




%% Calculate the OBS trend global map
if i_calc==1
    for it_trend=1:1 %length(yr_start_trend_box_obs)
       
        switch season
            case 'DJF'
                yr_start = yr_start + 1; %start a year later because data starts in Jan 1985
        end
        yr_start_trend_used_box_obs2=yr_start; yr_end_trend_used_box_obs2=yr_end;                
        
        [trend_dat_map{it_trend}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start, yr_end, ...
            it_trend, years_obs_map, obs_annual_map, p_conf);
    end
    
else
    load(savedat_name);  
    it_trend=1;
end


%% Plot
%figure('position',scrsz);
%figure
%set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%set(gcf,'position',[5 30 1252 590]);
%set(gcf,'position',[5 30 500 590]);

%yr_start=yr_start_trend_used_box_obs2; yr_end=yr_end_trend_used_box_obs2;
istart=find(dat_ukesm.years_ukesm_1d==yr_start);
iend=find(dat_ukesm.years_ukesm_1d==yr_end);

dat_modis = squeeze(trend_dat_map{it_trend}.coeffs(2,:,:));

switch obs_str
    case 'MODIS'
        dat_modis =  griddata(gcm_Plat2D_AMSRE,gcm_Plon2D_AMSRE,dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM);
end

    
 %var_UM = [obs_str ' ' var_str ' trend of ensemble mean between y' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and y' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
 var_UM = [obs_str ' ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(yr_end)];
        %title(tit_short);
    %tit_str_clean = ['UKESM, iens=' num2str(iens) ' ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];
    %subtitle_str = tit_str_clean;
    %add_str = [' ' units_str_trend];
    %subtitle_str=num2str(iens); add_str='';
    subtitle_str=var_UM; add_str='';   
    

%run plotting script
    %figure
    figure
    set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
    %set(gcf,'position',[5 30 1252 590]);
    set(gcf,'position',[5 30 500 590]);
    ioverride_proj_type=1;
    %proj_type_DRIVER='ortho'; %set below now
    irestrict_domain_DRIVER=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims);
    %caxis([-0.3 0.3]);
     %increase_font_size_map_figures;   %This creates a gap between map and
        %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?    
    xlabel(hc,units_str_trend); %label the colour bar
    title(var_UM);
    
     if iscreen_sig==1
        add_str=' screened for significance';
%         switch region_choice
%             case 'global'
%             otherwise
                
               % m_plot(gcm_Plon2D_AMSRE(itrend_not_sig_MODIS),gcm_Plat2D_AMSRE(itrend_not_sig_MODIS),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
               
                %2-tailed t-test with auto-correlation effect added in.
                m_plot(gcm_Plon2D_GENERIC(trend_dat_map{1}.itrend_not_sig2),gcm_Plat2D_GENERIC(trend_dat_map{1}.itrend_not_sig2),'ko','markersize',marker_size,'markerfacecolor','k'); %
                %2-tailed t-test with no auto-corr
                m_plot(gcm_Plon2D_GENERIC(trend_dat_map{1}.itrend_not_sig),gcm_Plat2D_GENERIC(trend_dat_map{1}.itrend_not_sig),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        %end
        
     end
        
    % Plot the box
    ioverride_box_colour=1;
    irotated_pole_box=0;
    itext_in_box=0;
    imap=1;
    col_str='k-';
    box_lwidth = 3;
    if iplot_boxes==1
        ACSIS_Robson_paper_plot_all_boxes
    end
    
    if i_calc==1
       %save(savedat_name,'-V7.3','trend_dat_map');        
    end
    
    if isave_plot==1        
        savename=[savedir_date plot_region_str '_' var_UM];
        clear opts
        %        opts.iplot_png=1;
        opts.iplot_eps=1;
        [savename_out] = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);
        
    end

