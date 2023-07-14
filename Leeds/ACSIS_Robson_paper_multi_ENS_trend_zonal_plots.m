ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on
% var_ukesm

iens_mean = 0; %whether to do a sub-plot of all members, or plot the ens mean
iens_single = 0; ens_single=6;
isave_plot = 1;
i_plot_all_boxes=0;
iplot_obs=0;
iplot_ens=0; %plot individual member profiles
iamip=0;

land_ocean = 'land+ocean';
land_ocean = 'ocean only';
land_ocean = 'land only';  

proj_type_DRIVER='ortho'; plot_region='NA'; %Spherical globe projection - "angle of view" is chosen in plot_global_maps at present 
proj_type_DRIVER='other'; plot_region='global'; %Full global map in miller projection

iscreen_sig=0;

iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;

yr_start_trend = 1985; yr_end_trend = 2014;
%yr_start_trend = 1850; yr_end_trend = 1970;
%yr_start_trend = 2003; yr_end_trend = 2014; clims=[-6 6];

savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
    var_ukesm '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];

if iamip==1
    %AMIP data - load from previously saved trend map
    switch var_ukesm
        case 'rsut'
            amip_dat=load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_rsut_1985_to_2014.mat');
        otherwise
            amip_dat=load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1-AMIP_' var_ukesm '_1985_to_2014.mat']);
    end
end

switch var_ukesm
    case 'ts'
        lor='SouthEast';
    otherwise
        lor='NorthEast';
end



if iens_mean==1
    nplots=1;
    iens_inds=1;
elseif iens_single==1
    iens_inds=ens_single;
else
    nplots=size(dat_ukesm.dat_annual_ens,1);
    iens_inds = 1:nplots;
    subplotting=1; %
    xsub=3; ysub=3; %no. rows, columns
end



% Load the trend data for the ensemble
 load(savefile);
 %var_UM = ['UKESM ' var_str ' trend of ensemble mean between y' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and y' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_trend];
 var_UM = ['UKESM ' var_str ' trend of ensemble mean between y' num2str(yr_start_trend) ' and y' num2str(yr_end_trend) '; ' units_str_trend];
 
 
 
clear trend_vs_lat
for iens=iens_inds
        
    dat_zonal_in = squeeze(trend_map{iens}.coeffs(2,:,:));  
    if iens==1    
        clear opts; opts.screen_type_optional = land_ocean;
        %[dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_zonal_in,gcm_Plat2D_UM,gcm_Plon2D_UM);
        [dat_modis,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_zonal_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
        opts.lmask_in_optional = lmask_out;
    else
        %[dat_tmp,lmask_out] = SW_Robson_obs_land_mask_FUNC(dat_zonal_in,gcm_Plat2D_UM,gcm_Plon2D_UM,lmask_out);
        
        [dat_modis,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_zonal_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);        
    end            
    
    trend_vs_lat(:,iens) = meanNoNan(dat_modis,2);            
    
end

[trend_vs_lat_ens_mean,N,trend_vs_lat_ens_std] = meanNoNan(trend_vs_lat,2);
[trend_vs_lat_ens_min,imin_ens] = min(trend_vs_lat,[],2);
[trend_vs_lat_ens_max,imax_ens] = max(trend_vs_lat,[],2);
lats = gcm_Plat2D_UM(:,1);

%% AMIP model
if iamip==1
    dat_zonal_in = squeeze(amip_dat.trend_map_ens_mean.coeffs(2,:,:));
    [dat_modis,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_zonal_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
    trend_vs_lat_amip = meanNoNan(dat_modis,2);
else
    trend_vs_lat_amip = NaN * ones(size(trend_vs_lat_ens_mean));
end

%% obs

if iplot_obs==1
    
    obs_loadname = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Deep-C_trend_map_dat.mat'];
    obs_dat_in = load(obs_loadname);
    
    %For Deep-C SW this is ok
    lat_obs2 = gcm_Plat2D_UM;
    lon_obs2 = gcm_Plon2D_UM;
    
    dat_zonal_in = squeeze(obs_dat_in.trend_dat_map{1}.coeffs(2,:,:));
    
    clear opts; opts.screen_type_optional = land_ocean;
    [dat_modis,lmask_out] = SW_Robson_obs_land_mask_FUNC2(dat_zonal_in,lat_obs2,lon_obs2,opts);
    trend_vs_lat_obs = meanNoNan(dat_modis,2);
    
end

%% plot
%figure('position',scrsz);
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%set(gcf,'position',[5 30 1252 590]);
%set(gcf,'position',[5 30 500 590]);
set(gcf,'position',[5 30 1200 590]);

if iplot_ens==1
    for iens=iens_select %iens_inds
        p=plot(trend_vs_lat(:,iens),lats,'-','linewidth',3,'color',[0.7 0.7 0.7]); hold on
        %This line stops the line being listed in the legend :-
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

clear leg_str; ileg=1;
plot(trend_vs_lat_ens_mean,lats,'b-','linewidth',3); hold on
leg_str{ileg}='UKESM ens mean'; ileg=ileg+1;

plot(trend_vs_lat_ens_min,lats,'b--','linewidth',1);
leg_str{ileg}='UKESM ens min'; ileg=ileg+1;
plot(trend_vs_lat_ens_max,lats,'b--','linewidth',1);
leg_str{ileg}='UKESM ens max'; ileg=ileg+1;
%Plotting +/- 2*std_dev looks similar
%plot(trend_vs_lat_ens_mean-2*trend_vs_lat_ens_std,lats,'b-','linewidth',1);
%plot(trend_vs_lat_ens_mean+2*trend_vs_lat_ens_std,lats,'b-','linewidth',1);

if iamip==1
    plot(trend_vs_lat_amip,lats,'r-','linewidth',3); hold on
    leg_str{ileg}='UKESM AMIP'; ileg=ileg+1;
end

if iplot_obs==1
    plot(trend_vs_lat_obs,lats,'k-','linewidth',3);
    leg_str{ileg}='Deep-C obs'; ileg=ileg+1;
end

increase_font_size_map_figures;   %This creates a gap between map and
%colorbar...
%fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?

xlabel([var_str ' trend (' units_str_trend ')']);
ylabel('Latitude (degrees)');
legend(leg_str,'location',lor);

grid on


tit_short = [var_str ' trend ' num2str(yr_start_trend) ' to ' num2str(yr_end_trend) ' (' land_ocean ')'];
title(tit_short);




        
        

for iens=[] %iens_inds
    if (iens_mean~=1 & iens_single~=1)
        %hs{iens}=subplot(xsub,ysub,iens);
        %isub=isub+1;
        [jsubplot,isubplot]=ind2sub([xsub ysub],iens);
        ygap=0.06; %gap between plots
        ygap=0.025; %gap between plots
        ypos_start = 0.15; %ypos of first plot from bottom (to get enough clearance for colour bar)
        yheader = 0.07; %gap at the top of the plot
        sub_height = 1/ysub * (1 - ypos_start - yheader - (ysub-1)*ygap);               
        ypos = ypos_start + (ysub-isubplot)*(sub_height+ygap);
        
        xgap=0.0; %gap between plots
        xpos_left = 0.23; %xpos of first plot from left
        xpos_right = xpos_left; %gap at the right of the plot
        sub_width = 1/xsub * (1 - xpos_left - xpos_right - (xsub-1)*xgap);               
        xpos = xpos_left + (jsubplot-1)*(sub_width+xgap);
        
        
        %set(hs{iens},'Position',[(jsubplot-1)*1/xsub ypos 0.9/xsub sub_height]);        
        set(hs{iens},'Position',[xpos ypos sub_width sub_height]);
        
        %axes(hs{iens});title(['iens=',num2str(iens)]);
        %axes(hs{iens});title(['']);
    end
end

if isave_plot==1
    savename2 = ['zonal_UKESM1 ' var_str ' trends y' num2str(yr_start_trend) ' and y' num2str(yr_end_trend) '_' land_ocean];
    savename=[savedir_date plot_region '_' savename2];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    [savename_out] = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)

    save([savename_out '_DATA.mat'],'trend_vs_lat','trend_vs_lat_amip','lats','-V7.3');
    if iplot_obs==1
        save([savename_out '_DATA.mat'],'trend_vs_lat_obs','-APPEND');
    end
end
