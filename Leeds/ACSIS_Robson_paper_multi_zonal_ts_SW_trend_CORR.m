ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on
% var_ukesm

%Create and save the zonal data using ACSIS_Robson_paper_multi_ENS_trend_zonal_plots.m


isave_plot = 0;

dat_pair = 'sw,ts';            
dat_pair = 'sw,scldncl';

%land_ocean = 'land+ocean';
%land_ocean = 'ocean only';
%land_ocean = 'land only';

land_ocean_multi = {'ocean only','land only'};
patt_str={'b-','r-'};

%figure('position',scrsz);
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%set(gcf,'position',[5 30 1252 590]);
%set(gcf,'position',[5 30 500 590]);
set(gcf,'position',[5 30 1200 590]);

clear leg_str; ileg=1;

for iland=1:length(land_ocean_multi)
    
    land_ocean = land_ocean_multi{iland};
    
    
    yr_start_trend = 1985; yr_end_trend = 2014;
    %yr_start_trend = 1850; yr_end_trend = 1970;
    %yr_start_trend = 2003; yr_end_trend = 2014; clims=[-6 6];
    
    
    yr_start_trend_str = num2str(yr_start_trend);
    yr_end_trend_str = num2str(yr_end_trend);
    load_dir = '/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_20200325T024432/';
    
    switch dat_pair
        case 'sw,ts'
            loadfile_x1 = [load_dir 'global_zonal_UKESM1_SW_TOA_up_trends_y' yr_start_trend_str '_and_y' yr_end_trend_str '_' land_ocean '_DATA.mat'];
            loadfile_x2 = [load_dir 'global_zonal_UKESM1_Surface_Temperature_trends_y' yr_start_trend_str '_and_y' yr_end_trend_str '_' land_ocean '_DATA.mat'];
            xlab = ['Correlation SWTOA trends vs ts trends'];
        case 'sw,scldncl'
            loadfile_x1 = [load_dir 'global_zonal_UKESM1_SW_TOA_up_trends_y' yr_start_trend_str '_and_y' yr_end_trend_str '_' land_ocean '_DATA.mat'];
            loadfile_x2 = [load_dir 'global_zonal_UKESM1_N_d_trends_y' yr_start_trend_str '_and_y' yr_end_trend_str '_' land_ocean '_DATA.mat'];            
            xlab = ['Correlation SWTOA trends vs N_d trends'];
    end
            
            savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                var_ukesm '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
   
    
    %% Calc the correlation corr between ts trend and SWTOA trend across the ensemble
    x1_dat = load(remove_character(loadfile_x1,' ','_'));
    x2_dat = load(remove_character(loadfile_x2,' ','_'));
    
    clear corr_vs_lat
    for ilat=1:length(x1_dat.lats)
        x1 = x1_dat.trend_vs_lat(ilat,:);
        x2 = x2_dat.trend_vs_lat(ilat,:);
        corr_vs_lat(ilat) = corr(x1(:),x2(:));
    end
    
    %% plot
    
    
    plot(corr_vs_lat,lats,patt_str{iland},'linewidth',3); hold on
    leg_str{ileg}=land_ocean; ileg=ileg+1;
    
    %plot(trend_vs_lat_ens_min,lats,'b--','linewidth',1);
    %leg_str{ileg}='UKESM ens min'; ileg=ileg+1;
    %plot(trend_vs_lat_ens_max,lats,'b--','linewidth',1);
    %leg_str{ileg}='UKESM ens max'; ileg=ileg+1;
    
    
    %plot(trend_vs_lat_amip,lats,'r-','linewidth',3); hold on
    %leg_str{ileg}='UKESM AMIP'; ileg=ileg+1;
    
    %if iplot_obs==1
    %    plot(trend_vs_lat_obs,lats,'k-','linewidth',3);
    %    leg_str{ileg}='Deep-C obs'; ileg=ileg+1;
    %end
    
end

increase_font_size_map_figures;   %This creates a gap between map and
%colorbar...
%fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?

xlabel(xlab);
ylabel('Latitude (degrees)');
legend(leg_str,'location',lor);

grid on


%tit_short = ['Correlation SWTOA vs ts trends ' num2str(yr_start_trend) ' to ' num2str(yr_end_trend) ' (' land_ocean ')'];
tit_short = ['Correlation SWTOA vs ts trends ' num2str(yr_start_trend) ' to ' num2str(yr_end_trend)];
title(tit_short);




if isave_plot==1
    savename2 = ['zonal_UKESM1 trend corr y' num2str(yr_start_trend) ' and y' num2str(yr_end_trend)];
    savename=[savedir_date plot_region '_' savename2];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    [savename_out] = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
    
    %     save([savename_out '_DATA.mat'],'trend_vs_lat','trend_vs_lat_amip','lats','-V7.3');
    %     if iplot_obs==1
    %         save([savename_out '_DATA.mat'],'trend_vs_lat_obs','-APPEND');
    %     end
end

