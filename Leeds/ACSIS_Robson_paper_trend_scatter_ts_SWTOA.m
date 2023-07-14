isave_pdf=1;

yr_start_trend=1985;
yr_end_trend=2014;

data_to_plot = 'Ens mean';
data_to_plot = 'Individual ens members';
%data_to_plot = 'AMIP run';
%data_to_plot = 'observations';

corr_vars='ts,ts_obs';
corr_vars='ts,sw';

land_ocean = 'land+ocean';
land_ocean = 'ocean only';
%land_ocean = 'land only';

regions={'1','3','8','4'};
regions={'4'};
%regions={'10'};
%regions={'11'}; %latest US outflow region
%regions={'12'}; %US mainland and east coast for emissions.
%regions={'13'}; %NA region up to 50S instead of 60S to avoid sea-ice region.
%regions={'0'}; %Global
regions={'00'}; %Global -60 to 60
%regions={'01'}; %Global -55 to 60
%regions={'02'}; %Global -50 to 55
%regions={'03'}; %Southern Ocean -60 to -40
%regions={'14'};

box_region=regions{1};

iplot_individual_ens=1; %needed for script below
ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.

PDF_script_name = 'DRIVER_template_2D_histo_SW_ts_trends';
xlabelstr = 'ts trend (K yr^{-1})';
ylabelstr = 'SWTOA trend (W m^{-2} yr^{-1})';
plot_name = '2D histogram SW trends';

ilim_x=0; ilim_y=0; ilim_c=0;
iplot_fit = 0; %default - might be overwritten below

iobs=0;
switch data_to_plot
    case 'Ens mean'
        iens_individual = 0;
        dat_str = 'ens_mean';
    case 'Individual ens members'
        iens_individual = 1;  
        dat_str = 'ens';
        clear x_mean_ens y_mean_ens xpdf_ens ypdf_ens mid_Xbins_ens mid_Ybins_ens   
        
        switch corr_vars
            case 'ts,sw'                
                var_01 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                    'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);
                var_02 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                    'ts_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);
                PDF_script_name = 'DRIVER_template_2D_histo_SW_ts_trends';
            case 'sw,Nd'
                var_01 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                    'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);
                var_02 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                    'scldncl_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);
            case 'ts,ts_obs'
                %keep var01 as the ensemble
                var_02 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1-AMIP_' ...
                    'ts_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);                
                var_01 = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                    'ts_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat']);
                iobs=1;
                PDF_script_name = 'DRIVER_template_2D_histo_ts_ts_trends';                
                iplot_fit = 1;
        end
        
    case 'AMIP run'
        iens_individual = 0;
        dat_str = 'AMIP';
    case 'observations'        
        iens_individual = 0;
        dat_str='obs';
end

if iens_individual==0
    nplots_scatter=1;
    subplotting_DRIVER=0;
    xsub=1; ysub=1;
else
    nplots_scatter=length(var_01.trend_map);
    subplotting_DRIVER=1; %
    xsub=3; ysub=3;%no. rows, columns     
end

 max_nsub = xsub*ysub;

ilat = find(dat_ukesm.gcm_Plat2D_UM(:,1)>LAT_val(1) & dat_ukesm.gcm_Plat2D_UM(:,1)<LAT_val(2));
ilon = find(dat_ukesm.gcm_Plon2D_UM(1,:)>LON_val(1) & dat_ukesm.gcm_Plon2D_UM(1,:)<LON_val(2));






%%

clear Pfit   
nfig=0;
for iens=1:nplots_scatter
    nsub = rem(iens,max_nsub);
    if nsub==0
        nsub=max_nsub;
    end    
    if nsub==1
        nfig=nfig+1; %this is when a new fig is made        
    end    
    
        
    %nsub=iens;
    isub=iens;
    a=xsub;
    b=ysub;
    %subplot positions are controlled in ~/matlab/work/graphs/plotTimeHeightVap3.m
    % open in vim not matlab since very large and slows the editor down.
    
    switch data_to_plot
        case 'observations'
            
            %% Observed SWTOA trends
            %Obs data (Deep-C)
            obs_loadname = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Deep-C_trend_map_dat.mat'];
            obs_dat_in = load(obs_loadname);
            Y_driver_in = squeeze(obs_dat_in.trend_dat_map{1}.coeffs(2,:,:));
            clear opts; opts.screen_type_optional = land_ocean;
            [Y_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(Y_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
            %% Observed ts trends (using AMIP here)
            amip_dat=load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1-AMIP_ts_1985_to_2014.mat']);
            X_driver_in = squeeze(amip_dat.trend_map{1}.coeffs(2,:,:)); % .* land_mask_mean;
            clear opts; opts.screen_type_optional = land_ocean;
            [X_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(X_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
    case 'AMIP run'
            
            %% AMIP SWTOA trends    
            rsut_amip = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_rsut_1985_to_2014.mat']);                        
            Y_driver_in = squeeze(rsut_amip.trend_map{1}.coeffs(2,:,:));
            clear opts; opts.screen_type_optional = land_ocean;
            [Y_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(Y_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
            %% Observed ts trends (using AMIP here)
            ts_amip = load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1-AMIP_ts_1985_to_2014.mat']);
            X_driver_in = squeeze(ts_amip.trend_map{1}.coeffs(2,:,:)); % .* land_mask_mean;
            clear opts; opts.screen_type_optional = land_ocean;
            [X_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(X_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
                    
            %         case 'Ens mean'
            %             Y_driver = squeeze(trend_map_ens_mean.coeffs(2,:,:)) .* land_mask_mean;
            %             ylabelstr='UKESM1 ens. mean';
        case 'Individual ens members'           
            if iobs==1
                X_driver_in = squeeze(var_02.trend_map{1}.coeffs(2,:,:));
            else
                X_driver_in = squeeze(var_02.trend_map{iens}.coeffs(2,:,:));
            end            
            clear opts; opts.screen_type_optional = land_ocean;
            [X_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(X_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
            Y_driver_in = squeeze(var_01.trend_map{iens}.coeffs(2,:,:));
            clear opts; opts.screen_type_optional = land_ocean;
            [Y_driver_in,lmask_out] = SW_Robson_obs_land_mask_FUNC2(Y_driver_in,gcm_Plat2D_UM,gcm_Plon2D_UM,opts);
            
            
            %             Y_driver = squeeze(trend_map{iens}.coeffs(2,:,:)) .* land_mask_mean;
            %             ylabelstr='UKESM1';
            %         case 'AMIP run'
            %             savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
            %                 'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
            %             amip_dat = load(savefile);
            %             Y_driver = squeeze(amip_dat.trend_map_ens_mean.coeffs(2,:,:)) .* land_mask_mean;
            %             ylabelstr='UKESM1-AMIP';
            %
            %
            %
    end
    

    
    X_driver = NaN*ones(size(X_driver_in));
    Y_driver = NaN*ones(size(Y_driver_in));
    X_driver(ilat,ilon) = X_driver_in(ilat,ilon);
    Y_driver(ilat,ilon) = Y_driver_in(ilat,ilon);
    
    %Run script to plot
    eval(PDF_script_name);    
    
    %plot one-to-one line
    plot([Xbins(1) Xbins(end)],[Ybins(1) Ybins(end)],'w');
    plot([0 0],[Ybins(1) Ybins(end)],'w--');
    plot([Xbins(1) Xbins(end)],[0 0],'w--');
    
    if iplot_fit==1
        %Do linear fit to data
        Pfit{iens}=polyfit(X,Y,1);
        f=polyval(Pfit{iens},Xbins);
        plot(Xbins,f,'w-.');
    end
    
    %Save the y_mean, x_mean data and the PDFs for y
    eval(['x_mean_' dat_str '{iens} = X_mean_accurate;']);
    eval(['y_mean_' dat_str '{iens} = Y_mean_accurate;']);
    eval(['xpdf_' dat_str '{iens} = NpY;']);
    eval(['ypdf_' dat_str '{iens} = NpX;']);
    eval(['mid_Xbins_' dat_str '{iens} = mid_Xbins;']);
    eval(['mid_Ybins_' dat_str '{iens} = mid_Ybins;']);
    eval(['pdf2D_' dat_str '{iens} = qh;']);
    
    
    tstr = ['r=' num2str(corr_coeffXY,'%1.3f') ', iens=' num2str(iens)];
    
    
    
    if subplotting_DRIVER==1
        %t=get(gca,'title');
        %tstr=get(t,'String');
        %tstr{2}=[tstr{2} ' iens=' num2str(iens)];
        %tstr=[tstr ' iens=' num2str(iens)];
        title(tstr);
        %title('');
        %title(['r=' num2str(corr_coeffXY,'%1.3f')]);
        
        %set(gcf,'position',[5 30 500 590]);
        
        %jsubplot is the x-value (across screen)
        if isubplot ~= xsub
            xlabel('');
            set(gca,'XTickLabel','');
        end
        
    else
        title({[box_region_str ', ' tstr ', ' data_to_plot],[],[box_region_str ', ' data_to_plot]});        
    end   
    
    if ilim_x==1
       set(gca,'xlim',xlims_pdf); 
    end
    if ilim_y==1
       set(gca,'ylim',ylims_pdf); 
    end
    if ilim_c==1
       caxis(clims_pdf); 
    end
    
    %set(gcf,'color','w');
    savename=[savedir_date plot_name ' ' box_region_str ' region ' num2str(nfig)];
    if isave_pdf==1      
        if iens==max_nsub | iens==nplots_scatter
            clear opts
            %        opts.iplot_png=1;
            opts.iplot_eps=1;
            
            savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
        end
    end
    iplot_box_on_map=0;
    sat_data = X_driver;
    um_data = Y_driver;
    
    var_Latex = 'SWtrend';
    DRIVER_calc_biases_for_regions
    
    switch data_to_plot
        case 'Individual ens members'            
            NA_biases{iens} = region_biases{1};
            corr_vals_ens(iens) = corr_coeffXY;
        case 'observations'
            corr_val_obs = corr_coeffXY;
        case 'AMIP run'
            corr_val_amip = corr_coeffXY;
            
    end
    
end


clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
%saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);

ACSIS_Robson_paper_ts_sw_CORR_Box_Whiskier_generic


