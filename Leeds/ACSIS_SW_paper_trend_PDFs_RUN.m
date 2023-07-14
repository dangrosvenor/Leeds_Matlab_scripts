regions={'1','3','8','4'};
regions={'4'};
%regions={'10'};
%regions={'11'}; %latest US outflow region
%regions={'12'}; %US mainland and east coast for emissions.
%regions={'13'}; %NA region up to 50S instead of 60S to avoid sea-ice region.
%regions={'0'}; %Global
%regions={'00'}; %Global -60 to +60
%regions={'01'}; %Global -55 to 60
%regions={'02'}; %Global -50 to 55
%regions={'03'}; %Southern Ocean -60 to -40
%regions={'14'}; %Eastern NA
%regions={'15'}; %Southern Atlantic

land_ocean = 'land+ocean';
land_ocean = 'ocean only';
%land_ocean = 'land only'; 

season='Annual';

iplot_amip=1;
isave_pdf=0;
iplot_obs=0;

for ibox=1:length(regions)
    
    box_region = regions{ibox};
    
    % ---
    ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.
    % ---
    
    LAT_val_DRIVER = LAT_val; LON_val_DRIVER = LON_val; %set in ACSIS_Robson_paper_choose_regional_box2.m
    LAT_UM = gcm_Plat2D_UM(:,1);
    LON_UM = gcm_Plon2D_UM(1,:)';
    
    yr_start_trend=1985; yr_end_trend=2014;
    ylabelstr=[var_str ' trend (' units_str_trend ')'];
    
    switch var_ukesm
        case {'SW_up_TOA','rsut'}
            % Define the bins
            %Ybins_DRIVER = [-0.01 10.^[log10(20):0.15:log10(2500)]]; ichoose_Ybins=1;
            Ybins_DRIVER = [-1:0.02:1];
        case {'ts'}
            Ybins_DRIVER = [-0.2:0.005:0.2];
            
    end
            
    switch var_ukesm
        case 'SW_up_TOA'
            savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                var_ukesm '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
        otherwise
            savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1_' ...
                var_ukesm '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
    end
    map_dat_in = load(savefile);
    Nens = length(map_dat_in.trend_map);
            
    %AMIP data - load from previously saved trend map
    switch var_ukesm
        case 'rsut'
            amip_dat=load('/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_rsut_1985_to_2014.mat');
        otherwise
            amip_dat=load(['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_UKESM1-AMIP_' var_ukesm '_1985_to_2014.mat']);
    end
    
    if iplot_obs==1
        %Obs data (Deep-C)
        obs_loadname = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/Deep-C_trend_map_dat.mat'];
        obs_dat_in = load(obs_loadname);
            
        %Obs data 
        Y_driver_obs = obs_dat_in.trend_dat_map{1}.coeffs(2,:,:);
    end                
    
    %set ens mean data
    Y_driver_ens_mean = map_dat_in.trend_map_ens_mean.coeffs(2,:,:);
    
    %AMIP data - load from previously saved trend map
    Y_driver_amip = amip_dat.trend_map_ens_mean.coeffs(2,:,:);

        
    
    clear x_save_trend_pdf y_save_trend_pdf h Y_mean_overall_ens
    for iens=1:Nens               
        
        %individual ens member data
        Y_driver = map_dat_in.trend_map{iens}.coeffs(2,:,:);
                
        %Run the script
        ACSIS_SW_paper_trend_PDFs
        close(gcf)
        
        x_save_trend_pdf{iens} = xdat(1).x;
        y_save_trend_pdf{iens} = ydat(1).y;
        Y_mean_overall_ens{iens} = Y_mean_overall;
        
    end
    
    
    %Run the PDF script for the ens mean
    Y_driver = Y_driver_ens_mean;
    ACSIS_SW_paper_trend_PDFs
    close(gcf);
    x_save_trend_pdf_ens_mean = xdat(1).x;
    y_save_trend_pdf_ens_mean = ydat(1).y;
    Y_mean_overall_ens_mean = Y_mean_overall;
    
    %Run the script for the AMIP model    
    Y_driver = Y_driver_amip;
    ACSIS_SW_paper_trend_PDFs
    close(gcf)
    x_save_trend_pdf_amip = xdat(1).x;
    y_save_trend_pdf_amip = ydat(1).y;
    Y_mean_overall_amip = Y_mean_overall;
    
    if iplot_obs==1
        %Run the PDF script for the obs
        Y_driver = Y_driver_obs;
        ACSIS_SW_paper_trend_PDFs
        close(gcf)
        x_save_trend_pdf_obs = xdat(1).x;
        y_save_trend_pdf_obs = ydat(1).y;
        Y_mean_overall_obs = Y_mean_overall;
    end
    
end

%% Plot final PDF
figure('color','w');
set(gcf,'position',[3         197        856         422]);
ileg=1; clear leg_str
clear hline

ind_ens=0;
for iens=1:Nens
    if ind_ens==1
        figure('color','w');
        set(gcf,'position',[3         197        856         422]);
        ileg=1; clear leg_str
        clear hline
    end
    
    hline{iens}=plot(x_save_trend_pdf{iens}',y_save_trend_pdf{iens}');
    
    if iens>1 & ind_ens==0
        %This line stops the line being listed in the legend :-
        set(get(get(hline{iens},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    elseif ind_ens==1
        leg_str{ileg}=['iens=' num2str(iens)]; ileg=ileg+1;
    else
        leg_str{ileg}='Ensemble members'; ileg=ileg+1;
    end
    
    set(hline{iens},'color',[0.5 0.5 0.5]);
    hold on
    
    if ind_ens==1 | iens==Nens
        
        hline2=plot(x_save_trend_pdf_ens_mean',y_save_trend_pdf_ens_mean','linewidth',3);
        set(hline2,'color',[0 0 1]);
        leg_str{ileg}='Ensemble mean'; ileg=ileg+1;
        
        if iplot_amip==1
            hline3=plot(x_save_trend_pdf_amip',y_save_trend_pdf_amip','linewidth',3);
            set(hline3,'color',[1 0 0]);
            leg_str{ileg}='AMIP model'; ileg=ileg+1;
        end
        
        if iplot_obs==1
            hline4=plot(x_save_trend_pdf_obs',y_save_trend_pdf_obs','linewidth',3);
            set(hline4,'color',[0 0 0]);
            leg_str{ileg}='Observations'; ileg=ileg+1;
        end
        
        set(gca,'fontsize',18);
        
        grid on
        xlabel(ylabelstr);
        ylabel('No. datapoints');
        tit_str = [box_region_str ' ' land_ocean]
        title(tit_str);
        legend(leg_str);
        
    end
    
end



if isave_pdf==1
    savename=[savedir_date 'PDF ' var_ukesm ' ' tit_str];
    %savename=[savedir_date titlenam ' ' ens_str '_1850_start'];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    
    savename_out = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
end




