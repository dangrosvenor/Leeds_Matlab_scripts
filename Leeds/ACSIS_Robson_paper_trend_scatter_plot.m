box_region='4'; %N Atlantic basin
data_to_plot = 'Ens mean';
data_to_plot = 'Individual ens members';
%data_to_plot = 'AMIP run';
            
ACSIS_Robson_paper_choose_regional_box2 %run script - also chooses ylims, etc.

switch data_to_plot
    case 'Ens mean'        
        iens_individual = 0;
    case 'Individual ens members'        
        iens_individual = 1;
    case 'AMIP run'        
        iens_individual = 0;        
end

if iens_individual==0
    nplots_scatter=1;
    subplotting_DRIVER=0;
else
    nplots_scatter=length(trend_map);
    subplotting_DRIVER=1; %
    xsub=3; ysub=3;%no. rows, columns
end


ilat = find(dat_ukesm.gcm_Plat2D_UM(:,1)>LAT_val(1) & dat_ukesm.gcm_Plat2D_UM(:,1)<LAT_val(2));
ilon = find(dat_ukesm.gcm_Plon2D_UM(1,:)>LON_val(1) & dat_ukesm.gcm_Plon2D_UM(1,:)<LON_val(2));



%Make mask to screen out land points
if iscreen_land==1
    load_type = 'merged netCDF';
    var_UM = 'Land_mask'; %Max height of cloud with in-cloud LWC>=0.05 g/kg (in metres); my diag, not COSP
    um_case='u-bf666'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; icoarse=0; ivar_dir=1; %global ACSIS PD emissions
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    dat_global = UM_load_merged_netCDF(dirUM,var_UM,run_type,load_type);
    
    
    land_mask_mean2 = squeeze(ones(size(trend_map{1}.coeffs(2,:,:))));
    land_mask_mean2(dat_global.dat==1)=NaN;
  
    
%     land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
%     land_mask_ens = permute(land_mask_ens,[3 1 2]);
    
else
    land_mask_mean2 = squeeze(ones(size(trend_map{1}.coeffs(2,:,:))));    
    %land_mask_ens = repmat(land_mask_mean,[1 1 size(dat_ukesm.dat_annual_ens,1)]);
    %land_mask_ens = permute(land_mask_ens,[3 1 2]);
    
end

land_mask_mean = NaN*ones(size(land_mask_mean2));
land_mask_mean(ilat,ilon)=land_mask_mean2(ilat,ilon);

%figure('position',scrsz);
%set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%%set(gcf,'position',[5 30 1252 590]);
%set(gcf,'position',[5 30 800 590]);


%%

for iens=1:nplots_scatter
    nsub=iens;
    isub=iens;
    a=xsub;
    b=ysub;
    %subplot positions are controlled in ~/matlab/work/graphs/plotTimeHeightVap3.m
    % open in vim not matlab since very large and slows the editor down.
    
    X_driver = squeeze(trend_dat_map{it_trend}.coeffs(2,:,:)) .* land_mask_mean;
    switch data_to_plot
        case 'Ens mean'
            Y_driver = squeeze(trend_map_ens_mean.coeffs(2,:,:)) .* land_mask_mean; 
            ylabelstr='UKESM1 ens. mean';
        case 'Individual ens members'
            Y_driver = squeeze(trend_map{iens}.coeffs(2,:,:)) .* land_mask_mean; 
            ylabelstr='UKESM1';
        case 'AMIP run'
            savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
            amip_dat = load(savefile);            
            Y_driver = squeeze(amip_dat.trend_map_ens_mean.coeffs(2,:,:)) .* land_mask_mean; 
            ylabelstr='UKESM1-AMIP';
            
            
            
    end
    %Y_driver(inan)=NaN;
    %Z_driver = forcing(inds_PI_clear_low);
    
    xlabelstr = 'Observations';
    
    %Run script to plot
    DRIVER_template_2D_histo_SW_trends
    plot([Xbins(1) Xbins(end)],[Ybins(1) Ybins(end)],'w');    
    plot([0 0],[Ybins(1) Ybins(end)],'w--');
    plot([Xbins(1) Xbins(end)],[0 0],'w--');
    
    title([box_region_str ', r=' num2str(corr_coeffXY,'%1.3f')]);
    
    if iens_individual==1
        t=get(gca,'title');
        tstr=get(t,'String');
        %tstr{2}=[tstr{2} ' iens=' num2str(iens)];
        tstr=[tstr ' iens=' num2str(iens)];
        %title(tstr);
        %title('');
        title(['r=' num2str(corr_coeffXY,'%1.3f')]);
        
        %set(gcf,'position',[5 30 500 590]);
        
        %jsubplot is the x-value (across screen)
        if isubplot ~= xsub
            xlabel('');
            set(gca,'XTickLabel','');
        end
        
    end
    
    %set(gcf,'color','w');
    savename=[savedir_date '2D histogram SW trends ' box_region_str ' region'];
    iplot_box_on_map=0;
    sat_data = X_driver;
    um_data = Y_driver;
    
    var_Latex = 'SWtrend';
    DRIVER_calc_biases_for_regions
    
    NA_biases{iens} = region_biases{1};
    
end


clear opts
%        opts.iplot_png=1;
opts.iplot_eps=1;
saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts);



