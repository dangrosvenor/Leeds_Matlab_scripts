%Convert this to a function, so that can pass in different dat_ukesm
%variables (dat_ukesm_nudged, etc.)

ACSIS_Robson_paper_choose_clims_etc %run script to choose clims, units, etc. based on
% var_ukesm

% ioverride_clims = 1;
% if ioverride_clims==1
%    clims = [-0.2 0.2]; 
% end

i_calc = 1; %whether to re-calculate the trend maps, or load from the .mat file
iens_mean = 1; %whether to do a sub-plot of all members, or plot the ens mean
iens_single = 0; ens_single=1;
if i_calc==1
    MIP_maps=MIP; %to prevent accidental overwrite of data - since it uses the loaded data here.
else
    MIP_maps = 'CMIP';
    %MIP_maps = 'AMIP';
    MIP_maps = 'DAMIP';
    %MIP_maps = 'NUDGED';
    %MIP_maps = 'DAMIP hist-aer SW calculated';
end

p_conf = 95;

isave_plot = 0;
iplot_maps=1; %whether to actually plot, or just calculate the trends
i_plot_all_boxes=1;


            
iplot_delta = 0; %whether to plot the delta (or if set to zero the trend).

proj_type_DRIVER='ortho'; plot_region_str='NA'; %Spherical globe projection - "angle of view" is chosen in plot_global_maps at present 
proj_type_DRIVER='other'; plot_region_str='global'; %Full global map in miller projection

LAT_val_DRIVER_override = [-1e9 1e9]; LON_val_DRIVER_override = [-1e9 1e9];

iscreen_sig=0;

iplot_mgrid_lines_DRIVER=1; %whether to plot the grid lines for maps using m_grid
ioverride_ticks_DRIVER=1;

yr_start_trend = 1985; yr_end_trend = 2014;
%yr_start_trend = 1990; yr_end_trend = 2014;
%yr_start_trend = 1983; yr_end_trend = 2009; %PATMOS period
%yr_start_trend = 1985; yr_end_trend = 2009; %PATMOS and Deep-C overlap period
%yr_start_trend = 1870; yr_end_trend = 1970;
%yr_start_trend = 1980; yr_end_trend = 2005; %To match those in Zhou.
yr_start_trend = 1850; yr_end_trend = 1970;
yr_start_trend = 1850; yr_end_trend = 2005;
%yr_start_trend = 1971; yr_end_trend = 2014;
%yr_start_trend = 2003; yr_end_trend = 2014; clims=[-6 6];
%yr_start_trend = 1985; yr_end_trend = 2001; %temp for u-by844 before the run finished

iadd_DAMIP=0;

var_str = var_ukesm; %pass in
switch MIP_maps
    case 'AMIP'
        model_str_map = 'UKESM1-AMIP';
        switch var_ukesm
            case 'SW_up_TOA'
                %savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                %'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
                var_str = 'rsut';
        end
        iens_mean = 1;
        
    case 'DAMIP'
        iadd_DAMIP=1;
        model_str_map = model_str; %'HadGEM3-GC3.1-DAMIP';
        switch var_ukesm
            case 'SW_up_TOA'
                %savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                %'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
                var_str = 'rsut';
                
             case 'rsutcs'
                var_str = var_ukesm;
                %savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                %var_str '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
                
            case 'SWTOA Calc'
                var_str = 'rsut_calc';
                   
        end
        
    case 'DAMIP hist-aer SW clear-sky'   
        model_str_map = MIP;
        var_str = 'rsutcs'
        iadd_DAMIP=1; %used in ACSIS_Robson_paper_choose_regional_box2 when plotting boxes.
                
     case 'DAMIP hist-aer SW calculated'   
        model_str_map = MIP;
        var_str = 'rsutcs_calc'
        iadd_DAMIP=1; %used in ACSIS_Robson_paper_choose_regional_box2 when plotting boxes.   
        
    case 'CMIP'
        model_str_map = 'UKESM1';
        %savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' model_str_map '_' ...
        %    var_ukesm '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
        
    case 'NUDGED'
        model_str_map = ['NUDGED_' nudged_run_str];
        switch var_ukesm
            case 'SW_up_TOA'
                %savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/ensemble_trend_maps_' ...
                %'rsut_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];
                var_str = 'rsut'; 
        end
        iens_mean = 1;    
        
    otherwise
        model_str_map = MIP;
end

ens_str2 = 'ensemble'; %default
subplotting=0;
if iens_mean==1
    nplots=1;
    iens_inds=1;
    ens_str2 = 'ensemble_mean';
elseif iens_single==1
    iens_inds=ens_single;
else
    nplots=size(dat_ukesm.dat_annual_ens,1);
    %eval_str=['nplots=size(' dat_var '.dat_annual_ens,1);'];
    %eval(eval_str);
    iens_inds = 1:nplots;
    subplotting=1; %
    xsub=3; ysub=3; %no. rows, columns
    xsub=2; ysub=2; %no. rows, columns
end

savefile = ['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/' ens_str2  '_trend_maps_' model_str_map '_' ...
        var_str '_' num2str(yr_start_trend) '_to_' num2str(yr_end_trend) '.mat'];

%figure('position',scrsz);
figure
set(gcf,'color','w'); %set background colour of the figure to white for better plots when screen grabbing.
%set(gcf,'position',[5 30 1252 590]);
%set(gcf,'position',[5 30 500 590]);
set(gcf,'position',[5 30 1200 590]);




%clear coeffs t_trend
for iens=iens_inds    
    
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
    
    if i_calc==1        
        if iens==1
            clear trend_map
            dat = squeeze(dat_ukesm.dat_annual(:,:,:));
            [trend_map_ens_mean] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend, yr_end_trend, ...
                1, dat_ukesm.years_ukesm_1d, dat, p_conf);
            map_dat_in.trend_map_ens_mean = trend_map_ens_mean;
        end
        
        if iens_mean==1
            trend_map{iens} = NaN;
        else
            dat = squeeze(dat_ukesm.dat_annual_ens(iens,:,:,:));
            [trend_map{iens}] = ACSIS_Robson_paper_compute_trends_FUNC(yr_start_trend, yr_end_trend, ...
                1, dat_ukesm.years_ukesm_1d, dat, p_conf);
            %    %[coeffs{iens},t_trend{iens}] = trend_linear_fit_nd_data(x,squeeze(dat_ukesm.dat_annual_ens(iens,istart:iend,:,:)),1);
        end
        map_dat_in.trend_map = trend_map;
    else
        %if ~exist('trend_map')
            map_dat_in = load(savefile);
        %end
    end
    
    if iplot_maps==1
        
        if (iens_mean~=1 & iens_single~=1)
            [jsubplot,isubplot]=ind2sub([xsub ysub],iens);
            hs{iens}=subplot(xsub,ysub,iens);
        else
            isubplot=1; jsubplot=1;
        end
    
    %% Map of linear MODEL trend
    yr_start=yr_start_trend_used; yr_end=yr_end_trend_used;
    istart=find(dat_ukesm.years_ukesm_1d==yr_start);
    iend=find(dat_ukesm.years_ukesm_1d==yr_end);
    nyears_trend = yr_end - yr_start;
    if iens_mean==1
        dat_modis = squeeze(map_dat_in.trend_map_ens_mean.coeffs(2,:,:)); 
        ens_str='';
    else
        dat_modis = squeeze(map_dat_in.trend_map{iens}.coeffs(2,:,:));
        ens_str = ['iens=' num2str(iens)];
    end
    
    if iplot_delta==1
        dat_modis = dat_modis * nyears_trend; %plot the delta rather than the trend 
        clims_plot = clims_delta;
    else
        clims_plot = clims;
    end
    
    if iscreen_sig==1
        switch region_choice
            case 'global'
                %dat_modis(trend_map{iens}.itrend_not_sig)=NaN; %make NaN for now, but can put a dot on, etc.
        end
    end
    
    if iplot_delta==1
        units_str_plot = units_str;        
    else
        units_str_plot = units_str_trend;        
    end
    
    var_UM = [model_str_map ' ' var_str ' trend of ensemble mean between ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' and ' num2str(dat_ukesm.years_ukesm_1d(iend)) '; ' units_str_plot];
    %tit_str_clean = ['UKESM, iens=' num2str(iens) ' ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend))];
    %subtitle_str = tit_str_clean;
    %add_str = [' ' units_str_trend];
    %subtitle_str=num2str(iens); add_str='';
    subtitle_str=''; add_str='';
    
    %run plotting script
    %figure
    ioverride_proj_type=1;
    ioverride_LAT_plots=0;
    %proj_type_DRIVER='ortho'; %set at top of script.
    %proj_type_DRIVER='other';
    irestrict_domain_DRIVER=0;
    igeneric_plot=0;
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands_global
    lb_map = lbmap(16,'brownblue'); colormap(flipdim(lb_map,1));
    caxis(clims_plot);
    %increase_font_size_map_figures;   %This creates a gap between map and
        %colorbar...
    fontsize_figure(gcf,gca,18); %Might not increase fonts of everything? E.g., lat lon labels?
    %caxis([-0.3 0.3]);
   
    xlabel(hc,units_str_plot); %label the colour bar

    if iens_mean ~= 1 & iens_single ~= 1
        
        
        
        %gca = hs{iens}; Hc1s = find_peer_colorbars_of_an_axes(gca);
        Hc1s = find_peer_colorbars_of_an_axes(hs{iens});
        if iens>1
            delete(Hc1s)  %just need one big colorbar
        else
            pos=get(Hc1s,'position');
            pos(2)=0.11;
            pos(3)=1-2*pos(1);
            set(Hc1s,'position',pos);
        end
        
    end
        
    if iens_mean==1 | (isubplot==1 & jsubplot==2) | iens_single==1
        tit_short = [model_str_map ' ' var_str ' trend ' num2str(dat_ukesm.years_ukesm_1d(istart)) ' to ' num2str(dat_ukesm.years_ukesm_1d(iend)) ' ' ens_str];
        title(tit_short);
    else
        title('');
    end
    
    
    
    
    
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
        %plot_box_on_map
    end
    
    if iscreen_sig==1
        add_str=' screened for significance';
        switch region_choice
            case 'global'
            otherwise
                %2-tailed t-test with auto-correlation effect added in.
                m_plot(gcm_Plon2D_UM(trend_map{iens}.itrend_not_sig2),gcm_Plat2D_UM(trend_map{iens}.itrend_not_sig2),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
                %2-tailed t-test with no auto-corr
                m_plot(gcm_Plon2D_UM(trend_map{iens}.itrend_not_sig),gcm_Plat2D_UM(trend_map{iens}.itrend_not_sig),'ko','markersize',marker_size,'markerfacecolor','k'); %m_plot works using lon,lat
        end
        
    end 
    
end
    
end

if iplot_maps==1
for iens=iens_inds
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
end

if isave_plot==1 & iplot_maps==1
    savename=[savedir_date plot_region_str '_' var_UM];
    clear opts
    %        opts.iplot_png=1;
    opts.iplot_eps=1;
    [savename_out] = saveas_ps_fig_emf(gcf,[savename],'',0,0,0,'',[],0,opts)
end

if i_calc==1
   save(savefile,'trend_map_ens_mean','trend_map','yr_start_trend','yr_end_trend','-V7.3');
end