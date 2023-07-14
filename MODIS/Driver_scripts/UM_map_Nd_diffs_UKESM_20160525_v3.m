% Compare the UM Nd data to MODIS data

clear UKESM  %will use this to save the outputs.

% UKESM data
UKESM.Nd_UKESM_savefile = '/home/disk/eos8/d.grosvenor/UM/UKESM/Nd_UKESM_seasonal_clim_2000-2007.mat';
UM=load(UKESM.Nd_UKESM_savefile);
%N.B. info on UM jobs, etc. in UM_map
% Lat and lon are in UM_grid         

% MODIS data
UKESM.Nd_file_UKESM = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2014_allCF_CTT273_region_mask_20150416T063530.mat';
MOD=load(UKESM.Nd_file_UKESM);

%MODIS names are e.g. Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_DJF

UKESM.MOD_str = 'Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_';

is=1;
season_str{is}='DJF'; is=is+1;
season_str{is}='MAM'; is=is+1;
season_str{is}='JJA'; is=is+1;
season_str{is}='SON'; is=is+1;

for iseason=1:4
    
    %Interpolate the UM data onto the MODIS grid
    for imodel=1:3
       UKESM.Nd_UKESM_grid{imodel,iseason} = griddata(UM.UM_grid.gcm_Plon2D_UM, UM.UM_grid.gcm_Plat2D_UM, UM.Nd_UKESM{imodel,iseason}, MOD.MLON',MOD.MLAT );                
       mod_dat = eval(['MOD.' UKESM.MOD_str season_str{iseason}]);
       UKESM.Nd_UKESM_bias_abs{imodel,iseason} = Nd_UKESM_grid{imodel,iseason} - mod_dat;
       UKESM.Nd_UKESM_bias_prc{imodel,iseason} = 100*(Nd_UKESM_grid{imodel,iseason} - mod_dat) ./ mod_dat;    
       
       UKESM.Nd_UKESM_bias_abs_mean{imodel,iseason} = meanNoNan(Nd_UKESM_bias_abs{imodel,iseason}(:),1);
       UKESM.Nd_UKESM_bias_abs_mean_prc{imodel,iseason} = 100*Nd_UKESM_bias_abs_mean{imodel,iseason} ./ meanNoNan(mod_dat(:),1);
       UKESM.Nd_UKESM_bias_prc_mean{imodel,iseason} = meanNoNan(Nd_UKESM_bias_prc{imodel,iseason}(:),1); 
       
       % Overall spatial RMS value from the absolute values
       UKESM.Nd_UKESM_bias_abs_RMS{imodel,iseason} = sqrt( meanNoNan(Nd_UKESM_bias_abs{imodel,iseason}(:) .^2,1) );
       % Overall spatial RMS values from the percentage biases
       UKESM.Nd_UKESM_bias_prc_RMS{imodel,iseason} = sqrt( meanNoNan(Nd_UKESM_bias_prc{imodel,iseason}(:) .^2,1) );
    end
    
    
end




%% Options for PDF plotting

Ybins_DRIVER = [0:20:800];
logbin_norm_driver = 0; %Whether hav lognormal bins and normalise by these
i_plot_norm_driver = 1; %whether to normalise
i_div_bin_widths_driver = 1; %Whether to divide by the bin widths (also normalises)

pdf_type_driver='normal';  %normal or cumulative PDF
%pdf_type_driver='cumulative';



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
%UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/Iceland_Anja/';
UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/UKESM/';

clear fileUM xdat_import ydat_import
idat=1;
%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkh) 100cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkj) 400cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;idat=idat+1;
%fileUM{idat} = 'xkmph_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkmph)'; pole_lat=70; pole_lon=284; idat=idat+1;

% - 12th Nov case, as of analysis started July 2015

%fileUM{idat} = '/xlhg-u/xlhgu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
%fileUM{idat} = '/xlyd-x/xlydx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Nd_fixed_act-10';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlyd-y/xlydy_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-new';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlyd-p/xlydp_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-p/xlydp_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-z/xlydz_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-SWNd';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-z/xlydz_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0];
%marker_styleUM(idat).m='^'; idat=idat+1;

%UM_map.fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'CASIM-Ndvar(increased)';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  UM_map.flag{idat} = 'load_mat';  UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;

%UM_map.fileUM{idat} = '/u-ad234/u-ad234_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad234'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
%UM_map.fileUM{idat} = '/u-ad571/u-ad571_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad571, low aerosol, volcano ON'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
%UM_map.fileUM{idat} = '/u-ad572/u-ad572_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad572, low aerosol, volcano OFF'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
UM_map.fileUM{idat} = '/u-ab642/u-ab642_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA7-N96 AMIP'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;
UM_map.fileUM{idat} = '/u-ab754/u-ab754_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA6-N96 AMIP'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;
UM_map.fileUM{idat} = '/u-ac043/u-ac043_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA7-N96 ORCA1.0'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;

    
    
UM_map.iset_min_clim=1;
UM_map.clim_min=0;
UM_map.iset_max_clim=1;
UM_map.clim_max=650;


UM_map.iplot_markers=0;
    
% Plot the map(s)
if iplot_maps==1      
    vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);
end

%% Plot diff if requested

if plot_diff==1 & length(vars_map_out.P_save)>1
    
    dat_modis = vars_map_out.P_save{1} - vars_map_out.P_save{2};


    %--- run the file to set up the defaults
    plot_global_maps_defaults

    irestrict_domain=UM_map.irestrict_domain_DRIVER; %whether to restrict the domain or not

    thresh_LAT = UM_map.LAT_val_DRIVER;
    thresh_LON = UM_map.LON_val_DRIVER;

    %--- set some options for these particular plot loops
    set_screening = {'none'};
    modis_data_plot = 'Map of 2D data from outside driver script';

    %         iset_min_clim=1;
    %         clim_min=0;
    %         iset_max_clim=1;
    %         clim_max=300;
    %
    %         isave_plot=0;
    %         iplot_markers=0;



    %Set various things
    %        time = time_matlab(it_driver);
    %Round to the nearest minute as sometimes get 18:59:59
    %        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM');
    time_str = vars_map_out.time_str;
    %        if iz>-1
    %            z_str = ['at ' num2str(z_um(iz)/1000,'%.2f') ' km '];
    %        else
    %            z_str ='';
    %        end
    z_str = vars_map_out.z_str;

    titlenam_driver = ['SW difference for ' time_str ' for ' UM_map.labs_UM(1).l ' minus ' UM_map.labs_UM(2).l ' at ' z_str];
    units_str_plot = UM_map.var_units_str;

    mod_data_type='AMSRE';
    gcm_str_select='UM';
    %        daynum_timeseries3_UM = [1:length(time)];
    %        gcm_time_UTC_UM = [1:length(time)];

    %        gcm_Plon2D_UM = lons2D_driver;
    %        gcm_Plat2D_UM = lats2D_driver;

    %        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(lats2D_driver,lons2D_driver);



    %        i_dpcolor=1;
    ifull_swath=0;
    igcm_screen=0;




    %--- Apply override flags
    ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
    % iocean_only=1;
    ioverride_time_selection=0; %Override the times to include
    ioverride_plotglobal_loc=1; %Override the location of the plot window
    ioverride_years_time_screen=0; %Override years for screening?

    %---  Run plot script and save
    plot_global_maps
    if UM_map.isave_plot==1
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
        close(gcf);
    end

    
    
end







%%% PDF plot

if plot_pdf_diff==1 & length(vars_map_out.P_save)>1
    
    for idat_driver=1:2
    
          %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
        
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';
        
        logflag=0;
        dlogflag=0;
        

                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
%        y_axis_vals = 'UM LWP';
         y_axis_vals = 'General y-axis no ilat simple'; %'General GCM-style';
         datatype = 'gcm_data';
         
        ylabelstr = ['SW TOA W m^{-2}'];
        Ybins = Ybins_DRIVER; ichoose_Ybins=1;    
        Y_driver = vars_map_out.P_save{idat_driver};
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
        axis1D = 'y';                                
        
        logbin_norm = logbin_norm_driver;
        i_plot_norm=i_plot_norm_driver;
        i_div_bin_widths=i_div_bin_widths_driver;
        pdf_type = pdf_type_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3
        close(gcf);
        waterVapourMay2005
        close(gcf);
        
        %store the PDF data
         switch pdf_type_driver
            case 'normal'
                ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
                xdat_import(idat_driver) = xdat_norm;
            case 'cumulative'
                ydat_import(idat_driver) = ydat_cum;
                xdat_import(idat_driver) = xdat_cum;
         end
        
       
        labs_import(idat_driver).l = UM_map.labs_UM(idat_driver).l;
        xlab_import = xlab;
        ylab_import = ylab;



    end
    
    
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'SW PDFs';
xlab='SW TOA (W m^{-2})';
ylab = ylab_import;
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=0;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square
ichoose_styles=1;

line_pattern = line_pattern_DRIVER;  line_colour=line_colour_DRIVER; marker_style=marker_style_DRIVER;


%---  Main script to do plots and save
DRIVER_lineplot_watervap

switch pdf_type_driver
    case 'cumulative'
        set(gca,'xscale','log');
%        set(gca,'xlim',[10 800]);
    case 'normal'
%        set(gca,'xscale','log');
%        set(gca,'xlim',[9 800]);
end
        
        
        
        
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
   
     









                            

                            
                            

                            
                
                    


end
  
%% UKESM .mat file contents :-
%
% whos('-file',Nd_UKESM_savefile)
%   Name              Size              Bytes  Class     Attributes
%   Nd_UKESM          3x4             2655552  cell                
%   Ndatap_UKESM      3x4             2655552  cell                
%   UM_grid           1x1              890832  struct  
%   UM_map            1x1                5346  struct              

    
    
%% MODIS data file contents :-
% 
% whos('-file',Nd_file_UKESM)
%   Name                                                         Size               Bytes  Class     Attributes
%   Date_text                                                    1x67                 134  char                
%   MLAT                                                         1x180               1440  double              
%   MLON                                                         1x360               2880  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_ANNUAL              180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_DJF                 180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_JJA                 180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_MAM                 180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_Ndays_ANNUAL        180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_Ndays_DJF           180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_Ndays_JJA           180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_Ndays_MAM           180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_Ndays_SON           180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_SON                 180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_std_dev_ANNUAL      180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_std_dev_DJF         180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_std_dev_JJA         180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_std_dev_MAM         180x360             518400  double              
%   Nd_SZA_lt_70_allCF_AND_CTT_gt_273_clim_std_dev_SON         180x360             518400  double              
%   Notes                                                        1x1571              3142  char                
%   Notes2                                                       1x691               1382  char                
