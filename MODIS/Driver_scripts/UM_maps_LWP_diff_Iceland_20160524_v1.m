% LWP map plots for all times for UM
% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];
try %to enable the clearing of override flags
    
if ~exist('ioverride_LWP_diff') | ioverride_LWP_diff==0
    clear UM_map
    %Can use for either LWP or RWP
    UM_map.varname = 'LWP';
    UM_map.varname = 'RWP';
end


icoarsen=0; icoarsen_diff=1; %Whether to coarse grain the maps or diff plots (done within UM_maps*_FUNC.m for maps, below for diffs)
iprc_diff=1; %Plot the percentage LWP difference
thresh_LWP=10; %Threshold LWP (smoothed) for which to calculate % differences
 %Used 50 g/m2 for LWP and 10 for RWP
plot_pdf_diff=0;

UM_map.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
UM_map.thresh_LWP_mask = 20;

UM_map.isave_plot=0;

if ~exist('ioverride_LWP_diff') | ioverride_LWP_diff==0
    clear UM_map
    %Can use for either LWP or RWP
    UM_map.varname = 'LWP';
    UM_map.varname = 'RWP';
end

UM_map.VAR_NAME_STR=''; %The name part in the filename to replace VAR_NAME with
            %OR, if the VAR_NAME is not in the supplied filename then will
            %just use the supplied one - make sure it is correct!
UM_map.VAR_NAME_STR=UM_map.varname;

UM_map.flag2='load_mat'; %replacing the old flag for whether is a mat file or not - set to '' if not.
UM_map.var_units_str = 'g m^{-2}';


UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_map.LAT_val_DRIVER = [-90 90]; UM_map.LON_val_DRIVER =[-180 180]; %Not needed unless restricting domain below
UM_map.irestrict_domain_DRIVER=0;

%Choose a height index for 3D data if needed
UM_map.iz = -1; %=3km for Iceland case (can do z_um = nc{'hybrid_ht'}(:); to check (using previously opened nc file)

%UM_map.tol_mins=1;
%time_single = datenum('13-Nov-2008 19:00'); %UM_map.time_range(1) = time_single - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;
%time_single = datenum('01-Sep-2014 21:00'); %UM_map.time_range(1) = time_single - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;
time_single = datenum('01-Sep-2014 12:00'); %UM_map.time_range(1) = time_single - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;

UM_map.time_range(1) = time_single;




%% For PDF plotting

Ybins_DRIVER = [-0.01:50:1500];
%Ybins_DRIVER = [-0.01 10.^[0.1:0.1:log10(2000)] ]; %RWP
%Ybins_DRIVER = [-0.01:4:50];
logbin_norm_driver = 0; %Whether hav lognormal bins and normalise by these
i_plot_norm_driver = 1; %whether to normalise
i_div_bin_widths_driver = 1; %Whether to divide by the bin widths (also normalises)

pdf_type_driver='normal';  %normal or cumulative PDF
pdf_type_driver='cumulative';




%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/Iceland_Anja/';
clear fileUM xdat_import ydat_import
idat=1;
%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkh) 100cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkj) 400cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;idat=idat+1;
%fileUM{idat} = 'xkmph_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkmph)'; pole_lat=70; pole_lon=284; idat=idat+1;

% - 12th Nov case, as of analysis started July 2015

%fileUM{idat} = '/xlhg-u/xlhgu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
%fileUM{idat} = '/xlyd-x/xlydx_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd_fixed_act-10';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlyd-y/xlydy_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-new';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlyd-p/xlydp_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-p/xlydp_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-z/xlydz_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-SWNd';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-z/xlydz_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0];
%marker_styleUM(idat).m='^'; idat=idat+1;

%UM_map.fileUM{idat} = '/xmmz-u/xmmzu_LWP_.pp.nc.mat'; UM_map.labs_UM(idat).l = 'CASIM-Ndvar(increased)';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-v/xmmzv_LWP_.pp.nc.mat';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0
%    0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-x/xmmzx_LWP_.pp.nc.mat';  UM_map.labs_UM(idat).l ='CASIM-Ndvar-0.025';  UM_map.flag{idat} = 'load_mat';  UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0.8 0]; UM_map.marker_styleUM(idat).m='o'; idat=idat+1;        
%UM_map.fileUM{idat} = '/xmmz-w/xmmzw_LWP_.pp.nc.mat';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-10';  UM_map.flag{idat} = 'load_mat';  UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;
UM_map.fileUM{idat} = '/u-ad571/u-ad571_VAR_NAME_.pp.nc.mat'; UM_map.labs_UM(idat).l = 'u-ad571, low aerosol, volcano ON'; UM_map.pole_lat=25; UM_map.pole_lon=165;
    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
UM_map.fileUM{idat} = '/u-ad572/u-ad572_VAR_NAME_.pp.nc.mat'; UM_map.labs_UM(idat).l = 'u-ad572, low aerosol, volcano OFF'; UM_map.pole_lat=25; UM_map.pole_lon=165;
    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
    
UM_map.iset_min_clim=1;
UM_map.clim_min=0;
UM_map.iset_max_clim=1;
UM_map.clim_max=300;


UM_map.iplot_markers=0;
    
    

%% Plot the map
    
%Boutle_case_12Nov2008_UM_maps_20141126T041216_FUNC(UM_map)   
vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);


gcm_Plat2D_edges_UM = vars_map_out.gcm_Plat2D_edges_UM;
gcm_Plon2D_edges_UM = vars_map_out.gcm_Plon2D_edges_UM;
gcm_Plat2D_UM = vars_map_out.gcm_Plat2D_UM;
gcm_Plon2D_UM = vars_map_out.gcm_Plon2D_UM;



%% Now plot the difference

        dat_modis = vars_map_out.P_save{1} - vars_map_out.P_save{2};

        
        if icoarsen_diff==1
            coarsen_method = '2d smoothing filter'; N = 128; M=N;  %Smoothing
            coarsen_method = 'gridpoint blocks';N = 64; M=N; %Coarse graining to NxM gridpoint blocks (rather than smoothing)            
%            coarsen_method = 'lat lon blocks'; dlat_target=0.25; dlon_target=0.25;
            
            coarse_grain_UM; %script for this (works on dat_modis)

        end
        
        

        if iprc_diff==1
            diff_dat = dat_modis;
            dat_modis = vars_map_out.P_save{2}; %No volcano case
            if icoarsen_diff==1
                %reset these since have already smoothed once
                gcm_Plat2D_edges_UM = vars_map_out.gcm_Plat2D_edges_UM;
                gcm_Plon2D_edges_UM = vars_map_out.gcm_Plon2D_edges_UM;
                gcm_Plat2D_UM = vars_map_out.gcm_Plat2D_UM;
                gcm_Plon2D_UM = vars_map_out.gcm_Plon2D_UM;

                coarsen_method = '2d smoothing filter'; N = 128; M=N;  %Smoothing
                coarsen_method = 'gridpoint blocks';N = 64; M=N; %Coarse graining to NxM gridpoint blocks (rather than smoothing)
                %            coarsen_method = 'lat lon blocks'; dlat_target=0.25; dlon_target=0.25;

                coarse_grain_UM; %script for this (works on dat_modis)
            end
            
            %Percentage change relative to smoothed no volcano data
            dat_modis(dat_modis<thresh_LWP)=NaN;
            dat_modis = 100*diff_dat./dat_modis;

            titlenam_driver = [UM_map.varname ' % difference for ' time_str ' for ' UM_map.labs_UM(1).l ' minus ' UM_map.labs_UM(2).l ' at ' z_str];
            units_str_plot = '%';

        else

            titlenam_driver = [UM_map.varname ' difference for ' time_str ' for ' UM_map.labs_UM(1).l ' minus ' UM_map.labs_UM(2).l ' at ' z_str];
            units_str_plot = UM_map.var_units_str;

        end
        

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
        

         
        mod_data_type='AMSRE';
        gcm_str_select='UM';
%        daynum_timeseries3_UM = [1:length(time)];
%        gcm_time_UTC_UM = [1:length(time)];
        
%        gcm_Plon2D_UM = lons2D_driver;
%        gcm_Plat2D_UM = lats2D_driver;                
                           
%        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(lats2D_driver,lons2D_driver);

         gcm_time_matlab_UM = 1;
         gcm_time_UTC_UM = 1;
         daynum_timeseries3_UM = 1;
         modisyear_timeseries3_UM = 1;


        
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

 
 
%% PDFs

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
titlenam = [UM_map.varname ' PDFs'];
xlab=[UM_map.varname  ' (' UM_map.var_units_str ')'];



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


clear ioverride_LWP_diff
catch caught_error
    clear ioverride_LWP_diff
    rethrow(caught_error)
end