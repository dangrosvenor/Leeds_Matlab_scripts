% Loops through all 4 seasons for all the years and makes a climatology for
% multiple models. Saves in Nd_UKESM_savefile in the form of 
%  Name              Size              Bytes  Class     Attributes
% 
%   Nd_UKESM          3x4             2655552  cell                
%   Ndatap_UKESM      3x4             2655552  cell                
%   UM_map            1x1                5346  struct 
%
% Where Nd_UKESM{imodel,iseason} is a cell array for the different models
% selected.

Nd_UKESM_savefile = '/home/disk/eos8/d.grosvenor/UM/UKESM/Nd_UKESM_seasonal_clim_2000-2007.mat';

%Nd_file_UKESM = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2014_allCF_CTT273_region_mask_20150416T063530.mat';



clear UM_map

UM_map.isave_plot=0;

iplot_maps=1;
plot_diff=0;
plot_pdf_diff=0;

UM_map.varname = 'Nd_zweight'; %The name of the variable in the nc file
UM_map.VAR_NAME_STR='Nd_GCM'; %The name part in the filename to replace VAR_NAME with
            %OR, if the VAR_NAME is not in the supplied filename then will
            %just use the supplied one - make sure it is correct!
            
%UM_map.flag2=''; %replacing the old flag for whether is a mat file or not - set to '' if not.
UM_map.flag2='load_mat'; %replacing the old flag for whether is a mat file or not - set to '' if not.

UM_map.var_units_str = 'cm^{-3}';

UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov

UM_map.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

UM_map.time_tol='same month'; %set to this to match by month and year

%Setting time in the time loop below
%tol_mins=1;
%time_single = datenum('01-May-2007 00:00'); 
%UM_map.time_range(1) = time_single;  % - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;


UM_map.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
UM_map.thresh_LWP_mask = 20;


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

UM_map.noplot=1; %suppress the plotting part - just read in the data.

% time_single = datenum(['01-Jul-2004 00:00']); datestr(time_single)
% UM_map.time_range(1) = time_single;
% vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);
% 


%% Loop over time for creating seasonal data

% Loop over time - create seasonal climatologies
% Requeseted the files from 2000-2008, but at first hadn't used the 360 day
% calendar. Now that I have used 360 day calendar, we have data for 16th of every month
% from 2000-2007 (Jan to Dec).
% Need to decide what to do RE DJF ,etc. Prob is ok to use Jan and Feb 2000
% and Dec 2007 since am just doing a climatology
% (7 years). Should be enough for Nd, which doesn't very that much I think
% - can check this.
% So, just using DJF of each year, e.g. Dec, Jan and Feb of 2000.

years_clim = [2000:2007];
im=1; clear months
months{im}='Dec'; im=im+1;
months{im}='Jan'; im=im+1;
months{im}='Feb'; im=im+1;
months{im}='Mar'; im=im+1;
months{im}='Apr'; im=im+1;
months{im}='May'; im=im+1;
months{im}='Jun'; im=im+1;
months{im}='Jul'; im=im+1;
months{im}='Aug'; im=im+1;
months{im}='Sep'; im=im+1;
months{im}='Oct'; im=im+1;
months{im}='Nov'; im=im+1;



Ndatap = length(years_clim)*3; %3 months for each year

for iseason_clim=1:4
    ifirst=1; %whether we are starting a new season for averaging
    for iyear_clim=1:length(years_clim)
        for imonth_clim=1:3
%             %Special case for Dec of each year - want to use the current
%             %year, whereas for all other months will use the next year
%             if iseason_clim==1 & imonth_clim==1
%                 year = years_clim(iyear_clim);
%             else
%                 year = years_clim(iyear_clim+1);                 
%             end

            year = years_clim(iyear_clim);
            im = (iseason_clim-1)*3 + imonth_clim;
            time_single = datenum(['01-' months{im} '-' num2str(year) ' 00:00']); datestr(time_single)
            UM_map.time_range(1) = time_single;        
            vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);            
            
            for imodel=1:3
                if ifirst==1
%                    dat_season{imodel}=zeros(size(vars_map_out.P_save{imodel}));
                    dat_season{imodel}=zeros([length(years_clim)*3 size(vars_map_out.P_save{imodel})]);                    
                end
%                dat_season{imodel} = dat_season{imodel} + vars_map_out.P_save{imodel} / Ndatap;
                i_ind = (iyear_clim-1)*3 + imonth_clim;
%                dat_season{imodel}(iyear_clim,imonth_clim,:,:) = vars_map_out.P_save{imodel};                
                dat_season{imodel}(i_ind,:,:) = vars_map_out.P_save{imodel};                  
            end
            ifirst=0;
        end 
        
        
    end
    
    for imodel=1:3
        Nd_UKESM{imodel,iseason_clim} = dat_season{imodel};
        [Nd_UKESM{imodel,iseason_clim}, Ndatap_UKESM{imodel,iseason_clim}] = meanNoNan(dat_season{imodel},1);        
        %cut out data if there were any months with NaN for that location
        %(essentially the same as doing a straight mean)
        i_cut = find(Ndatap_UKESM{imodel,iseason_clim} < length(years_clim)*3);
        Nd_UKESM{imodel,iseason_clim}(i_cut) = NaN;
    end
            
end
    

save(Nd_UKESM_savefile,'Nd_UKESM','Ndatap_UKESM','UM_map');






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
  
