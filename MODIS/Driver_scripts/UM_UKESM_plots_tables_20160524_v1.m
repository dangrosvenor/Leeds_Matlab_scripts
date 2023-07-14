
%The stats are stored in here.
clear savefile_gridded_UKESM
UKESM_dir = '/home/disk/eos8/d.grosvenor/UM/UKESM/';
savefile_gridded_UKESM{1} = ['UKESM_vs_MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2008_allCF_CTT273_region_mask_20150416T063530.mat'];
savefile_gridded_UKESM{2} = ['UKESM_vs_MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2014_allCF_CTT273_region_mask_20150416T063530.mat'];
savefile_gridded_UKESM{3} = ['UKESM_vs_MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2008_CF80_CTT273_20150416T063530.mat'];
savefile_gridded_UKESM{4} = ['UKESM_vs_MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2014_CF80_CTT273_20150416T063530.mat'];



for iUKESM=1:length(savefile_gridded_UKESM)
    load_file{iUKESM} = [UKESM_dir savefile_gridded_UKESM{iUKESM}];
    load(load_file{iUKESM});

    var_label = 'RMS of % bias';
    var_label = 'Mean of % bias';
    var_label = 'Mean of abs bias';
%    var_label = 'Mean of abs bias as %';
%    var_label = 'RMS of abs bias';
    var_label = 'Correlation coefficient';
    
    switch var_label
        case 'RMS of % bias'
            var_str = 'Nd_UKESM_bias_prc_RMS';
        case 'Mean of % bias'
            var_str = 'Nd_UKESM_bias_prc_mean';
        case 'Mean of abs bias'
            var_str = 'Nd_UKESM_bias_abs_mean';
        case 'Mean of abs bias as %'
            var_str = 'Nd_UKESM_bias_abs_mean_prc';
        case 'RMS of abs bias'
            var_str = 'Nd_UKESM_bias_abs_RMS';
        case 'Correlation coefficient'
            var_str = 'Nd_UKESM_corr';
    end
    
    fprintf(1,'\n%s\n',savefile_gridded_UKESM{iUKESM});
    fprintf(1,'%s; %s; %s; %s; %s;\n',var_label,'DJF','MAM','JJA','SON'); 
    for i=1:3
        fprintf(1,'%s; ', UKESM.UM_map.labs_UM(i).l);
        for j=1:4

            fprintf(1,'%f; ', eval(['UKESM.' var_str '{i,j};']) )
        
        end
        fprintf(1,'\n');
    end
end




return

%% Do some plots, just for just CF>80 2000-2008 MODIS data

isave_plot=0;
    
%Set what to plot
plot_case_UKESM = 'UKESM map';
%plot_case_UKESM = 'MODIS';
plot_case_UKESM = 'UKESM difference';

%Reload the UKESM structure for CF>80 2000-2008 MODIS data (3rd file)
iUKESM=3;
load_file{iUKESM} = [UKESM_dir savefile_gridded_UKESM{iUKESM}];
load(load_file{iUKESM});

% MODIS data file(s)
MOD_dir = ['/home/disk/eos1/d.grosvenor/saved_misc_mat_files/'];
% The MODIS file used is already stored in UKESM
MOD=load([MOD_dir UKESM.Nd_file_MODIS]);

is=1;
season_str{is}='DJF'; is=is+1;
season_str{is}='MAM'; is=is+1;
season_str{is}='JJA'; is=is+1;
season_str{is}='SON'; is=is+1;

switch plot_case_UKESM
    case 'MODIS'
        nModels=1;
    otherwise
        nModels=3;    
end
        

for iseason=1:4

    %Interpolate the UM data onto the MODIS grid
    for imodel=1:nModels
        
        data_str = ['MOD.' UKESM.MOD_str season_str{iseason}];
        mod_dat = eval(data_str);



% gcm_Plat2D_edges_UM = UKESM.UM_grid.gcm_Plat2D_edges_UM;
% gcm_Plon2D_edges_UM = UKESM.UM_grid.gcm_Plon2D_edges_UM;
% gcm_Plat2D_UM = UKESM.UM_grid.gcm_Plat2D_UM;
% gcm_Plon2D_UM = UKESM.UM_grid.gcm_Plon2D_UM;

% gcm_Plat2D_edges_UM = UKESM.UM_grid.gcm_Plat2D_edges_UM;
% gcm_Plon2D_edges_UM = UKESM.UM_grid.gcm_Plon2D_edges_UM;

%Mode is in MODIS grid at the moment.
[gcm_Plat2D_UM, gcm_Plon2D_UM] = meshgrid(MOD.MLAT,MOD.MLON); %UKESM.UM_grid.gcm_Plat2D_UM;
gcm_Plon2D_UM = gcm_Plon2D_UM';
gcm_Plat2D_UM = gcm_Plat2D_UM';

[gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);


 %--- run the file to set up the defaults
        plot_global_maps_defaults  



%% What to plot is set here.
switch plot_case_UKESM
    case 'UKESM difference'
        dat_modis = 100*(UKESM.Nd_UKESM_grid{imodel,iseason} - mod_dat) ./ mod_dat;
        titlenam_driver = [UKESM.UM_map.varname ' % difference for ' time_str ' ' season_str{iseason} ' for ' UKESM.UM_map.labs_UM(imodel).l ' minus ' data_str];
        iset_min_clim=1;
        clim_min=-100;
        iset_max_clim=1;
        clim_max=100;
    case 'UKESM map'
        dat_modis = UKESM.Nd_UKESM_grid{imodel,iseason};
        titlenam_driver = [UKESM.UM_map.varname ' for ' time_str  ' ' season_str{iseason} ' for ' UKESM.UM_map.labs_UM(imodel).l ' (cm^{-3})'];
        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=260;        
    case 'MODIS'
        dat_modis = mod_dat;              
        titlenam_driver = ['MODIS Nd data (' data_str ') for ' time_str ' ' season_str{iseason} ];
        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=260;
end

       
        
        irestrict_domain=0; %UM_map.irestrict_domain_DRIVER; %whether to restrict the domain or not
        
%        thresh_LAT = UM_map.LAT_val_DRIVER;
%        thresh_LON = UM_map.LON_val_DRIVER;
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

      
%         
%         isave_plot=0;
%         iplot_markers=0;
        

        
        %Set various things
        time_str='2000-2008';
%        time = time_matlab(it_driver);
          %Round to the nearest minute as sometimes get 18:59:59
%        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
%ti        time_str = vars_map_out.time_str;
%        if iz>-1
%            z_str = ['at ' num2str(z_um(iz)/1000,'%.2f') ' km '];
%        else
%            z_str ='';
%        end
%        z_str = vars_map_out.z_str;
        z_str='';
        
        units_str_plot = UKESM.UM_map.var_units_str;
         
        mod_data_type='AMSRE';
        gcm_str_select='UM';
%        daynum_timeseries3_UM = [1:length(time)];
        daynum_timeseries3_UM = [1:length(UKESM.UM_map.time_range)];
        gcm_time_UTC_UM = daynum_timeseries3_UM;
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

 if isave_plot==1
     saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
     close(gcf);
 end
 
 
    end
    
end

 

return
 
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