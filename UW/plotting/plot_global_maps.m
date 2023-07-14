try %Catching errors so that flags are reset even if the program
    %aborts due to an error. Commands after the catch statement are executed
    %if an error now occurs. The flags are then reset and the error is
    %"rethrown"

    % Uses:-
    % *** script for choosing the time indices:-
    % time_inds_modisL3_timeseries3    ***
    % *** script for doing the SCREENING    :-
    % modisL3_screening_timeseries3    ***

    % For restricting to a certain region use irestrict_domain=1 - search
    % for thresh_lat and set these to the region required

    %  -----------------   flags to check if doesn't work  --------------------
    %proj_type='polar';   proj_type='global oval'; 
    %iseaice=0;
    %detailed_BAS_coast=0;
    %ifull_swath=0;

    
    
    clear dat_modis2
    
    savenotes_filedir='/home/disk/eos1/d.grosvenor/notes/'; %path for dir in which to save notes about figures
    save_notes_filepath = [savenotes_filedir 'dummy.txt']; %set a blank dummy notes file so that old ones are not used.
    
    if exist('ioverride_plotglobal_thresh') & ioverride_plotglobal_thresh==1
        %clear ioverride_plotglobal_thresh  %don't clear yet - wait until last
        %time is used
    else  %otherwise use these values
        % ----------------------------------------------------------------------------------
        irestrict_domain=1 %whether to restrict the domain or not
        % ----------------------------------------------------------------------------------
        
         gcm_time_of_day_select=0; %is set in the plot type specification below. Set to 2 for time
    %screening based on local time (but will need to change the way that
    %time_inds_average is used
    
        inew_cticks=0; %
        
        i_ocean_only_SCREEN = 0; %Whether to exclude non-ocean points for MODIS (based on amsre landmask)
        
        cont_col_str ='k';
        iplot_mgrid_lines=1; %whether to plot the grid lines for maps using m_grid
        ioverride_ticks=0;
    
    
    end
    

    % -----------------   default  flags   --------------------
    
    if exist('gcm_str_select')
        gcm_str=gcm_str_select
    else
        disp('WARNING - no gcm_str_select is set');
    end

    if ~exist('colormap_choose')  %only define if not already set
        colormap_choose=jet; %default
    end
    %    colormap_savejet_no_white;
    %    colormap_choose=jetNEW; %jet colormap with the white removed - now use
    %    black for NaN data, so this is not necessary

    iseaice=0; %is it a seaice plot?
    iplot_global=1;
    igcm_screen=0;
    iplot_mpace_flight_path=0; %is only set for level 2 plots


    contour_label=1; %Default =1 - whether to label contours (for icontour==1)
    
    if ~exist('ilabel_colorbar')
        ilabel_colorbar=0; %If =1 puts a label below the colorbar (need to specify the string, col_bar_lab_str)
    end
  

    
    detailed_BAS_coast=0; %flag to plot a detailed coastline from BAS (run convert_Antarctic_esri_file_lat_lon_mmap to
    %convert a shapefile to m-map lat/lon co-ords)
          
    ifull_swath=1; %flag to plot all of the swath for L2 
    %(not just the part with 5km lat,lon) by using the extrapolated fields
    
    i_dpcolor=1; %flag to use m_dpcolor insteat of m_pcolor - this plots all of the data by supplying
    %the edges for m_pcolor (i.e. one size bigger than the data) and
    %padding the data with NaNs at the end
    
    
    time_series_type='daily'; %default type of time data




    if exist('ioverride_plotglobal_thresh') & ioverride_plotglobal_thresh==1
        %if the override flag is set then use the values set elsewhere
        %will clear the flag after the last time it is used
    else  %otherwise use these values
        %type of projection
        iset_min_clim=0;  %flags for setting the colorscale limits or not
        iset_max_clim=0;

        proj_type='polar';
        proj_type='global oval'; %select this for anything but polar
%        stereographic

        %        data_select='latest_modis_data';
        data_select='specific_modis_data';
        ifilter_ndays=0 %flag for whether to cut out points for which there aren't many days
        icontour=0;
        
        noplot=0;
        
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/';
               
    end

   
    
    plot_mpace_flight_path=0; %default - set the flag above not this one



%get the screensize
    scrsz=get(0,'ScreenSize');
    %posit=[9 50 scrsz(3)/1.01 scrsz(4)/1.2];

    if ~exist('comp')
        comp='UWchallenger';
    end

    switch comp
        case 'uni'
            %    posit=[9 60 scrsz(3)/1.2 scrsz(4)/1.14];
            posit=[9 60 scrsz(3)/1.7 scrsz(4)/1.8];
        case 'lacieLap'
            posit=[9 50 scrsz(3)/1.46 scrsz(4)/2.07];
        case 'UWchallenger'
            posit=[9 50 scrsz(3)/1.46 scrsz(4)/1.46];
            fsize=12;
    end

    if exist('iover_ride_plot_global') & iover_ride_plot_global==1
        clear iover_ride_plot_global  %clear for the next time
    else    %otherwise set flags as the defaults
        inew_figure=1;
        supress_colorbar=0; %Setting to one stops the colorbar from appearing
        i_increase_font_size_map_figures_OFF = 0;  % setting to one stops this from running
    end

%choose the plot to be made
    switch data_select
        case 'latest_modis_data'
            P=test;

            MODIS_varname2_plot=MODIS_varname2;
            units_str_plot=units_str;
            modis_day_str_plot=modis_day_str;
            modis_year_str_plot=modis_year_str;

        case 'specific_modis_data'

            if exist('ioverride_plotglobal_thresh') & ioverride_plotglobal_thresh==1
                %clear ioverride_plotglobal_thresh  %reset on last use
            else  %otherwise use these values

                iasc_desc=0; %whether to use ascending (=1), descending (=2) or mean of both (=0) nodes

%% Select the data to plot
                
%indicate the type of data
                % L3 plots timeseries3 plots (data for every day)
                mod_data_type='timeseries3';
                
                


                
                %  ---- Note, for the ones below will want to add in some
                %  extra settings for mod_data_type='timeseries3 lambert'; 
                %  if want a lambert plot (if want a restricted region plot,
                %  etc.)
                
%                modis_data_plot= 'Bootstrap droplet number error of mean';

                modis_data_plot='Number of droplets cell values time mean'; %prob better to use 'Number of droplets cell values time mean - specific days'
                  %If want to select certain days fo the year, etc.

                %                 modis_data_plot='LWC at cloud top'; %from timeseries3
                %                 modis_data_plot='Cloud Depth'; %from timeseries3
                %                 modis_data_plot='Effective Radius'; %from timeseries3
                %                                 modis_data_plot='Optical Depth'; %from timeseries3

                %                modis_data_plot= 'Cloud Fraction Liq+Ice+Undetermined';
                %                modis_data_plot= 'Cloud Fraction Combined';
%                                modis_data_plot= 'Cloud Fraction Liquid';
                %                modis_data_plot= 'Cloud Fraction Ice';
                %                modis_data_plot= 'Cloud Fraction Undetermined';


                %                modis_data_plot='Number of Droplets Joint Histogram'; %calcualted using MODIS_N_H_calc
                %                modis_data_plot='No. pixels'; %calcualted from histograms
                %                modis_data_plot='Various'; %calcualted from histograms
                %                modis_data_plot='Number of droplets histogram time mean';

%                                modis_data_plot='Number of good data days for droplet number from cell values time mean';              
                %                modis_data_plot='Std dev LWP annual mean'; %standard deviation of LWP (W)
                %                modis_data_plot='LWP annual mean'; %standard deviation of
                %                LWP (W)
                %                modis_data_plot='Std dev Nd annual mean'; %standard deviation of Nd
                %                modis_data_plot='Max Nd annual';
                %                 modis_data_plot='Std dev Reff (\mum) annual mean'; %standard deviation of Reff


%                 modis_data_plot='Number of droplets at max LWC GCM'; mod_data_type='GCM';     gcm_time_of_day_select=1;
                % modis_data_plot='Number of droplets vertical mean over LWC GCM'; mod_data_type='GCM';     gcm_time_of_day_select=1;
%                 modis_data_plot='Cloud Fraction at max LWC GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;
%                modis_data_plot='LWP GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;
%                  modis_data_plot='RWP GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;
                  %TWP is LWP+RWP+IWP
%                  modis_data_plot='TWP GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;
                  %TLWP is LWP+RWP
%                  modis_data_plot='TLWP GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                  
%                  modis_data_plot='IWP GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                  
                  
%                  modis_data_plot='RWP ratio GCM (grid-box mean)'; mod_data_type='GCM';   gcm_time_of_day_select=2;
%                modis_data_plot='LWP GCM (grid-box mean, no CF screening)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                
%                modis_data_plot='In-Cloud LWP GCM (normalised by MODIS COSP CF)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                
%               modis_data_plot='In-Cloud LWP GCM (normalised by model CF)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                                
%                modis_data_plot='In-Cloud TLWP GCM (normalised by MODIS COSP CF)'; mod_data_type='GCM';   gcm_time_of_day_select=2;                                
%                modis_data_plot='Max LWC in Cloud Layer GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='Nd at max LWC within Cloud Layer GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                
%                modis_data_plot='LTS GCM'; mod_data_type='GCM';   gcm_time_of_day_select=2;                
%                modis_data_plot='qv700 GCM'; mod_data_type='GCM';   gcm_time_of_day_select=2;                
%                modis_data_plot='LTS ECMWF'; mod_data_type='GCM';   gcm_time_of_day_select=2; 
%                modis_data_plot='LTS bias from ECMWF'; mod_data_type='Monthly_data';   gcm_time_of_day_select=2;
%                modis_data_plot='qv700 ECMWF'; mod_data_type='GCM';   gcm_time_of_day_select=2;                                
               
%                modis_data_plot='Max Nd noCF GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='Max Nd in Cloud Layer GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='Mean Nd in Cloud Layer GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                
%                modis_data_plot='Max Nd in all non-screened layers GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
%                modis_data_plot='Max Nd, no CF screening minus CF screening'; mod_data_type='GCM';  gcm_time_of_day_select=2;

%                modis_data_plot='Max Nd noCF minus Max Nd in Cloud Layer GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='Max Nd COSP MODIS CF screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='LWP COSP MODIS CF screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                modis_data_plot='LWP COSP tau re screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                
%                modis_data_plot='Max Cloud Fraction after screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                 modis_data_plot='Max Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                 modis_data_plot='Max Low Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
%                 modis_data_plot='Max Mid Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;
%                 modis_data_plot='Max High Cloud Fraction no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2;                 
%                 modis_data_plot='Max Low Cloud Fraction no screening GCM bias from satellite'; mod_data_type='GCM';  gcm_time_of_day_select=2;

%                 modis_data_plot='IWP in Low Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                 modis_data_plot='IWP in Mid Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%                 modis_data_plot='IWP in High Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                 
%                 modis_data_plot='IWP in Mid+High Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                  
%                 modis_data_plot='LWP in Mid+High Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                  
%                 modis_data_plot='TWP in Mid+High Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                                   
%                 modis_data_plot='TWP in Mid Cloud no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                                                    
%                 modis_data_plot='TWP no screening GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;                                                                                     
%                modis_data_plot='Total surface precip rate GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2; igcm_screen=1; %screen for CTT
%                modis_data_plot = 'LWP removal timescale GCM'; mod_data_type='GCM';  gcm_time_of_day_select=2; igcm_screen=1; %screen for CTT
 %               modis_data_plot = '2d_pdf_lat_lon'; mod_data_type='GCM';  gcm_time_of_day_select=2; igcm_screen=1; %


%                modis_data_plot='Total Liquid Cloud Fraction COSP-MODIS GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;
%                modis_data_plot='Optical Depth Liquid COSP-MODIS GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;
%                modis_data_plot='Effective Radius Liquid COSP-MODIS GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;
%                modis_data_plot='Nd Liquid COSP-MODIS GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;
%                modis_data_plot='Low Cloud Fraction CALIPSO'; mod_data_type='CALIPSO';  gcm_time_of_day_select=0;  gcm_str='CALIPSO_monthly';
%                modis_data_plot='Low Cloud Fraction COSP-CALIPSO GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;

%                modis_data_plot='Mid Cloud Fraction CALIPSO'; mod_data_type='CALIPSO';  gcm_time_of_day_select=0;
%                modis_data_plot='Mid Cloud Fraction COSP-CALIPSO GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;

%                modis_data_plot='High Cloud Fraction CALIPSO'; mod_data_type='CALIPSO';  gcm_time_of_day_select=0;                
%                modis_data_plot='High Cloud Fraction COSP-CALIPSO GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=2; gcm_str=gcm_str_select;

%                modis_data_plot='Combined Cloud Fraction COSP-CALIPSO'; mod_data_type='CALIPSO';  gcm_time_of_day_select=0;

%                modis_data_plot='LWP AMSRE'; mod_data_type='AMSRE';  gcm_time_of_day_select=0; gcm_str='AMSRE';
                modis_data_plot='LWP AMSRE time3'; mod_data_type='timeseries3 lambert';  gcm_time_of_day_select=0; gcm_str='MODIS';                
%                modis_data_plot='In-cloud AMSRE LWP (normalised by MOD35 CF)'; mod_data_type='timeseries3 lambert';  gcm_time_of_day_select=0; gcm_str='MODIS';                
                modis_data_plot='AMSRE clear-sky bias map'; mod_data_type='timeseries3 lambert';  gcm_time_of_day_select=0; gcm_str='MODIS';                


%                modis_data_plot='Total Cloud Fraction COSP-CALIPSO GCM'; mod_data_type='CALIPSO';  gcm_time_of_day_select=1;
%                modis_data_plot='TOA SW Cloud Forcing COSP-CERES GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=1;                
%                modis_data_plot='TOA LW Cloud Forcing COSP-CERES GCM'; mod_data_type='GCM_COSP';  gcm_time_of_day_select=1;                                
%                modis_data_plot='Cloud-Sat warm rain precip rate'; mod_data_type='CLOUDSAT-PRECIP';  gcm_time_of_day_select=1; gcm_str='CLOUDSAT_PRECIP';                               
                modis_data_plot='Cloud-Sat warm rain fraction precip rate'; mod_data_type='CLOUDSAT-PRECIP';  gcm_time_of_day_select=1;
%                modis_data_plot='Cloud-Sat warm rain precip rate BIAS (GCM minus CLOUDSAT)'; mod_data_type='CLOUDSAT-PRECIP';  gcm_time_of_day_select=1;                                                
%                modis_data_plot='Cloud-Sat all rain precip rate'; mod_data_type='CLOUDSAT-PRECIP';  gcm_time_of_day_select=1;                                                
%                modis_data_plot='Cloud layer depth GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;

%                 modis_data_plot='Re3.7 minus re2.1 mockL3';  mod_data_type='timeseries3 lambert';
%                 modis_data_plot='Re1.6 minus re2.1 mockL3';  mod_data_type='timeseries3 lambert';
%                 modis_data_plot='Re3.7 mockL3';  mod_data_type='timeseries3 lambert';
%                 modis_data_plot='Re2.1 mockL3';  mod_data_type='timeseries3 lambert';

%                 modis_data_plot='Mock L3 Re2.1 minus L3'; mod_data_type='timeseries3 lambert';
%                  modis_data_plot='Mock L3 Tau minus L3'; mod_data_type='timeseries3 lambert';
%                  modis_data_plot='Re1.6 mockL3 minus POLDER'; mod_data_type='timeseries3 lambert';                  
                  
                 
% modis_data_plot='IWP GCM (old)'; mod_data_type='GCM';  gcm_time_of_day_select=1;
% modis_data_plot='LWP GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1; 
% modis_data_plot='IWP to TWP ratio GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1; 
%    modis_data_plot='Cloud Top Temperature GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%    modis_data_plot='Cloud Top Height GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;
%    modis_data_plot='Cloud Base Height GCM'; mod_data_type='GCM';  gcm_time_of_day_select=1;    
%    modis_data_plot='Cloud Top Temperature at max LWC (all layers) GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;    
%     modis_data_plot='Surface Temperature GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;        
%     modis_data_plot='Land Mask GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1; 
%     modis_data_plot='Surface Height GCM'; mod_data_type='GCM';
%     gcm_time_of_day_select=1;
%     modis_data_plot='Surface Height from Geopotential GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;     
%     modis_data_plot='Surface Pressure GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;  
%     modis_data_plot='Lat spacing GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;  
%     modis_data_plot='Lon spacing GCM'; mod_data_type='GCM'; gcm_time_of_day_select=1;  

% **** POLDER / PARASOL Reff ****
%       modis_data_plot='POLDER effective radius - using individual days'; mod_data_type='POLDER_daymeanALL';
%       modis_data_plot='POLDER effective radius - using colocated data'; mod_data_type='timeseries3 lambert';
%       modis_data_plot='POLDER effective radius - using mean over all years'; mod_data_type='POLDER_daymean';  
%        modis_data_plot='POLDER effective radius - using individual days from VOCALS save file'; mod_data_type='timeseries3 lambert';
% 
%       modis_data_plot='MODIS 1.6\mum minus POLDER effective radius - using mean POLDER over all years'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'
%       modis_data_plot='MODIS 2.1\mum minus POLDER effective radius - using mean POLDER over all years'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'
%       modis_data_plot='MODIS 3.7\mum minus POLDER effective radius - using mean POLDER over all years'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'

%      modis_data_plot='MODIS 1.6\mum minus POLDER effective radius - using colocated POLDER'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'
%      modis_data_plot='MODIS 2.1\mum minus POLDER effective radius - using colocated POLDER'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'
%      modis_data_plot='MODIS 3.7\mum minus POLDER effective radius - using colocated POLDER'; mod_data_type='POLDER_daymean'; gcm_str='MODIS'; mod_data_type='timeseries3 lambert'


% *** GOES ***

%       modis_data_plot='GOES Nd'; mod_data_type='L2 swath';  %gcm_time_of_day_select=1; 

%       
% **** Looking at specific day(s) from timeseries3 data *****
      
                % L3 plots timeseries3 plots - choosing specific time ranges for the mean
%               modis_data_plot='Number of droplets cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%               modis_data_plot='Number of droplets differences between datasets'; mod_data_type='timeseries3 lambert';               
%                modis_data_plot='Number of droplets cell values time MINIMUM - specific days'; mod_data_type='timeseries3 lambert';                
%                modis_data_plot='Max Solar Zenith Angle cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Min Solar Zenith Angle cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Max minus Min Solar Zenith Angle cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Number of droplets 1.6um cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Number of droplets 3.7um cell values time mean - specific days'; mod_data_type='timeseries3 lambert';                                
%                modis_data_plot='Number of droplets 3.7um minus 2.1um cell values time mean - specific days'; mod_data_type='timeseries3 lambert';                                
%                 modis_data_plot= 'Cloud_Fraction_Day_Mean'; mod_data_type='timeseries3 lambert'; igcm_screen=0; %MOD35 daytime cloud fraction
                                
%                 modis_data_plot='Nd ALL or low SZA - specific days'; mod_data_type='timeseries3 lambert';                                
%                modis_data_plot='Percentage Nd difference ALL vs low SZA - specific days'; mod_data_type='timeseries3 lambert';                                
%                modis_data_plot='Nd ALL SZA, but no low SZA - specific days'; mod_data_type='timeseries3 lambert';                                                

%                modis_data_plot='Number of good swaths for droplet number from mock L3, selected days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Mean number of swaths with good Nd per day from mock L3, selected days'; mod_data_type='timeseries3 lambert';
%                modis_data_plot='Max number of swaths in one day from mock L3, selected days'; mod_data_type='timeseries3 lambert';                
%               modis_data_plot='Min of the max daily SZA over period from mock L3, selected days'; mod_data_type='timeseries3 lambert';                
%               modis_data_plot='Max of the min daily SZA over period from mock L3, selected days'; mod_data_type='timeseries3 lambert';
%               modis_data_plot='Min of the min daily SZA over period from mock L3, selected days'; mod_data_type='timeseries3 lambert';               
%               modis_data_plot='Max of the mean daily SZA over period from mock L3, selected days'; mod_data_type='timeseries3 lambert';               
%                modis_data_plot='Number of different days with good Nd from mock L3, selected days'; mod_data_type='timeseries3 lambert';


%                      modis_data_plot='LWP cell values time mean (grid-box mean) - specific days'; mod_data_type='timeseries3 lambert';
%                      modis_data_plot='LWP cell values time mean (grid-box mean using MOD35) - specific days'; mod_data_type='timeseries3 lambert';
%                      modis_data_plot='LWP cell values time mean (in-cloud avearge) - specific days'; mod_data_type='timeseries3 lambert';
 
                      
                %                modis_data_plot='LWC max time mean - specific days'; mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Cloud depth cell values time mean - specific days'; mod_data_type='timeseries3 lambert';
%%                                modis_data_plot='Liquid Cloud Fraction time mean - specific days'; mod_data_type='timeseries3 lambert';
                                         %MOD06
%                                modis_data_plot='Ice Cloud Fraction time mean - specific days'; mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Combined Cloud Fraction time mean - specific days'; mod_data_type='timeseries3 lambert';                                
%                                modis_data_plot='Liquid Cloud Fraction Daily from mockL3 time mean - specific days'; mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Liquid Cloud Fraction Daily from L3 portion time mean - specific days'; mod_data_type='L3 time segment';
%                                modis_data_plot='Number of droplets cell values from L3 portion time mean - specific days'; mod_data_type='L3 time segment';
%                                modis_data_plot='Number of droplets cell values Daily from mockL3 time mean - specific days'; mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Number of droplets multi-year monthly means OLD'; mod_data_type='Monthly_data';
%                                modis_data_plot='Number of days multi-year monthly means'; mod_data_type='Monthly_data';
%                                modis_data_plot='Number of droplets multi-year monthly medians'; mod_data_type='Monthly_data';
%                               modis_data_plot='Number of droplets multi-year monthly means'; mod_data_type='Monthly_data';                
%                               modis_data_plot='LWP multi-year monthly means'; mod_data_type='Monthly_data';                                               
%                modis_data_plot='Sea ice areal coverage'; mod_data_type='Sea-ice'; iseaice=1;

                %                 modis_data_plot='Cloud Top Temperature Minimum time mean - specific days'; mod_data_type='timeseries3 lambert';
                %                 modis_data_plot='Cloud Top Pressure Minimum time max - specific days'; mod_data_type='timeseries3 lambert';
                %                 modis_data_plot='Cloud Top Pressure Minimum time mean - specific days'; mod_data_type='timeseries3 lambert';
                %                 modis_data_plot='Cloud Top Pressure Standard Deviation time mean - specific days'; mod_data_type='timeseries3 lambert';
                %                 modis_data_plot='Number of Droplets Joint Histogram mean from selected days'; mod_data_type='timeseries3 lambert';
                %                 modis_data_plot='Number of good data days for droplet number from cell values time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                 modis_data_plot='Mean Sensor ZA time mean, selected days';  mod_data_type='timeseries3 lambert';                
%                                 modis_data_plot='Mean SZA time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Max SZA time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                 modis_data_plot='Min SZA time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                  modis_data_plot='Max - Min SZA time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Max Sensor ZA time mean, selected days';  mod_data_type='timeseries3 lambert';
%                                modis_data_plot='Min Sensor ZA time mean, selected days';  mod_data_type='timeseries3 lambert';                                
%                modis_data_plot='Time max SZA from mockL3, selected days';  mod_data_type='timeseries3 lambert';
                                                 
%                                 modis_data_plot='Cloud Fraction Liquid Npixels time mean, selected days';  mod_data_type='timeseries3 lambert';
                                 
                %                 modis_data_plot='Cloud Depth Joint Histogram mean from selected days'; mod_data_type='timeseries3 lambert';
                %                  modis_data_plot='Cloud Depth timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
                %                  modis_data_plot='LWP Joint Histogram mean from selected days'; mod_data_type='timeseries3 lambert';
%                                  modis_data_plot='Cloud Top Temperature timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
                %                   modis_data_plot='Cloud Fraction mean from selected days'; mod_data_type='timeseries3 lambert';
%                                   modis_data_plot='Cloud Optical Depth timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
%                                   modis_data_plot='Cloud Effective Radius timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';

%                    modis_data_plot='MODIS Cloud Fraction minus CALIPSO CF'; mod_data_type='timeseries3 lambert';
                                 %can choose MOD06 or MOD35
%                    modis_data_plot='MODIS LWP minus AMSRE'; mod_data_type='timeseries3 lambert'; igcm_screen=0;
                                 %can choose MOD06 or MOD35                                 
                     %modis_data_plot='Fraction of days remaining after screening'; mod_data_type='timeseries3 lambert'; igcm_screen=0;

  %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7' in pdf2d_plot_commands)
%modis_data_plot='Cloud Effective Radius 1.6\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
%modis_data_plot='Cloud Effective Radius 2.1\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
%modis_data_plot='Cloud Effective Radius 3.7\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
%modis_data_plot='Cloud Effective Radius 3.7\mum minus 2.1\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';
%modis_data_plot='Cloud Effective Radius 1.6\mum minus 2.1\mum timeseries3 mean from selected days'; mod_data_type='timeseries3 lambert';


                %                 modis_data_plot='Cloud Top LWC Joint Histogram mean from selected days'; mod_data_type='timeseries3 lambert';
                % L3 single day plots


                %        modis_data_plot= 'Number from mode of 2D histogram - single day'; mod_data_type='daily';
                %        modis_data_plot= 'Number of Droplets';
                %        mod_data_type='daily';
                
                
%               modis_data_plot='SMOS soil moisture - specific days'; mod_data_type='SMOS'; %mod_data_type='timeseries3 lambert';                



%%  -------------------    L2 1km plots   -----------------------------------

                           modis_data_plot='Number of droplets (cm^{-3}) L2 swath'; mod_data_type='L2 swath';                        
                %                 modis_data_plot='N_d Uncertainty (cm^{-3}) L2 swath'; mod_data_type='L2 swath';
                %                 modis_data_plot='N_d % Uncertainty L2 swath'; mod_data_type='L2 swath';
%                                                 modis_data_plot='Cloud Fraction L2 swath'; mod_data_type='L2 swath';
                %                                 modis_data_plot='Cloud Depth (m) L2 swath'; mod_data_type='L2 swath';
                               %modis_data_plot='Optical Depth L2 swath'; mod_data_type='L2 swath';
%                               modis_data_plot='Latitude L2 swath'; mod_data_type='L2 swath';   
%                               modis_data_plot='Longitude L2 swath'; mod_data_type='L2 swath';                                  
                %                 modis_data_plot='Optical Depth Uncertainty L2 swath'; mod_data_type='L2 swath';
%                                modis_data_plot='Optical Depth % Uncertainty L2 swath'; mod_data_type='L2 swath';
%                               modis_data_plot='Reff (\mum) L2 swath'; mod_data_type='L2 swath';
                %                 modis_data_plot='Reff Uncertainty (\mum) L2 swath'; mod_data_type='L2 swath';
%                                 modis_data_plot='Reff % Uncertainty L2 swath'; mod_data_type='L2 swath';
%                                 modis_data_plot='Reff Difference (\mum) L2 swath'; mod_data_type='L2 swath';
                                 %modis_data_plot='Phase flag L2 swath 1621'; mod_data_type='L2 swath';
%                                 modis_data_plot='Various flags L2 swath'; mod_data_type='L2 swath';
                %                 modis_data_plot='Surface Temperature L2 swath'; mod_data_type='L2 swath';
%                               modis_data_plot='Cloud Top Temperature L2 swath'; mod_data_type='L2 swath';
%modis_data_plot='Cloud Top Temperature Joint L2 swath'; mod_data_type='L2 swath';                               
                %                                 modis_data_plot='Brightness Temperature L2 swath'; mod_data_type='L2 swath';
%                            modis_data_plot='Solar Zenith Angle'; mod_data_type='L2 swath';
%                            modis_data_plot='Scan Time'; mod_data_type='L2 swath';                            
%                            modis_data_plot='Sensor Zenith Angle'; mod_data_type='L2 swath';
%                              modis_data_plot='Liquid Water Path'; mod_data_type='L2 swath';

%%   ------------  L2 5km plots
                %              modis_data_plot='Various 5km flags L2 swath'; mod_data_type='L2 swath 5km';
                %              modis_data_plot='Sensor zenith angle'; mod_data_type='L2 swath 5km';
                %              modis_data_plot='Brightness temp diff 29-31'; mod_data_type='L2 swath 5km';
                %              modis_data_plot='Brightness temp diff 31-32'; mod_data_type='L2 swath 5km';
                %              modis_data_plot='Brightness temp band 31'; mod_data_type='L2 swath 5km';
                %           modis_data_plot='Cloud Top Height (km) L2 swath'; mod_data_type='L2 swath 5km';                   





                %%% filter thresholds and criteria (for removing points with certain
                %%% criteria)



                %     ****   mock L3 (reg lat-lon) plots  ****
                %                                   modis_data_plot='Mean Nd (cm^{-3}) mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Mean Nd from mean Tau & Re (cm^{-3}) mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Mean Re (\mum) mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Mean Tau mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='N datapoints mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Nd distribution skewness mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Number of ice or undetermined pixels mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Number of rejected pixels mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Pixel rejection ratio mock L3'; mod_data_type='mock L3';
                %                                   modis_data_plot='Mean Nd uncertainty (%) mock L3'; mod_data_type='mock L3';
                %                                    modis_data_plot='Cloud Depth (m) L2 swath mock L3'; mod_data_type='mock L3';
                %                                    modis_data_plot='Cloud Fraction L2 swath mock L3'; mod_data_type='mock L3';
                %                                    modis_data_plot='LWP homogeneity factor (Wood 2006) L2 swath mock L3'; mod_data_type='mock L3';
                %                                    modis_data_plot='LWP homogeneity factor (Cahalan 1994) L2 swath mock L3'; mod_data_type='mock L3';


            end

    end

    switch mod_data_type
        case {'timeseries3','daily','timeseries3 lambert','Monthly_data'}
            %moved filtering_data_get to load_saved_modis_vars.m as it is slow & only need to do once for each
            %timseries3 dataset
            %        filtering_data_get %get some of the fields used for cutting out some of the points
            
            gcm_time_of_day_select=0;
            
            %have to be careful about the direction of order of the arrays
            %- using dlat and dlon like this should account for this
            dlat=(mean(diff(MLAT)));
            dlon=(mean(diff(MLON)));
            
%make lat lon arrays - MLAT and MLON are the centre-points, so make the
%cell edges here and add one to the end            
            [Plon2D,Plat2D]=meshgrid(MLON,MLAT);

            Plat_lin=[MLAT-dlat/2 MLAT(end)+dlat/2];
            Plon_lin=[MLON-dlon/2 MLON(end)+dlon/2];            
            [Plon2D_edges,Plat2D_edges]=meshgrid(Plon_lin,Plat_lin);

            i_dpcolor=1;


            % *** choosing time periods ***
            %    daynum_timeseries3_MODIS = [1:355 358:366]; %days available in timeseries3 -
            %    %this is now set in the load script
            %time_mean_str = 'calculate'; %the default - otherwise reset below

            %    days_required_for_mean = 315;

            %if want all the days then just select 1:366 (even if there are fewer days
            %in the year)
            
            daynum_timeseries3=daynum_timeseries3_MODIS;

            % **** external script ***
            time_inds_modisL3_timeseries3
            % **** *************** ***

            
            if exist('aqua_terra_timeseries3')
                iaqua=strcmp('AQUA',aqua_terra_timeseries3);
                iterra=strcmp('TERRA',aqua_terra_timeseries3);
                if sum(iaqua)>0
                    aqua_terra_str='AQUA';
                    if sum(iterra)>0
                        aqua_terra_str='AQUA & TERRA';
                    end
                elseif sum(iterra)>0
                    aqua_terra_str='TERRA';
                end

            else
                aqua_terra_str='';
            end

            modisyears_str='';
            if length(modisyear_timeseries3)>1
                modisyears = unique(modisyear_timeseries3(time_inds_average));
            else
                modisyears = modisyear_timeseries3;
            end

            for iyear=1:length(modisyears)
                modisyears_str = [modisyears_str ' ' num2str(modisyears(iyear))];
            end

            %             if exist('ioverride_plotglobal_thresh') & ioverride_plotglobal_thresh==1
            %             else
            %                 ifilter_ndays=1 %flag for whether to cut out points for which there aren't many days
            %             end





            if strcmp(time_mean_str,'calculate')==1
                time_mean_str=datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy');
            end

            
        case {'L2 swath','GCM','GCM_COSP','CLOUDSAT-PRECIP','L3 time segment','CALIPSO','Sea-ice','AMSRE','POLDER_daymean','SMOS'}
            %filtering_data_L2

            switch mod_data_type
                case {'POLDER_daymean','POLDER_daymeanALL','SMOS'}
                    %need to set this before running
                    %time_inds_modisL3_timeseries3
                    switch mod_data_type
                        case 'POLDER_daymean'
                            daynum_timeseries3=[1:366];
                        case {'POLDER_daymeanALL','SMOS'}
                            daynum_timeseries3=eval(['daynum_timeseries3_' gcm_str]);
                    end
                     
                    gcm_time_of_day_select=0;
                    
                    %script
                    time_inds_modisL3_timeseries3
                    
                   
                    gcm_time_UTC=eval(['gcm_time_UTC_' gcm_str ';']);
                    
                     Plat2D_edges=eval(['gcm_Plat2D_edges_' gcm_str ';']);
                     Plon2D_edges=eval(['gcm_Plon2D_edges_' gcm_str ';']);
%                    Plat2D_edges = gcm_Plat2D_edges;
%                    Plon2D_edges = gcm_Plon2D_edges;


                     Plat2D = eval(['gcm_Plat2D_' gcm_str ';']);
                     Plon2D = eval(['gcm_Plon2D_' gcm_str ';']); 
                    
                case{'Sea-ice'}

                    Plat2D_edges=seaice_lat;
                    Plon2D_edges=seaice_lon;
                    
                    Plat2D=seaice_lat;
                    Plon2D=seaice_lon;
                    
                    i_dpcolor=0;


                        
                case{'CLOUDSAT-PRECIP'}
                    if exist(['daynum_timeseries3_' gcm_str])
                        daynum_timeseries3=eval(['daynum_timeseries3_' gcm_str]);
                    end
%                    gcm_time_UTC=eval(['gcm_time_UTC_' gcm_str ';']);
                    
                    time_series_type='monthly calipso precip';     
                    time_inds_modisL3_timeseries3
                    
                    %[Plon2D_edges,Plat2D_edges]=meshgrid(Plon,Plat);
                    Plat2D_edges = Plat2D_matt_edges;
                    Plon2D_edges = Plon2D_matt_edges;
                    Plat2D = Plat2D_matt_centres;
                    Plon2D = Plon2D_matt_centres;
                    ifull_swath=0;
                    i_dpcolor=1;
                    


                case {'GCM','GCM_COSP','CALIPSO','AMSRE'}
                    daynum_timeseries3=eval(['daynum_timeseries3_' gcm_str]);                  
                    gcm_time_UTC=eval(['gcm_time_UTC_' gcm_str ';']);
                    
                     Plat2D_edges=eval(['gcm_Plat2D_edges_' gcm_str ';']);
                     Plon2D_edges=eval(['gcm_Plon2D_edges_' gcm_str ';']);
%                    Plat2D_edges = gcm_Plat2D_edges;
%                    Plon2D_edges = gcm_Plon2D_edges;


                     Plat2D = eval(['gcm_Plat2D_' gcm_str ';']);
                     Plon2D = eval(['gcm_Plon2D_' gcm_str ';']);   
                     
                     
                   
                    %[Plon2D_edges,Plat2D_edges]=meshgrid(Plon,Plat);
                    
                                  


%                    Plat = gcm_Plat2D_edges;
%                    Plon = gcm_Plon2D_edges;
                    ifull_swath=0;
                    igcm_screen=0;
                    
%                    i_dpcolor=1;
                    
                    
                    switch mod_data_type

                        case {'GCM','GCM_COSP'}
                            modisyear_timeseries3=eval(['modisyear_timeseries3_' gcm_str]);
                            
                            if length(findstr('AM3',gcm_str))>0
%                            switch gcm_str
%                                case {'AM3','AM3_CLUBB','AM3_CLUBBv2','AM3_CLUBBv2_COSP_100km'} %may want to fix AM3 so that it uses dpcolor
                                    %are things aligned properly like this?
                                    i_dpcolor=0;                                    
 %                               otherwise
                            else
                                    i_dpcolor=1;
                            end

                        case 'CALIPSO'
                            i_dpcolor=1;
                            time_series_type = 'CALIPSO';
                            
                        case 'AMSRE'
                            i_dpcolor=1;
                            time_series_type = 'AMSRE';
                            
                        case 'POLDER_daymean'
                            time_series_type = 'POLDER';                            

                    end
                    
                    switch mod_data_type
                        case 'GCM'
                            switch gcm_str
                                case 'ERAInt_old'
                                    %don't do time inds if ECMWF
                                otherwise
                                    time_inds_modisL3_timeseries3
                            end
                            
                        otherwise   
                            switch gcm_str
                                case {'UM','CALIPSO_monthly','GENERIC'}
                                otherwise
                                    time_inds_modisL3_timeseries3
                            end
                    end
                     
                case 'L2 swath'
                    if ifull_swath==0
                        Plat2D=Plat_L2; %use the values read in from the file - these are the same as those worked out above
                        Plon2D=Plon_L2; %Plat_L2 and Plon_L2 are created in filtering_data_L2
                        %centre points           

                        Plat2D_edges = Plat3_L2;
                        Plon2D_edges = Plon3_L2;
                        
                        i_dpcolor=0;
            
                    else
                        Plat2D_edges=Plat2_L2; 
                        Plon2D_edges=Plon2_L2;
                        Plat2D=Plat2_L2;
                        Plon2D=Plon2_L2;
                    end
                    
                    
                case 'L3 time segment'
%                    daynum_timeseries3=eval(['daynum_timeseries3_' gcm_str]);
                    gcm_time_of_day_select=0;
                    time_inds_modisL3_timeseries3
                   
%                    Plat = MLAT_L32007_Arctic;
%                    Plon = MLON_L32007_Arctic;
                    ifull_swath=0;
                    
                                %have to be careful about the direction of order of the arrays
            %- using dlat and dlon like this should account for this
            dlat=(mean(diff(MLAT_L32007_Arctic)));
            dlon=(mean(diff(MLON_L32007_Arctic)));
            
%make lat lon arrays - MLAT and MLON are the centre-points, so make the
%cell edges here and add one to the end            
            [Plon2D,Plat2D]=meshgrid(MLON_L32007_Arctic,MLAT_L32007_Arctic);

            Plat_lin=[MLAT_L32007_Arctic-dlat/2 MLAT_L32007_Arctic(end)+dlat/2];
            Plon_lin=[MLON_L32007_Arctic-dlon/2 MLON_L32007_Arctic(end)+dlon/2];            
            [Plon2D_edges,Plat2D_edges]=meshgrid(Plon_lin,Plat_lin);

            i_dpcolor=1;
            
            
                    
%                    [Plon2D_edges,Plat2D_edges]=meshgrid(MLON_L32007_Arctic,MLAT_L32007_Arctic);
            end

            if iplot_mpace_flight_path==1
                plot_mpace_flight_path=1;
            end

        case 'mock L3'
            switch box_type
                case 'NxN pixel square';
                    %                    Plat2D_edges=LAT_mockL3;
                    %                    Plon2D_edges=LON_mockL3;
                    Plat2D_edges=LAT_mockL3_edge;
                    Plon2D_edges=LON_mockL3_edge;
                otherwise
                    [Plon2D_edges,Plat2D_edges]=meshgrid(LONS,LATS);
            end

            if iplot_mpace_flight_path==1
                plot_mpace_flight_path=1;
            end
        case 'L2 swath 5km'
            Plat2D_edges = lat;  %do the edges properly one day
            Plon2D_edges = lon;
            Plat2D = lat;
            Plon2D = lon;
    end

    if exist('thresh_str')
        thresh_str_old = thresh_str; %save the previous thresh_str
    end
    thresh_str='xxx';
%% choose thresholds and LAT LON and other screening ranges

    if exist('ioverride_plotglobal_loc') & ioverride_plotglobal_loc==1
        %    clear ioverride_plotglobal_thresh  %reset on last use
    else  %otherwise use these values
        thresh_ndays=0; %threshold no. days
%        thresh_ndays=5; %threshold no. days
%        thresh_ndays=1; %threshold no. days  
%        thresh_ndays=10; %threshold no. days          

        thresh_SZA=[65 90];
        thresh_SZA=[0 65];         
%        thresh_SZA=[0 90];        
        thresh_CF=0.8;
        thresh_CF=[0.8 1.00001];
%        thresh_CF=[0.1 1.00001];        
%        thresh_CF=[0 1.00001];        
        thresh_CF=[-0.01 0.8];   
        thresh_CF=[-0.01 0.01];                   
        thresh_CF=[-0.01 0.02];          
%        thresh_CF=[-0.01 0.05];           
%        thresh_CF=[-0.01 0.07];        
%        thresh_CF=[-0.01 0.09];           
%        thresh_CF=[-0.01 1.00001];


%        thresh_CF=[-0.01 10e-2];   %trying to look at AMSRE clear-sky bias
        
        thresh_NP=50;
%        thresh_NP=0;
        
        thresh_sensZA=45;
        thresh_sensZA=50;
        thresh_sensZA=40;
        thresh_sensZA=[0 41.4]; 
        thresh_sensZA=[0 90];         
        thresh_CTT = 273;
        thresh_CTT = [273-100 273+100];        
%        thresh_CTT = [273 273+100];    
%        thresh_CTT = [268 273+100];   
%        thresh_CTT = [273-10 273+100];  
        
        thresh_reff = [0 30];
        thresh_reff = [0 12];        
        
        thresh_CTP = [730 1000]; %Cloud top pressure (hPa)
        
        thresh_relAZ = [0 180];

        %    thresh_sensZA=80;
        thresh_maxSZA=[-1 81.4];
        thresh_minSZA=45;        
        thresh_stdSZA=[-1 0.5e9];
        thresh_stdSZA=[-1 0.5];
        
        thresh_dSZA=[-1 1e9];        
       thresh_dSZA=[-1 1];
        
        thresh_Nd_per_error = 100;
        thresh_Reff_per_error = 50;
        thresh_Reff_abs_error = 4;
        thresh_Reff = 30;
        
        gcm_CTT_thresh = [273.15 400];
        
        thresh_CTH = [-0.01 1e9]; %km
        thresh_CTH = [-0.01 3.2]; %km  - the approx height of 680 hPa (low cloud pressure threshold for ISCCP, CALIPSO, etc)
%        thresh_CTH = [3.2 20];
%        thresh_CTH = [6.5 20];        

        thresh_zeroCF = 0.05;
        
        thresh_stdW = [0 1e9];
%        thresh_stdW = [0 5];
        
        thresh_sigCTT = [0 1e9];
%        thresh_sigCTT = [0 1];  


      minfrac_CF = 0.9; %minimum fraction of the sampled points that had successful cloudy/clear/phase
      %determination (i.e. Npix/Nptot_mockL3 =
      %Cloud_Fraction_Liquid_Pixel_Counts./Cloud_Fraction_Liquid./Total_pixels
      % - restriction (2) as presented in the paper

      minfrac_NpNd = 0.9;        
       %Cloud_Fraction_Liquid_Pixel_Counts2.timeseries3./Cloud_Fraction_Liquid_Pixel_Counts.timeseries3
       %Fraction of points that remain after all previous filtering for
       %which we have an Nd retrieval. Restriction (4) in the SZA paper.
            
      thresh_NP_Nd = 50; %min no. of pixels required for an Nd, re, tau, etc measurement to count
      %(Uses Cloud_Fraction_Liquid_Pixel_Counts2). For this screening
      %usually
      %thresh_NP is the number of pixels that the swath must have covered -
      %i.e. the total number of pixels available. Only a portion of those 
      
        
%% Latitudes
        thresh_LAT = [70 74];
        %                      thresh_LAT = [72 73];


        %       thresh_LON = [-160 -135];
        %       thresh_LON = [-160 -140];
        %                       thresh_LON = [-152.5 -150];
        thresh_LAT = [67 75]; thresh_LON = [-180 -135]; %MPACE domain
        %         thresh_LAT = [70.33 71.33]; thresh_LON = [-154 -150]; %MPACE domain - zoomed into flight track
        % %        thresh_LAT = [70.33 71.33]; thresh_LON = [-158 -148]; %MPACE domain - zoomed into flight track even closer
        % %        thresh_LAT = [69 72]; thresh_LON = [-160 -146]; %MPACE
        %         thresh_LAT = [70 71.5]; thresh_LON = [-160 -146]; %MPACE
        %         thresh_LAT = [70 71.5]; thresh_LON = [-151 -146]; %MPACE - profiles to the east
        %         thresh_LAT = [70 71.5]; thresh_LON = [-153 -148]; %MPACE - profiles
        %         thresh_LAT = [70 71.5]; thresh_LON = [-153.5 -151.5]; %MPACE - profiles
        thresh_LAT = [60 64]; thresh_LON = [25 29]; %Sami's station (Puijo)
        thresh_LAT = [62.5 63.5]; thresh_LON = [25 29]; %Sami's station   (@ 62.909 N, 27.656 E)
%        thresh_LAT = [50 70]; thresh_LON = [10 40]; %Sami's station   
        
%                thresh_LAT = [-40 -10];  thresh_LON = [-100 -60];
                %VOCALS CAPT:- One used for CPT stuff in general 10th Jan
                %2013
                thresh_LAT = [-40 10];  thresh_LON = [-140 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                
                %To match Fig. 13 of Toniazzo 2011 paper 
%                 thresh_LAT = [-35 -10];  thresh_LON = [-90 -65]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                 
               %GOES domain to roughly match the area of the VOCALS GOES
               %images
%                thresh_LAT = [-33 -7];  thresh_LON = [-103 -67]; lon_ticks=[-10:5:-65]; lat_ticks=[-30:5:-5];%VOCALS CAPT
 

               %Region surroudning the UM domain for 12th Nov
                thresh_LAT = [-30 -10];  thresh_LON = [-90 -70];  lon_ticks=[-10:5:-65]; lat_ticks=[-30:5:-5];%VOCALS CAPT

                %new CPT domain extended up to California
%                thresh_LAT = [-45 45];  thresh_LON = [-150 -55]; 
                
                
%                thresh_LAT = [-60 60];  thresh_LON = [-180 -50]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT                %
                
%                thresh_LAT = [-50 50];  thresh_LON = [-160 -30]; lon_ticks=[-140:5:-50]; lat_ticks=[-40:10:10];%VOCALS CAPT
                
 %               thresh_LAT = [0 65]; thresh_LON = [-160 -60]; %Eastern US,
        %        UK western Europe
%                thresh_LAT = [0 85]; thresh_LON = [-180 180]; %N. hemisphere - whole globe    
        %        thresh_LAT = [-85 0]; thresh_LON = [-180 180]; %S. hemisphere - whole globe    
                
%                thresh_LAT = [0 65]; thresh_LON = [-180 180]; %N. hemisphere - whole globe
%                 thresh_LAT = [-65 0]; thresh_LON = [-180 180]; %S.hemisphere - whole globe
%                   thresh_LAT = [36 52]; thresh_LON = [-30 0]; %UK, Atlantic off SW coast of UK
%            thresh_LAT = [44 51]; thresh_LON = [-30 -6]; %UK, Atlantic off SW coast of UK

%                    thresh_LAT = [20 60]; thresh_LON = [-60 30]; %UK, Atlantic off SW coast of UK
                   
        %    thresh_LAT = [44 52]; thresh_LON = [-15 0]; %UK, Atlantic off
        %    SW coast of UK
        %    thresh_LAT = [36 65]; thresh_LON = [-15 60]; %mainland Europe

%                 thresh_LAT = [-10 60];  thresh_LON = [-36 60]; %Europe
%                 thresh_LAT = [30 80];  thresh_LON = [-90 80]; %UK and Iceland

%                 thresh_LAT = [-80 -40]; thresh_LON = [-200 -30];
%                 %Antarctica
%                 thresh_LAT = [-75 -60]; thresh_LON = [-80 -50];
%                 %Antarctica   
%                 thresh_LAT = [-72 -63]; thresh_LON = [-72 -55];
%                 %Antarctica    
%                 thresh_LAT = [-70 -65]; thresh_LON = [-72 -55]; %Used for paper. Antarctica Peninsula and Larsen C
%also set  detailed_BAS_coast=1 (run convert_Antarctic_esri_file_lat_lon_mmap to convert a shapefile to m-map lat/lon co-ords)                 
%                 thresh_LAT = [-69 -65]; thresh_LON = [-69 -64]; %Antarctica Peninsula - close up on Rothera 

%                 thresh_LAT = [-75 -40]; thresh_LON = [-80 -50]; %Antarctica   
%                 thresh_LAT = [-75 -40]; thresh_LON = [-90 -20];
%                 %Antarctica                    

%                 thresh_LAT = [-90 -50]; thresh_LON = [-180 180];  %Southern pole

%                 thresh_LAT = [40 80];  thresh_LON = [-90 80]; %UK and Iceland for Anja plots
%                 thresh_LAT = [45 75];  thresh_LON = [-60 40]; %UK and Iceland for Anja OMI vs model plots
%                 thresh_LAT = [48 80];  thresh_LON = [-70 40]; %UK and Iceland for Anja plots vs UM

%                 thresh_LAT = [0 40];  thresh_LON = [60 90]; %India     
%                 thresh_LAT = [-40 40];  thresh_LON = [-20 50]; %Africa                      

%                 thresh_LAT = [36.5 41.5];  thresh_LON = [-31 -25]; %Azores   


%requested region was 70-80N, -20 to 60E
%thresh_LAT = [70 80]; thresh_LON = [-20 60]; %Arctic summer requested
%region
%thresh_LAT = [65 85]; thresh_LON = [-20 60]; %Arctic summer requested region

%thresh_LAT = [-23 -16]; thresh_LON = [-90 -80]; %

%thresh_LAT = [-60 -10]; thresh_LON = [110 155]; %Australia
%thresh_LAT = [-69 10]; thresh_LON = [90 180]; %

%thresh_LAT = [0 60]; thresh_LON = [90 180]; %China



%thresh_LAT = [-65 -35]; thresh_LON = [-180 180]; %Southern Ocean
%thresh_LAT = [-60 -45]; thresh_LON = [50 100]; %Southern Ocean box
%thresh_LAT = [-65 -30]; thresh_LON = [0 240]; %Southern Ocean box

%Cape Grimm is 41S, 144.5E

% Hawaii - Amy's domain
%thresh_LAT = [12 25]; thresh_LON = [-170 -150]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii
%thresh_LAT = [0 25]; thresh_LON = [-180 -150]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii
%wider picture
%thresh_LAT = [-25 25]; thresh_LON = [-180 180]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Downwind of Hawaii


%Koike paper region
thresh_LAT = [25 40]; thresh_LON = [120 135]; lon_ticks=[-170:5:-150]; lat_ticks=[10:5:25]; %Japan/S. Korea

%thresh_LAT = [-90 90]; thresh_LON = [-180 180]; %global

%Spain for Rev Geo picture
thresh_LAT = [36 53]; thresh_LON = [-17 3]; lon_ticks=[-15:5:0]; lat_ticks=[40:5:50];
thresh_LAT = [40 52]; thresh_LON = [-16 0]; lon_ticks=[-15:5:0]; lat_ticks=[40:5:50];


if length(thresh_CF)==1
    thresh_CF(2)=1;
end
if length(thresh_CTP)==1
    thresh_CTP(2)=1100;
end


 end


dlon=10;
dlat=5;

%dlat=10;
%dlon=30;

lat_01 = ceil(thresh_LAT(1)/dlat)*dlat;
lon_01 = ceil(thresh_LON(1)/dlon)*dlon;

if abs(lat_01)<500

lon_ticks=[lon_01:dlon:thresh_LON(2)]; 
lat_ticks=[lat_01:dlat:thresh_LAT(2)]; 
end

      



        thresh_Ngood = 0;
        icolormap_cf_grey=0;




   








    if irestrict_domain==0
        if strcmp(mod_data_type,'timeseries3 lambert')==1
            mod_data_type='timeseries3';
        end
        if strcmp(mod_data_type,'Monthly_data')==1
            mod_data_type='Monthly_data_global';
        end
    end


    %data got from call to filtering_data_get (above)

    %cut some of the data out of the average
    %ihtot are the points to CUT OUT in this case (are set
    %to NaN so that they are not included in the means)

    ihtot=[]; thresh_str='xxx';
    ihtot_cf=[];

%% Selection of screening type    
    switch mod_data_type
        case {'timeseries3','daily','timeseries3 lambert'}
            if exist('ioverride_plotglobal_thresh') & ioverride_plotglobal_thresh==1
                %                clear ioverride_plotglobal_thresh %use the values set outside and reset the flag
            else  %otherwise use these values
%                screen_type='NP + CF + MAX sensZA';
           
                %                                    screen_type='NP + MAX sensZA';
                screen_type='NP + CF + MEAN sensZA';
%                screen_type='NP + CF + MAX solarZA';
                screen_type='NP + CF + MIN solarZA';                     
                screen_type='NP + CF + MEAN solarZA';
                screen_type='NP + CF';
                screen_type='NP + CF mockL3';
                screen_type='NP + CF mockL3 + meanCTH';
                screen_type='NP + CF mockL3 + meanCTH + meanCTT';
                screen_type='NP + CF mockL3 + meanCTH + meanCTT + meanReff';                
                screen_type='NP + CF mockL3 + meanCTH + meanCTT + meanReff + SZA';
                % screen_type='NP + CF + warm';
                % screen_type='NP + CF + min pressure';
                % screen_type='NP + CF + min pressure + min temp';
                % screen_type='NP + CF + mean pressure';
                % screen_type='NP + CF + (mean pressure - std_dev)'; %using same pressure threshold, but applying to the
                % mean CTP minus the std dev of the CTP (i.e. most data
                % will be at higher pressure than mean CTP minus std dev of
                % CTP
                % screen_type='NP + CF + (mean pressure - 2*std_dev)'; %using same pressure threshold, but applying to the

                % screen_type='NP + CF + mean pressure + mean temp';
                % screen_type='NP + warm';
%                screen_type='NP';
                %
%                screen_type ='NP + CF + MEAN sensZA + MEAN relAZ + MEAN solarZA + min_CTT';
                % screen_type='CF';
                screen_type ='CTH only';
%                screen_type ='CTH only with no zeroCF screening';
%                screen_type ='max CTH only with no zeroCF screening';
% --------------------------------------------------------------------------                
                %Screening for lots of things as of 9th Jan, 2014. For
                %level-3 data
%                screen_type = 'NP + CF_L3, no iceCF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW';
%                screen_type = 'NP + CF_L3, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff';                
%                screen_type = 'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff';
%                screen_type = 'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + max_CTH with no zeroCF screening + sigCTT +sigW + reff';                
                 screen_type = 'NP + CF_L3_MOD35';
% --------------------------------------------------------------------------

%                screen_type ='NP + CF_L3';
%                screen_type ='CTH only, include NaN CTHs';

%           screen_type = 'NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH';
%           screen_type = 'NP3 + CF3 + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH + sigCTT'; %

%screen_type = 'NP + CF_L3_MOD35, ice CF allowed + MEAN sensZA + MEAN relAZ + MEAN solarZA + mean_CTT + min_tau + mean_CTH with no zeroCF screening + sigCTT +sigW + reff'
%                screen_type='NP + CF_L3 + max_SZA + std_SZA';
%                screen_type='NP + CF_L3 + dSZA'; %calculate dSZA using dSZA=
                  %G14.dSZA=G14.Solar_Zenith_Maximum.timeseries3-G14.Solar_Zenith_Minimum.timeseries3;
                  %dSZA = time_match_data({daynum_timeseries3_MODIS,modisyear_timeseries3},G14.Date_Time_Swath.timeseries3,G14.dSZA);
                  %dSZA=flipdim(dSZA,1); %since the G14 data is backwards
                  %in LAT dimension
                  
               screen_type ='none';
            end


            %screening done in here now --->
            modisL3_screening_timeseries3
            % **************************************

        case 'Monthly_average'


        case {'L2 swath'}


            ihtot=[]; thresh_str='xxx';

            %                                ihtot = find( ~ ( cfL2>=thresh_CF ) ); thresh_str=['CF.GTE.' num2str(thresh_CF)];
            %                                ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error ) ); thresh_str=['Nd Perror.LT.' num2str(thresh_Nd_per_error)];
            %                                ihtot = find( ~ ( re_un<thresh_Reff_per_error ) ); thresh_str=['Reff Perror.LT.' num2str(thresh_Reff_per_error)];
            %                                ihtot = find( ~ ( re_un_abs<thresh_Reff_abs_error ) ); thresh_str=['Reff abs error.LT.' num2str(thresh_Reff_abs_error)];
            %                                ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & N>400) ); thresh_str=['Perror.LT.' num2str(thresh_Nd_per_error) '.AND.N.GT.400'];
            %                                ihtot = find( ~ ( re_un<thresh_Reff_per_error & (phase_flag==2) ) ); thresh_str=['Reff % error.LT.' num2str(thresh_Reff_per_error) ' AND liquid phase.'];  %                                 ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & re<thresh_Reff) ); thresh_str=['Perror.LT.' num2str(thresh_Nd_per_error) 'AND Re.LT.' num2str(thresh_Reff)];
            %            ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & (phase_flag==2) ) ); thresh_str=['Nd % error.LT.' num2str(thresh_Nd_per_error) ' AND liquid phase.'];
            %                                ihtot = find( ~ ( re_un_abs<thresh_Reff_abs_error & (phase_flag==3) ) ); thresh_str=['Reff abs error.LT.' num2str(thresh_Reff_abs_error) ' AND ice phase.'];
%                        ihtot = find( ~ ( phase_flag==2 ) ); thresh_str=[' liquid phase.'];
            %            ihtot = find( ~ ( phase_flag==3 & tau_bounds==0 )  ); thresh_str=[' ICE phase.AND.Tau within bounds'];
            %            ihtot = find( ~ ( phase_flag==2 & tau_bounds==0) ); thresh_str=[' liquid phase.AND.Tau within bounds'];

            %                                 ihtot = find( ~ ( phase_flag==2 & tau_bounds==0 & surface_flag==0) ); thresh_str=[' liquid phase.AND.Tau within bounds, ocean only'];
            %                                 ihtot = find( ~ ( CTT>thresh_CTT ) ); thresh_str=[' CTT.GT.' num2str(thresh_CTT)];
%4                        ihtot = find( ~ ( percent_error_Nd<thresh_Nd_per_error & phase_flag==2 & tau_bounds==0) ); thresh_str=['Perror.LT.' num2str(thresh_Nd_per_error) ' ,liquid phase.AND.Tau within bounds'];

             screen_type = 'Pixel screening L2';
             screen_type = 'Pixel screening L2 and CTT restriction';  
             screen_type = 'CTH L2 swath';
             screen_type = 'CTH + SZA L2 swath';             
             screen_type = 'none';
             
             %screening done in here now --->
             modisL3_screening_timeseries3
             % *******************************
             


        case 'mock L3'
            ihtot = find( ~ ( Np_mockL3>=thresh_Ngood ) ); thresh_str=[thresh_str_mock_L3 ' N datapoints after screening.GTE.' num2str(thresh_Ngood)];

        case 'GCM'
            %            screen_type='CF';
            %            screen_type='none';

            %script to scale by CF, screen for CF and LWC
            %gcm_process;
            
            if igcm_screen==1
                screen_type = 'gcm CTT-layer';
                screen_type = 'none';
                %screening done in here now --->
                modisL3_screening_timeseries3
                % **************************************

            else
                
                if ~exist('i_ice_screen')
                    i_ice_screen=1;
                end
                
                if ~exist('low_cf_thresh')
                    low_cf_thresh=[0.8 1];
                end
                
                if ~exist('icf_low_only')
                    icf_low_only=0;
                end
                
                 
                if ~exist('thresh_Nlevs')
                    thresh_Nlevs=0;
                end


                if i_ice_screen==0
                    %thresh_str=['CF.GTE.' num2str(thresh_CF(1)) '.AND.CF.LT.' num2str(thresh_CF(2))];
                    thresh_str=['CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
                        ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3} icf_low_only=' num2str(icf_low_only)];
                elseif i_ice_screen==1
                    thresh_str=['CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
                        ' icf_low_only=' num2str(icf_low_only) ' LIQUID only'];

                end

            end

            %             switch screen_type
            %                 case 'CF';
            %                     ihtot = find( ~ ( gcm_CF_maxliq(time_inds_average,:,:)>=thresh_CF(1) & gcm_CF_maxliq(time_inds_average,:,:)<thresh_CF(2) ) ); thresh_str=['CF.GTE.' num2str(thresh_CF(1)) '.AND.CF.LT.' num2str(thresh_CF(2))];
            %             end
    end

    %                        ihtot=[]; thresh_str='xxx';






%% Start of main switch for choosing what to plot
    switch modis_data_plot
        case 'Generic plot specified outside of script'
            %allow all the options to be overruled!
            
            %Run plot_global_maps_defaults (might want to change these)
            %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
            % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
            % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
            % mod_data_type='AMSRE';
            % Can set these times to one, or do properly depending on what
            % time screening you want, etc.
            
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                colormap_choose=jet; %default

                dat_modis=N_time3;

                dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

                %need to make sure we are applying the time_inds_average AFTER
                %NaN screening since the screening does NOT use time_inds_average

                dat_modis2 = dat_modis(:,:,time_inds_average);

                [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


                MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];

                units_str_plot='';
                title_info = '';
                icontour=0;
                inew_cticks=0;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;

                ifilter_ndays=0;
                thresh_ndays=15;

            end
            
        case 'AMSRE clear-sky bias map'
            %allow all the options to be overruled!
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                
                                
%                 filename = ['AMSRE_clear_sky_bias_9percent_VOCALS_region_20150918T065119'];
%                filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region_20150918T064749'];                         
%                filename = ['AMSRE_clear_sky_bias_5percent_VOCALS_region_20150918T065504'];        
%                filename = ['AMSRE_clear_sky_bias_3percent_VOCALS_region_20150918T065607'];
                 filename = ['AMSRE_clear_sky_bias_extrapolated_VOCALS_region_20150923T064911'];
                 
                 clear P_save_AMSRE_clear_sky

                 filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                 filename_savevars = [filedir_savevars filename '.mat'];
                 load(filename_savevars);
                 tag = '_AMSRE_clear_sky';
                 
                 icontour=1;
                 cont_dat_choose = 'calipso mid+highCF';
                 
                 

                colormap_choose=jet; %default

%                 dat_modis=P_save_AMSRE_clear_sky;
% 
%                 dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
% 
%                 %need to make sure we are applying the time_inds_average AFTER
%                 %NaN screening since the screening does NOT use time_inds_average
% 
%                 dat_modis2 = dat_modis(:,:,time_inds_average);
% 
%                 [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

                  P = P_save_AMSRE_clear_sky; 
                  Ndays2 = Ndays2_AMSRE_clear_sky;


                MODIS_varname2_plot = [modis_data_plot ' using ' filename];

                units_str_plot='g m^{-2}';
                ilabel_colorbar=1;
                col_bar_lab_str = units_str_plot;
               
                title_info = '';

                inew_cticks=0;


                iset_min_clim=1;
                clim_min=-20;
                iset_max_clim=1;
                clim_max=20;

                ifilter_ndays=0;
                thresh_ndays=15;

            end    
            
            
            if icontour==1
                
                % -- Use CALIPSO mid+high CF to define a certain region
% But need to grid the CALIPSO data onto the MODIS grid
            clh = meanNoNan(clhcalipso_monthly,1)/100;
            clm = meanNoNan(clmcalipso_monthly,1)/100;
            cll = meanNoNan(cllcalipso_monthly,1)/100;   

                cont_ints=[0.3 0.3];
                %    cont_ints=[0.35 0.35];
                %N.B. - this is changed belwow some fields

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso lowCF'
                        cont_dat = cll;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'MOD35 CF'
                        cont_dat = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
                        cont_ints=[0.5:0.05:0.8];
                    case 'MOD35 CF minus 1 std dev'
                        [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
                        cont_dat = mean_CF - std_CF;
                        cont_ints=[0.1:0.05:0.8];
                    case 'MOD06 CF'
                        cont_dat = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
                        cont_ints=[0.5:0.05:0.8];
                    case 'MOD06 CF minus 1 std dev'
                        [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
                        cont_dat = mean_CF - std_CF;
                        cont_ints=[0.1:0.05:0.8];
                end

            end
    
    
            
            
        case 'IWP in Low Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*iwp_isccp_low_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=1;

case 'IWP in Mid Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*iwp_isccp_mid_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=5;
            
        case 'IWP in High Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*iwp_isccp_high_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=10;

        case 'IWP in Mid+High Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*(iwp_isccp_high_' gcm_str '+iwp_isccp_high_' gcm_str ')']);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=10;

        case 'LWP in Mid+High Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*(lwp_isccp_high_' gcm_str '+lwp_isccp_high_' gcm_str ')']);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=1;
            
        case 'TWP in Mid+High Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*(lwp_isccp_high_' gcm_str '+lwp_isccp_mid_' gcm_str '+iwp_isccp_high_' gcm_str '+iwp_isccp_mid_' gcm_str ')']);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=100; %15

        case 'TWP in Mid Cloud no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['1e3*(lwp_isccp_mid_' gcm_str '+iwp_isccp_mid_' gcm_str ')']);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=100;
            
            
         case 'TWP no screening GCM'
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
%            dat_modis = eval(['1e3*(lwp_isccp_high_' gcm_str '+lwp_isccp_mid_' gcm_str '+iwp_isccp_high_' gcm_str '+iwp_isccp_mid_' gcm_str ')']);            
            
            dat_modis = eval(['1e3*(gcm_lwp_' gcm_str '+gcm_iwp_' gcm_str ')']);                        
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=200; %15
    
        case '2d_pdf_lat_lon'
            

                dat_modis = qh;
                
                %interpolate to the standard gcm lat lon grid.
                
                lat_bins_2D = (repmat(Ybins,[length(Xbins) 1]))';
                lon_bins_2D = repmat(Xbins,[length(Ybins) 1]); 
                
                X = eval(['gcm_Plon2D_' gcm_str]);
                Y = eval(['gcm_Plat2D_' gcm_str]);                
                
                P = griddata(lon_bins_2D,lat_bins_2D,dat_modis,X,Y); %interpolate GCM data onto same grid
                




%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            cont_dat_choose = 'qv';
            cont_dat_choose = 'LTS';
            switch cont_dat_choose
                case 'qv'
                    cont_dat = eval(['1e3*gcm_qv700_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:2:16];
                case 'LTS'
                    cont_dat = eval(['gcm_LTS_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:4:30];
            end
            
            inew_cticks=0;



            iset_min_clim=0;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=0;
%            clim_max=295;
            clim_max=1;
            
            

            
            
            
        case 'Max Low Cloud Fraction no screening GCM'
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0                
                ifilter_clhcalipso=1;
                %These are the thresholds for data that we want to INCLUDE
                %I.e. data outside of these will be discounted
                thresh_clh = [-0.01 0.3];                               
%                thresh_clh = [0.3 1.01];               
            end

            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['cf_isccp_low_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %don't use this anymore - use the time_inds_average2 method for
            %local times as a function of longitude
%            dat_modis = dat_modis(time_inds_average,:,:);
             dat_modis = dat_modis+time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            cont_dat_choose = 'qv';
            cont_dat_choose = 'LTS';
            switch cont_dat_choose
                case 'qv'
                    cont_dat = eval(['1e3*gcm_qv700_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:2:16];
                case 'LTS'
                    cont_dat = eval(['gcm_LTS_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:4:30];
            end
            
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=1;
            
            
        case 'Max Low Cloud Fraction no screening GCM bias from satellite'

            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0                
                ifilter_clhcalipso=1;
                %These are the thresholds for data that we want to INCLUDE
                %I.e. data outside of these will be discounted
                thresh_clh = [-0.01 0.3];                               
%                thresh_clh = [0.3 1.01];   

                   model_data_for_bias ='model_CF';
                   model_data_for_bias ='COSP_CALIPSO_CF';

                   sat_data_for_bias = 'CALIPSO';
%                   sat_data_for_bias = 'MOD35';  %Probably better to use CALISPO since have both day and night data.
                   % Despite the issues in high cloud regions - screen for
                   % mid+high CALIPSO CF regions (separately for day and
                   % night). Or, could do MOD35 for day and CALIPSO for
                   % night for the whole region - if only want to
                   % concentrate on low cloud then doesn't matter.
                   
                   iocean_only=1;
            end

            switch model_data_for_bias
                case 'model_CF'  
                    %            dat_modis = eval(['gcm_cf_' gcm_str ';']);
                    dat_modis = eval(['cf_isccp_low_' gcm_str]);

                case 'COSP_CALIPSO_CF'
                    if cosp_flag==1
                        %Use COSP-CALIPSO for bias calc instead of model CF
                        dat_modis = eval(['cllcalipso_' gcm_str '/100']);
                    else  %make NaN just to give a NaN mean bias etc. for simplicity
                        dat_modis = eval(['NaN*cf_isccp_low_' gcm_str]);
                    end
            end
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %don't use time_inds_average anymore - use the time_inds_average2 method for
            %local times as a function of longitude
%            dat_modis = dat_modis(time_inds_average,:,:);
             dat_modis = dat_modis+time_inds_average2;

            [P2,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' using ' model_data_for_bias];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            cont_dat_choose = 'qv';
            cont_dat_choose = 'LTS';
            switch cont_dat_choose
                case 'qv'
                    cont_dat = eval(['1e3*gcm_qv700_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:2:16];
                case 'LTS'
                    cont_dat = eval(['gcm_LTS_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:4:30];
            end
            
            inew_cticks=0;



            iset_min_clim=0;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=0;
%            clim_max=295;
            clim_max=1;
            
            
            
            %calculate a bias relative to CALIPSO

                interp_to_calipso_grid   %creates the 2D grid for CALIPSO
                gcm_cf_grid = griddata(Plon2D,Plat2D,P2,X,Y); %interpolate GCM data onto same grid
                Plat2D = Y;
                Plon2D = X;

                Plon2D_edges = X2;
                Plat2D_edges = Y2;


                 if min(times_required==[12:15])==1
                     day_or_night = 'day';
                 elseif min(times_required==[0:3])==1
                     day_or_night = 'night';
                 else
                     error('Problem with day or night designation');
                 end
                 
                 switch day_or_night
                     case 'day'
                         cll = meanNoNan(cllcalipso_monthly,1)/100;
                         clh2 = meanNoNan(clhcalipso_monthly,1)/100;
                         clm2 = meanNoNan(clmcalipso_monthly,1)/100;
                         clh = clm2 + clh2 - clm2.*clh2;  %random overlap assumption
                     case 'night'
                         cll = meanNoNan(cllcalipso_monthly_NIGHTTIME,1)/100;
                         clh2 = meanNoNan(clhcalipso_monthly_NIGHTTIME,1)/100;
                         clm2 = meanNoNan(clmcalipso_monthly_NIGHTTIME,1)/100;
                         clh = clm2 + clh2 - clm2.*clh2;  %random overlap assumption
                 end
                 
                                 %remove data where there is significant high level CALIPSO
            %cloud
             if ifilter_clhcalipso==1              
                iclh = find(clh<=thresh_clh(1) | clh>thresh_clh(2));
                gcm_cf_grid(iclh) = NaN;
                
                MODIS_varname2_plot = [MODIS_varname2_plot ' for CALIPSO high CF.GT.' num2str(thresh_clh(1)) '.AND.LTE.' num2str(thresh_clh(2))];
             end
             

             
             switch sat_data_for_bias
                 case 'CALIPSO'
                     P = gcm_cf_grid-cll;
                     MODIS_varname2_plot = [MODIS_varname2_plot ' using CALISPO'];
                 case 'MOD35'
                     dat_modis_MOD35 = Cloud_Fraction_Day_Mean.timeseries3; 
                     
                     modisL3_screening_timeseries3
                     
                     dat_modis_MOD35(ihtot)=NaN;
                     dat_modis_MOD35 = dat_modis_MOD35(:,:,time_inds_average);
                     [P22,Npoints] = meanNoNan(dat_modis_MOD35,3); %time mean (all times)
            
                     modis_cf_grid = griddata(MLON,MLAT,P22,X,Y);
                     
                     P = gcm_cf_grid - modis_cf_grid;
                     
                     MODIS_varname2_plot = [MODIS_varname2_plot ' using MOD35 and for ' thresh_str];

             end
             
             if iocean_only==1
                  %the landmask from CAMCLUBBv2.
                  load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
                 iland=find(gcm_landmask_cal>0.01 | isnan(gcm_landmask_cal)==1 );
                 P(iland)=NaN;

             end

                mean_bias_gcm = meanNoNan(P(:),1);
                RMSE_gcm = sqrt(meanNoNan(P(:).^2,1));

             

            
            
        
        case 'Max Mid Cloud Fraction no screening GCM'
            %this uses the original CF as read from the model output - i.e.
            %no screening for min lwc, ice, etc
            
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['cf_isccp_mid_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %don't use this anymore - use the time_inds_average2 method for
            %local times as a function of longitude
%            dat_modis = dat_modis(time_inds_average,:,:);
             dat_modis = dat_modis+time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=0.3;
        
case 'Max High Cloud Fraction no screening GCM'
                        
%            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            dat_modis = eval(['cf_isccp_high_' gcm_str]);            
            
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %don't use this anymore - use the time_inds_average2 method for
            %local times as a function of longitude
%            dat_modis = dat_modis(time_inds_average,:,:);
             dat_modis = dat_modis+time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=1;
        
            
        case 'Max Cloud Fraction no screening GCM'
            
            dat_modis = eval(['gcm_cf_' gcm_str ';']);
            
            dat_modis = squeeze(max(dat_modis,[],2));
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=1;
            
          case 'LWP GCM'
              %this is the original unchanged LWP as loaded from the output
              %(no screening).
            
            dat_modis = eval(['1e3*gcm_lwp_' gcm_str ';']);
%            dat_modis = 1e3*gcm_lwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=100;
            
            
          case 'IWP to TWP ratio GCM'
            
             iwp = eval(['1e3*gcm_iwp_' gcm_str]);
             lwp = eval(['1e3*gcm_lwp_' gcm_str]);   
             
            twp_thresh = 1;
            dat_modis = iwp./(iwp+lwp);
            dat_modis(iwp+lwp<twp_thresh)=NaN;
%            dat_modis = 1e3*gcm_iwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for TWP GTE ' num2str(twp_thresh) ' g m^{-2} '];
            
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using
            %            normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=0.5;

            

         case 'IWP GCM'
            
            dat_modis = eval(['1e3*gcm_iwp_' gcm_str ';']);
%            dat_modis = 1e3*gcm_iwp_all_clouds;        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' all_clouds field'];
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];


            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=0;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=100;

            
        case 'Cloud Top Temperature GCM'
            
            dat_modis = eval(['gcm_CTT_layer_' gcm_str ';']);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
%            clim_min=270;
             clim_min=280;
            iset_max_clim=1;
%            clim_max=295;
            clim_max=292;

 case 'Cloud Top Temperature GCM'
            
            dat_modis = gcm_CTT_layer;
            %is the temperature of the uppermost continuous layer that is identified
            %although there is screening that is applied before
            %indentifying
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=270;
            iset_max_clim=1;
            clim_max=295;
            
            
            
 case 'Surface Temperature GCM'
            
            dat_modis = eval(['gcm_tsurf_' gcm_str]);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=275;
            clim_min=285;            
            iset_max_clim=1;
            clim_max=305;
            clim_max=294;            
            
            
            case 'Surface Height from Geopotential GCM'
            
            dat_modis = gcm_surfgeo/9.80616;
            %geopotential should be integral of g over height. So assuming
            %g=9.81 at sea-level and that it doesn't vary with height much
            %can prob just divide by this to get height?
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];


            units_str_plot='m';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=4000;
            
            
            case 'Land Mask GCM'
            
            dat_modis = gcm_landmask;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            P = dat_modis(:,:);
            Ndays = 1;

%            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=275;
            iset_max_clim=0;
            clim_max=305;
            
            
            
 case 'Surface Height GCM'
            
            dat_modis = gcm_zsurf;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            P = dat_modis(:,:);
            Ndays = 1;

%            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='m';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=4000;       
            
            
 case 'Lat spacing GCM'
     %CAM5 has regular spacing througout - just plot for AM3
            
            P = diff(gcm_Plat2D_edges,[],1);
            P(end+1,:) = NaN;

%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            clim_min=500;
            clim_min=350;            
            iset_max_clim=0;
            clim_max=1250;
%            clim_max=2500;       

 case 'Lon spacing GCM'
     %CAM5 has regular spacing througout - just plot for AM3
            
            P = diff(gcm_Plon2D_edges,[],2);
            
         

            P(:,end+1) = NaN;

%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;


            iset_min_clim=1;
            clim_min=minALL(P);
%            clim_min=500;
%            clim_min=350;            
            iset_max_clim=1;
            clim_max=maxALL(P);
%            clim_max=2500; 

          

            
            
            
 case 'Cloud Top Height GCM'
            
            dat_modis = eval(['gcm_CTH_' gcm_str]);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='m';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            clim_min=500;
            clim_min=350;            
            iset_max_clim=1;
            clim_max=1250;
%            clim_max=2500;            
            
 case 'Cloud Base Height GCM'
            
            dat_modis = gcm_CBH;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
%            thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
thresh_str='';
%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='m';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=2500;

            
        case 'Cloud layer depth GCM'
            dat_modis = (gcm_CTH-gcm_CBH);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);
thresh_str=[' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];

%            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
%                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];

            units_str_plot='m';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=150;
            iset_max_clim=1;
            clim_max=350;




        case 'Max Cloud Fraction after screening GCM'

            dat_modis = eval(['gcm_CF_max_screened_' gcm_str]);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
%             if icf_low_only==0
%                 thresh_str=[thresh_str 'NO pressure screening - all clouds'];
%             else
%                 thresh_str=[thresh_str 'P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'];
%             end
% 
%             if i_ice_screen==1
%                 thresh_str=[thresh_str ' Liquid only screening'];
%             end
% 
% 
% 
%             thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;
         
        case 'TOA SW Cloud Forcing COSP-CERES GCM'    
            dat_modis = - (toa_SW_outgoing_ceres - toa_SW_outgoing_clearsky_ceres);
            %i.e. the difference between the all-sky (clear+cloudy) and the
            %clear-sky SW at TOA
            %minus sign as the convention must be for the net to mean the
            %net incoming - checks out with other plots that show a
            %negative SWCF
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);  
            P = dat_modis;              

%            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='W m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=1;   
            
         case 'TOA LW Cloud Forcing COSP-CERES GCM'    
            dat_modis = -(toa_LW_outgoing_ceres - toa_LW_outgoing_clearsky_ceres);
            %i.e. the difference between the all-sky (clear+cloudy) and the
            %clear-sky LW at TOA
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);  
            P = dat_modis;              

%            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='W m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=1;     
            
        case 'Cloud-Sat warm rain precip rate, CF normalised'

            dat_modis = rain_warm_ocean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            P = dat_modis(time_inds_average,:,:);             

            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
%            thresh_str=[thresh_str gcm_str ' '];
            
            
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='mm hr^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.2;        
            
            
          
        
            
        case 'Cloud-Sat warm rain precip rate'
            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];
            
            inormalise_precip=0;
            asc_desc = 'descending';            
%            asc_desc = 'ascending';
%            asc_desc = '';            
            
            if inormalise_precip==1    
                normalise_data_get
                
               normalise_data = LWP_monthly_coarse/1000;  %convert LWP to kg/m2              
                
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=2.5;
            
            else
                normalise_data = 1;

                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=0.2;
                
                inew_cticks=1; %keep iset_min_clim as zero if using this
                % see ctick_range_choose for settings.
                
                units_str_plot='mm hr^{-1}'; 
            end
            
            dat_modis = rain_warm_ocean;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
            
            dat_modis = dat_modis ./ normalise_data;
            dat_modis = dat_modis(time_inds_average,:,:);
            
            switch asc_desc
                case 'ascending'                   
                    dat_modis = dat_modis(1:2:end,:,:);
%                    time_inds_average = time_inds_average(1:2:end);
                    iasc_desc=1; %set this in case want to run case 115 of waterVap (to get the right
                    %time inds in time_inds_modisL3_timeseries3
                case 'descending'
                    dat_modis = dat_modis(2:2:end,:,:);    
%                    time_inds_average = time_inds_average(2:2:end);
                    iasc_desc=2;
            end
            
% %set this near the start as needed for time_inds
%             switch iasc_desc
%                 case 1  %'asc'
%                     dat_modis = dat_modis(1:2:end,:,:);
%                 case 2  %'desc'
%                     dat_modis = dat_modis(2:2:end,:,:);
%             end
            



                              

            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str=[asc_desc];
            
%            thresh_str=[thresh_str gcm_str ' '];
            
            
           


            

           
            title_info = '';

            ifilter_ndays=0;

            icontour=0;


   case 'Cloud-Sat warm rain fraction precip rate'
            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];
            
            inormalise_precip=0;
            if inormalise_precip==1                
                normalise_data2 = cllcalipso_monthly(13:60,:,:)/100;
                
                %want to average the lon dimension from 2 to 4 degrees
                %so make an array half the size in lon
                normalise_data3 = NaN*ones([size(normalise_data2,1),size(normalise_data2,2),size(normalise_data2,3)/2]);
                
                for inorm=1:2:size(normalise_data2,3)
                    ind_norm = (inorm)/2+0.5;
                    normalise_data3(:,:,ind_norm) = 0.5*(normalise_data2(:,:,inorm)+normalise_data2(:,:,inorm+1));
                 
                end
                
                
                normalise_data3(normalise_data3<0.01)=NaN;
                
                normalise_data = NaN*ones([size(normalise_data2,1)*2,size(normalise_data2,2),size(normalise_data2,3)/2]);
                
                for inorm=1:size(normalise_data3,1)
                    ind_norm = (inorm-1)*2+1;
                    normalise_data(ind_norm,:,:) = normalise_data3(inorm,:,:);
                    normalise_data(ind_norm+1,:,:) = normalise_data3(inorm,:,:);                    
                end
                

                
                MODIS_varname2_plot = [MODIS_varname2_plot ' normalised by CF'];
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=0.2;
            
            else
                normalise_data = 1;

                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=01;
            end
            
            dat_modis = rain_warm_ocean./rain_ocean;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            thresh_min_precip = 0.01;
            dat_modis(rain_ocean<thresh_min_precip)=NaN;

            P = dat_modis(time_inds_average,:,:);                                  


            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
%            thresh_str=[thresh_str gcm_str ' '];
            
            
           


            

            units_str_plot='mm hr^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



 case 'Cloud-Sat warm rain precip rate BIAS (GCM minus CLOUDSAT)'
           
            
             gcm_str = 'AM3';
%             gcm_str = 'CAM5';     

 MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' '];
            
        
            
            
            
            inormalise_precip=0;
            if inormalise_precip==1                
                normalise_data2 = cllcalipso_monthly(13:60,:,:)/100;
                
                %want to average the lon dimension from 2 to 4 degrees
                %so make an array half the size in lon
                normalise_data3 = NaN*ones([size(normalise_data2,1),size(normalise_data2,2),size(normalise_data2,3)/2]);
                
                for inorm=1:2:size(normalise_data2,3)
                    ind_norm = (inorm)/2+0.5;
                    normalise_data3(:,:,ind_norm) = 0.5*(normalise_data2(:,:,inorm)+normalise_data2(:,:,inorm+1));
                 
                end
                
                
                normalise_data3(normalise_data3<0.01)=NaN;
                
                normalise_data = NaN*ones([size(normalise_data2,1)*2,size(normalise_data2,2),size(normalise_data2,3)/2]);
                
                for inorm=1:size(normalise_data3,1)
                    ind_norm = (inorm-1)*2+1;
                    normalise_data(ind_norm,:,:) = normalise_data3(inorm,:,:);
                    normalise_data(ind_norm+1,:,:) = normalise_data3(inorm,:,:);                    
                end
                

                
                MODIS_varname2_plot = [MODIS_varname2_plot ' normalised by CF'];
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=0.2;
            
            else
                normalise_data = 1;

                iset_min_clim=1;
                clim_min=-0.02; % -0.02
                iset_max_clim=1;
                clim_max=0.02;  % 0.02
            end
            
            normalise_data=1;
            
                
           
            
            

            eval(['dat_modis_GCM = 1e3*3600*gcm_precT_' gcm_str './normalise_data;']);
            
            dat_modis_GCM = meanNoNan(dat_modis_GCM,1);



            dat_modis_cloudsat = rain_warm_ocean;
            dat_modis_cloudsat = meanNoNan(dat_modis_cloudsat,1);
%            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
%            dat_modis(ihtot)=NaN;

%            P = dat_modis(time_inds_average,:,:);          
                         
%            [P,Ndays2] = meanNoNan(dat_modis,1);

            Plat_reg = [-40:1:40];
            Plon_reg=  [-180:1:-60];
            [XI,YI]=meshgrid(Plon_reg,Plat_reg);
            
            dat_modis_GCM_reg = eval(['griddata(gcm_Plon2D_edges_' gcm_str ',gcm_Plat2D_edges_' gcm_str ',dat_modis_GCM,XI,YI)']);
            dat_modis_cloudsat_reg = griddata(Plon2D,Plat2D,dat_modis_cloudsat,XI,YI);
            
            [P] = dat_modis_GCM_reg - dat_modis_cloudsat_reg; 
            
            [Plon2D, Plat2D] = meshgrid(Plon_reg,Plat_reg);
            
            Plat_reg_edges = [Plat_reg-0.5 Plat_reg(end)+0.5];
            Plon_reg_edges = [Plon_reg-0.5 Plon_reg(end)+0.5];
            
            [Plon2D_edges, Plat2D_edges] = meshgrid(Plon_reg_edges,Plat_reg_edges);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
%            thresh_str=[thresh_str gcm_str ' '];
            
            
           


            

            units_str_plot='mm hr^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;

     
            
            
 case 'Cloud-Sat all rain precip rate'
            dat_modis = rain_ocean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            P = dat_modis(time_inds_average,:,:);  

            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
%            thresh_str=[thresh_str gcm_str ' '];
            
            
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='mm hr^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.2;                    
        
            
         case 'Total Cloud Fraction COSP-CALIPSO GCM'


%            dat_modis = low_cf_COSP_CAL/100;
            dat_modis = tot_CF_calipso/100;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);  
            P = dat_modis;              

%            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;         
            
        case 'Low Cloud Fraction CALIPSO'
            calipso_daynight_label= 'DAYTIME';     
            calipso_daynight_label= 'DAYTIME2';                 
%            calipso_daynight_label= 'NIGHTTIME';
%            calipso_daynight_label= 'AVERAGE';            


%            dat_modis = low_cf_COSP_CAL/100;
%            dat_modis = low_CF_calipso/100;

switch calipso_daynight_label
    case 'DAYTIME'
        dat_modis = cllcalipso_monthly/100;  %using the cll* CMOR files - not 100% that they are day only??
                                            %UPDATE - just checked and it is daytime and indentical to daytime2
                                            
    case 'DAYTIME2'
        dat_modis = cllcalipso_monthly_DAYTIME2/100;  %using teh MapLowMidHigh* files that are labelled day      
    case {'NIGHTTIME'}
        dat_modis = cllcalipso_monthly_NIGHTTIME/100;
    case 'AVERAGE' 
        dat_modis = cllcalipso_monthly_AVERAGE/100;
end
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average


%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' calipso_daynight_label ' for ' years_calipso_str];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;   
            
 case 'Combined Cloud Fraction COSP-CALIPSO'
            calipso_daynight_label= 'DAYTIME';            
%            calipso_daynight_label= 'NIGHTTIME';
%            calipso_daynight_label= 'AVERAGE';            


%            dat_modis = low_cf_COSP_CAL/100;
%            dat_modis = low_CF_calipso/100;

switch calipso_daynight_label
    case 'DAYTIME'
        dat_modis = (cllcalipso_monthly+clmcalipso_monthly+clhcalipso_monthly)/100;
    case {'NIGHTTIME'}
        dat_modis = (cllcalipso_monthly_NIGHTTIME+clmcalipso_monthly_NIGHTTIME+clhcalipso_monthly_NIGHTTIME)/100;
    case 'AVERAGE' 
        dat_modis = (cllcalipso_monthly_AVERAGE+clmcalipso_monthly_AVERAGE+clhcalipso_monthly_AVERAGE)/100;
end
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average


%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' calipso_daynight_label ' for ' years_calipso_str];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;   
       
            
    case 'LWP AMSRE'
        
        if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
        
            amsre_daynight_label= 'daytime';            
            amsre_daynight_label= 'nighttime';
%            amsre_daynight_label= 'average';  
%            amsre_daynight_label= 'diurnal range';  
            
            iadd_correct_for_clear_sky_bias = 0; %Whether to correct AMSRE for the daytime clear-sky  bias
             %as calculated using MODIS CF
             
        end
            
            iadd_MODIS_AMSRE_bias_contour=0; %whether to plot the contour where MODIS vs AMSRE biases were high
            
% For the time filtering this searches through all of the requested months and years 
% looking for matches in month_amsre and year_amse.
% I think that month_amsre etc are set up for lwp_amsre and not
% lwp_amsre_time3. So if have filled in the Terra entries using the AMSRE
% LWP then it shouldn't include them here (which is good).

%This needs to be set to tell saveload_MULTI_lon_transect_data whether to
%save as day or night data
am3_dataset = [amsre_daynight_label];


amsre_case = 'normal';

%took out the option to plot time3 here, since don't think it would work
%because of the options set earlier based on sensor type, etc. So use the
%'LWP AMSRE time3' case.
switch amsre_case
    case 'normal'

        switch amsre_daynight_label
            case 'daytime'
                dat_modis = 1000*lwp_amsre(:,:,:,1);
            case {'nighttime'}
                dat_modis = 1000*lwp_amsre(:,:,:,2);
            case 'average'
                dat_modis = 1000*meanNoNan(lwp_amsre,4);
            case 'diurnal range'
                dat_modis = 1000 * (lwp_amsre(:,:,:,2) - lwp_amsre(:,:,:,1) );
        end
%        dat_modis = permute(dat_modis,[2 1 3]); %seems to be the right way
%        if loaded from load_amsre_saved_data
        
   

end

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' amsre_daynight_label];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average


%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3);
            Npoints = Ndays2;
            
            if iadd_correct_for_clear_sky_bias==1
                filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
%                filename = ['AMSRE_clear_sky_bias_5percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to                
                filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                filename_savevars = [filedir_savevars filename '.mat'];
                load(filename_savevars);
                tag = '_AMSRE_clear_sky';
                
                
                
                AMSRE_clear_sky = griddata(Plon2D_AMSRE_clear_sky,Plat2D_AMSRE_clear_sky,P_save_AMSRE_clear_sky,Plon2D,Plat2D); %for AMSRE

                P = P - AMSRE_clear_sky;  %is in g/m2
                %subtract the clear-sky bias 

                MODIS_varname2_plot = [MODIS_varname2_plot ' clear-sky correction ' filename ' '];


            end
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           




            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;    
            
            if iadd_MODIS_AMSRE_bias_contour==1
                        %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=1;
                            %load in the bias data
                            save_file='~/MODIS_vs_AMSRE_LWP_bias_saved.mat';
                            load(save_file,'Psave'); %don't want to load screen_edits
                            P_LWP_bias_allCF_reLT14 = Psave{4}; %Precentage bias                                                                                    
                            cont_dat = abs(P_LWP_bias_allCF_reLT14);
                            
                            %load the grid, but don't overwrite MLAT if
                            %exists in memory
                            grid_AMSRE_bias = load(save_file,'MLAT','MLON');
                            

                            cont_dat =  griddata(grid_AMSRE_bias.MLON,grid_AMSRE_bias.MLAT',cont_dat,Plon2D,Plat2D); 
                            cont_ints=[20 20];
                            
            end
                            

    case 'LWP AMSRE conditional sampling'
        
        if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
        
            amsre_daynight_label= 'daytime';            
            amsre_daynight_label= 'nighttime';
%            amsre_daynight_label= 'average';  
%            amsre_daynight_label= 'diurnal range';  
            
            iadd_correct_for_clear_sky_bias = 0; %Whether to correct AMSRE for the daytime clear-sky  bias
             %as calculated using MODIS CF
             
             LWP_thresh = 50; %ignore any data with LWP less than this (g/m2).
             
        end
            
            iadd_MODIS_AMSRE_bias_contour=0; %whether to plot the contour where MODIS vs AMSRE biases were high
            
% For the time filtering this searches through all of the requested months and years 
% looking for matches in month_amsre and year_amse.
% I think that month_amsre etc are set up for lwp_amsre and not
% lwp_amsre_time3. So if have filled in the Terra entries using the AMSRE
% LWP then it shouldn't include them here (which is good).

%This needs to be set to tell saveload_MULTI_lon_transect_data whether to
%save as day or night data
am3_dataset = [amsre_daynight_label];


amsre_case = 'normal';

%took out the option to plot time3 here, since don't think it would work
%because of the options set earlier based on sensor type, etc. So use the
%'LWP AMSRE time3' case.
switch amsre_case
    case 'normal'

        switch amsre_daynight_label
            case 'daytime'
                dat_modis = 1000*lwp_amsre(:,:,:,1);
            case {'nighttime'}
                dat_modis = 1000*lwp_amsre(:,:,:,2);
            case 'average'
                dat_modis = 1000*meanNoNan(lwp_amsre,4);
            case 'diurnal range'
                dat_modis = 1000 * (lwp_amsre(:,:,:,2) - lwp_amsre(:,:,:,1) );
        end
%        dat_modis = permute(dat_modis,[2 1 3]); %seems to be the right way
%        if loaded from load_amsre_saved_data
        
   

end

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' amsre_daynight_label];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
            

            dat_modis(dat_modis<LWP_thresh) = NaN;

%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3);
            Npoints = Ndays2;
            
            if iadd_correct_for_clear_sky_bias==1
                filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to
                filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                filename_savevars = [filedir_savevars filename '.mat'];
                load(filename_savevars);
                tag = '_AMSRE_clear_sky';
                
                AMSRE_clear_sky = griddata(Plon2D_AMSRE_clear_sky,Plat2D_AMSRE_clear_sky,P_save_AMSRE_clear_sky,Plon2D,Plat2D); %for AMSRE

                P = P - AMSRE_clear_sky;  %is in g/m2
                %subtract the clear-sky bias 

                MODIS_varname2_plot = [MODIS_varname2_plot ' clear-sky correction ' filename ' ONLY for LWP GTE ' num2str(LWP_thresh) 'g m^{-2} '];


            end
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           




            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;    
            
            if iadd_MODIS_AMSRE_bias_contour==1
                        %add the contour of the region where
                            %uncertainties in the Nd obs are large - so
                            %those where either the Bias cf. AMSRE is large
                            %, or those where the change with CF>0.8 and
                            %all CF is large.
                            icontour=1;
                            %load in the bias data
                            save_file='~/MODIS_vs_AMSRE_LWP_bias_saved.mat';
                            load(save_file,'Psave'); %don't want to load screen_edits
                            P_LWP_bias_allCF_reLT14 = Psave{4}; %Precentage bias                                                                                    
                            cont_dat = abs(P_LWP_bias_allCF_reLT14);
                            
                            %load the grid, but don't overwrite MLAT if
                            %exists in memory
                            grid_AMSRE_bias = load(save_file,'MLAT','MLON');
                            

                            cont_dat =  griddata(grid_AMSRE_bias.MLON,grid_AMSRE_bias.MLAT',cont_dat,Plon2D,Plat2D); 
                            cont_ints=[20 20];
                            
            end
                            
                            

            
            
            

     case 'LWP AMSRE time3'
            amsre_daynight_label= 'daytime';            
%            amsre_daynight_label= 'nighttime';
%            amsre_daynight_label= 'average';  


%This needs to be set to tell saveload_MULTI_lon_transect_data whether to
%save as day or night data
am3_dataset = [amsre_daynight_label];

 ifilter_ndays=1;
 





switch amsre_daynight_label
    case 'daytime'
        dat_modis = 1000*lwp_amsre_time3(:,:,:,1);
    case {'nighttime'}
        dat_modis = 1000*lwp_amsre_time3(:,:,:,2);
    case 'average' 
        dat_modis = 1000*meanNoNan(lwp_amsre_time3,4);
end

%dat_modis = permute(dat_modis,[2 1 3]);
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average


%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3);

        
            
            thresh_str=[thresh_str gcm_str ' '];
            
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' amsre_daynight_label ' ' thresh_str];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            

            
            

            icontour=0;
            inew_cticks=0;

            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                ifilter_ndays=1;

                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=150;
            end
            
            
%            if ifilter_ndays==1
%                 MODIS_varname2_plot = [MODIS_varname2_plot ' for ndays>' num2str(thresh_ndays) ' only '];
%            end
            
            
    case 'In-cloud AMSRE LWP (normalised by MOD35 CF)'
            amsre_daynight_label= 'daytime';            
% only can do for daytime as only have daytime MODIS CF 


%This needs to be set to tell saveload_MULTI_lon_transect_data whether to
%save as day or night data
am3_dataset = [amsre_daynight_label];



switch amsre_daynight_label
    case 'daytime'
        dat_modis = 1000*lwp_amsre_time3(:,:,:,1);
%     case {'nighttime'}
%         dat_modis = 1000*lwp_amsre(:,:,:,2);
%     case 'average' 
%         dat_modis = 1000*meanNoNan(lwp_amsre,4);
end

cf = Cloud_Fraction_Day_Mean.timeseries3;
CF_thresh=0.01;
cf(cf<CF_thresh)=NaN;

%normalise by the CF to get the in-cloud LWP
dat_modis = dat_modis./cf;


%dat_modis = permute(dat_modis,[2 1 3]);
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average


%            P = dat_modis;                

            [P,Ndays2] = meanNoNan(dat_modis(:,:,time_inds_average),3);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' amsre_daynight_label];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;             

            
            
         case 'Nd Liquid COSP-MODIS GCM'
             %looks like the COSP tau and re variables are grid-box
             %averages - so need to divide by the CF to get in-cloud values
             %Doesn't say in the nc file, but checked by comparing the
             %means for small CF ranges, e.g. 0.6-0.65 and 0.9-0.95.
             %Without dividing by the CF the ratio of these came out as equal to
             %0.925/0.625 and was =1 when dividing by CF.
             cf_cosp = eval(['liqCF_modis_' gcm_str])/100;            
             CF_gcm_thresh=0.01;
             CF_gcm_thresh=0.8;

             thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
%             cf_cosp(cf_cosp<CF_gcm_thresh)=NaN;
             tau_cosp(cf_cosp<CF_gcm_thresh)=NaN;
             re_cosp(cf_cosp<CF_gcm_thresh)=NaN;

             tau_cosp = eval(['liqTau_modis_' gcm_str]) ./ cf_cosp;
             re_cosp = eval(['liqRe_modis_' gcm_str]) ./ cf_cosp; %(m)
             
             t_top2 = 282; %hardwire this for the moment
            
           Wflag='calc';
           [N,H,W,k,Q,cw]=MODIS_N_H_func(tau_cosp,re_cosp,Wflag,NaN,t_top2);
                                       
            dat_modis = N;
             
             
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            


            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;    
         
             inew_cticks=1; %don't set the clims when using this
        
        
        case 'Optical Depth Liquid COSP-MODIS GCM'
            a = eval(['liqTau_modis_' gcm_str]);
            b = eval(['liqCF_modis_' gcm_str])/100;
            
            CF_gcm_thresh=0.01;
            CF_gcm_thresh=0.8;
            
            thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
            
            a(b<CF_gcm_thresh)=NaN;        
             
            dat_modis = a./b;
             
             
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            


            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=30;    
            
            
           case 'Effective Radius Liquid COSP-MODIS GCM'
            a = eval(['liqRe_modis_' gcm_str])*1e6;
            b = eval(['liqCF_modis_' gcm_str])/100;
            
            CF_gcm_thresh=0.01;
            CF_gcm_thresh=[0.8 1.01];
            
            thresh_str=['liqCF.GTE.' num2str(CF_gcm_thresh) ' '];
            
            a(b<CF_gcm_thresh(1))=NaN;        
            a(b>=CF_gcm_thresh(2))=NaN;                    
             
            dat_modis = a./b;
             
             
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            


            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='\mum';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=20;    
            
     
            
            
            
        case 'Total Liquid Cloud Fraction COSP-MODIS GCM'
            dat_modis = eval(['liqCF_modis_' gcm_str '/100']);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis(time_inds_average,:,:) + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;                 
            
            
     
            
 case 'Low Cloud Fraction COSP-CALIPSO GCM'


%            dat_modis = low_cf_COSP_CAL/100;
%            dat_modis = low_CF_calipso/100;

            dat_modis = eval(['cllcalipso_' gcm_str '/100']);
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            cont_dat_choose = 'qv';
%            cont_dat_choose = 'LTS';
            switch cont_dat_choose
                case 'qv'
                    cont_dat = eval(['1e3*gcm_qv700_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:2:16];
                case 'LTS'
                    cont_dat = eval(['gcm_LTS_' gcm_str]);
                    cont_dat(ihtot)=NaN;
                    cont_dat = cont_dat + time_inds_average2;
                    cont_dat = meanNoNan(cont_dat,1);
                    cont_ints=[0:4:30];
            end
            
            
            inew_cticks=0;

            
            


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;                 
            
 case 'High Cloud Fraction CALIPSO'

     calipso_daynight_label= 'DAYTIME';
%     calipso_daynight_label= 'NIGHTTIME';
%     calipso_daynight_label= 'AVERAGE';

     switch calipso_daynight_label
         case 'DAYTIME'
             dat_modis = clhcalipso_monthly/100;
         case 'NIGHTTIME'
             dat_modis = clhcalipso_monthly_NIGHTTIME/100;
         case 'AVERAGE' 
             dat_modis = clhcalipso_monthly_AVERAGE/100;
     end



            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            P = dat_modis;              

            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
           


             MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' calipso_daynight_label ' for ' years_calipso_str];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.3; 
            clim_max=1.0;             
            
            
 case 'Mid Cloud Fraction COSP-CALIPSO GCM'

%            dat_modis = low_cf_COSP_CAL/100;
%            dat_modis = low_CF_calipso/100;

            dat_modis = eval(['clmcalipso_' gcm_str '/100']);
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.3;                 

            
 case 'High Cloud Fraction COSP-CALIPSO GCM'


%            dat_modis = low_cf_COSP_CAL/100;
%            dat_modis = low_CF_calipso/100;

            dat_modis = eval(['clhcalipso_' gcm_str '/100']);
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %new local time based time screening
            %have taken the approach of a 3D matrix that contains all NaNs
                %except at the required times when it is zero. Then this can be added to
                %dat_modis, which will make the data that we don't want to
                %include NaN, but will leave the rest unaffected
            dat_modis = dat_modis + time_inds_average2;

%            P = dat_modis;                

%            [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);
            [P,Ndays2] = meanNoNan(dat_modis,1);            

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            %P = low_CF_calipso/100;
            
            thresh_str=[thresh_str gcm_str ' '];
           


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str gcm_years_loaded_str ' '];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.3;           
            
            
            
            case 'Mid Cloud Fraction CALIPSO'

                
                calipso_daynight_label= 'DAYTIME';
                calipso_daynight_label= 'NIGHTTIME';
%                calipso_daynight_label= 'AVERAGE';

                 switch calipso_daynight_label
                     case 'DAYTIME'
                         dat_modis = clmcalipso_monthly/100;
                     case 'NIGHTTIME'
                         dat_modis = clmcalipso_monthly_NIGHTTIME/100;
                     case 'AVERAGE'
                         dat_modis = clmcalipso_monthly_AVERAGE/100;
                 end



            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);  
%            P = dat_modis;              
 [P,Ndays2] = meanNoNan(dat_modis(time_inds_average,:,:),1);

%            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
            
            thresh_str=[thresh_str gcm_str ' '];
           


             MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' calipso_daynight_label ' for ' years_calipso_str];

            units_str_plot='';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=0.3;  
            

            
        case 'Total surface precip rate GCM'  
            if eval([ 'exist(''gcm_precT_' gcm_str ''')'])
                iuse_saved_data=1;
            else
                iuse_saved_data=0;
            end
            normalise_precip=0; %whether to normalise by the cloud fraction to get the 
            %in-cloud precip rate
%            gcm_str = 'CAM5';
%            gcm_str = 'AM3';    

            inew_cticks=0;
            
            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];
            
            if normalise_precip==1
                ioverride_normalise_flags=1;
                
                calc_or_load = 'load';
                normalise_case = 'LWP';
                data_type_for_normalise = 'GCM';
        
                normalise_data_get;
                normalise_data(normalise_data<0.01)=NaN;
                normalise_data = normalise_data/1000; %convert to kg/m2
               
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=2.5;
                
                units_str_plot='mm m^3 kg^{-1} hr^{-1}';      
            else
                inew_cticks=1; %keep iset_min_clim as zero if using this
                % see ctick_range_choose for settings.
                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=0.2;
                normalise_data = 1;
                units_str_plot='mm hr^{-1}';      
            end
            

            if iuse_saved_data==1
%                normalise_data=1;
                eval(['dat_modis = 1e3*3600*gcm_precT_' gcm_str './normalise_data;']);
%                dat_modis = normalise_data;
            else
                dat_modis = 1e3*3600*gcm_precT./normalise_data; %convert from m/s of precip to mm/hr
                %do we need to divide by the cloud fraction to get the in-cloud
                %precip?
            end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);
            
            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);


%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];




            title_info = '';

            ifilter_ndays=0;

            icontour=0;


        case 'In-Cloud LWP GCM (normalised by MODIS COSP CF)'        
            cf = eval(['liqCF_modis_' gcm_str])/100;
            CF_gcm_thresh=0.01;
            cf(cf<CF_gcm_thresh)=NaN;

            dat_modis = eval(['1e3*gcm_lwp_' gcm_str './cf']);
         
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;

            
        case 'In-Cloud TLWP GCM (normalised by MODIS COSP CF)'        
            cf = eval(['liqCF_modis_' gcm_str])/100;
            CF_gcm_thresh=0.01;
            cf(cf<CF_gcm_thresh)=NaN;

            dat_modis = eval(['1e3*gcm_TLWP_' gcm_str './cf']);
         
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;

            
            

        case 'LWP GCM (grid-box mean, no CF screening)'
             dat_modis = eval(['1e3*gcm_lwp_' gcm_str]);


         

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;

            
        case 'LWP GCM (grid-box mean, no CF screening) conditional sampling'
               if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                   LWP_thresh = 50; %only look at points where LWP is bigger than this (g/m2)
               end
               
             dat_modis = eval(['1e3*gcm_lwp_' gcm_str]);
        

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);


            dat_modis(dat_modis<LWP_thresh) = NaN;

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' ONLY for LWP.GTE.' num2str(LWP_thresh) ' gm^{-2} '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;
            


        case 'TLWP GCM (grid-box mean, no CF screening) conditional sampling'
             dat_modis = eval(['1e3*gcm_TLWP_' gcm_str]);
             
              if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                   LWP_thresh = 50; %only look at points where LWP is bigger than this (g/m2)
               end
        

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis(dat_modis<LWP_thresh) = NaN;

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' ONLY for TLWP.GTE.' num2str(LWP_thresh) ' gm^{-2} '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;


            
            
            
            
            
        case 'LWP GCM (grid-box mean)'
            dat_modis = eval(['1e3*gcm_lwp_minthreshCF_' gcm_str]);


            %do we need to divide by the cloud fraction to get the in-cloud CF??
            %What are units of Nd? Need to ask Chris G.


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;
            
 case 'RWP GCM (grid-box mean)'
            dat_modis = eval(['1e3*rwp_' gcm_str]);         


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=30;        
            
            
            
        case 'TLWP GCM (grid-box mean)'
            dat_modis = eval(['1e3*gcm_TLWP_' gcm_str]);           


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;    
            
  case 'IWP GCM (grid-box mean)'
            dat_modis = eval(['1e3*gcm_iwp_' gcm_str]);               


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;                
            
            
case 'TWP GCM (grid-box mean)'
            dat_modis = eval(['1e3*gcm_TWP_' gcm_str]);         


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g m^{-2}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;    
            
            
            

case 'RWP ratio GCM (grid-box mean)'
    twp = eval(['1e3*gcm_lwp_' gcm_str]) + eval(['1e3*rwp_' gcm_str]);
            dat_modis = eval(['1e3*rwp_' gcm_str])./twp;   
            
%            dat_modis = twp;


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='(unitless)';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;          

        case 'qv700 GCM'
             dat_modis = eval(['1e3*gcm_qv700_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='g kg^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=16;

          case 'LWP removal timescale GCM'
              %CAM5 precip rate is in m/s. Multiply by 1e3 to get mm/s,
              %which is equivalent to kg/m2/s
               LWP = eval(['gcm_lwp_' gcm_str]);  %kg/m2
            precip = eval(['1e3*gcm_precT_' gcm_str]); %kg/m2/s
            %only calcuate when precip greater than a threshold - otherwise
            %the zero precip times will swamp the answer
            %Perhaps better plotted vs precip rate? Or a PDF might be
            %better
            precip_thresh = 1e-5;
            precip(precip<precip_thresh)=NaN;
            
            dat_modis = LWP./precip/3600;
               
         

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='hrs';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=6;
            
            
            
        case 'LTS GCM'
             dat_modis = eval(['gcm_LTS_' gcm_str]);


         

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=30;

       case 'LTS bias from ECMWF'
           %run interpolate_to_2deg first - or actually this script does
           %the LTS bias plot too, so just run that.
             P = LTS_am3_2 - LTS_ec_2;   

          %  dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

         %   dat_modis = dat_modis+time_inds_average2;
            


            
%            [P,Ndays2] = meanNoNan(dat_modis,1);            



            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];
%             MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' num2str(hour_ec) ' UTC ' gcm_years_loaded_str ' '];            

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=-3;
            iset_max_clim=1;
            clim_max=3;            
            
            
  case 'LTS ECMWF'
             dat_modis = eval(['ecLTS_mon_Man_ALL']);        

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis+time_inds_average2;
            
%            dat_modis = dat_modis(time_inds_average,:,:);

% hour_ec = [0 6 12 18];
% clear ihour_ec
% for ihr=1:length(hour_ec)
%     ihour_ec(ihr,:) = find(hours_ec==hour_ec(ihr));
% end

%probably doesn't make sense to plot for individual hours as local time
%will vary - need to sort this out

%month_ec = 

%            dat_modis = dat_modis + time_inds_average2;

%mean over all months for now
%            [P,Ndays2] = meanNoNan(meanNoNan(dat_modis([ihour_ec(:)],:,:,:),1),1);
            
            [P,Ndays2] = meanNoNan(dat_modis,1);            

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];
%             MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' num2str(hour_ec) ' UTC ' gcm_years_loaded_str ' '];            

            units_str_plot='K';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=30;       
            
            
            
 case 'qv700 ECMWF'
             dat_modis = eval(['squeeze(1e3*ecqv_mon_Man_ALL(:,3,:,:))']);        

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

             dat_modis = dat_modis+time_inds_average2;

% hour_ec = [0 6 12 18];
% clear ihour_ec
% for ihr=1:length(hour_ec)
%     ihour_ec(ihr,:) = find(hours_ec==hour_ec(ihr));
% end

%probably doesn't make sense to plot for individual hours as local time
%will vary - need to sort this out

%month_ec = 

%            dat_modis = dat_modis + time_inds_average2;

%mean over all months for now
%            [P,Ndays2] = meanNoNan(meanNoNan(dat_modis([ihour_ec(:)],:,:,:),1),1);
            
            
            
            [P,Ndays2] = meanNoNan(dat_modis,1);

%             if i_ice_screen==0
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2)) '.AND.P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.'...
%                     num2str(thresh_P(2)/100) ' hPa icf__low_only=' num2str(icf_low_only)];
%             else
%                 thresh_str=['min column CF.GT.' num2str(low_cf_thresh(1)) '.AND.LTE.' num2str(low_cf_thresh(2))...
%                     ',icf__low_only=' num2str(icf_low_only) ' ICE screening'];
%             end

%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];
%             MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' num2str(hour_ec) ' UTC ' gcm_years_loaded_str ' '];            

            units_str_plot='g kg^{-1}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=16;       
            

            

        case 'Max Nd in all non-screened layers GCM'
            %             dat_modis = gcm_Nd_maxALL;
            dat_modis = eval(['gcm_Nd_max_screen_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

cosp_screen_str = '';

if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
            screen_using_cosp = 're'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF>0.1'; %screen based on whether have a COSP re value and COSP CFliq>10%            
            screen_using_cosp = 'none'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF range';
end

            switch screen_using_cosp
                case 're'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20) )']);
                    cosp_screen_str = ' only when have COSP re value';
                    %This is also equivalent to have a minimum in-cloud COSP
                    %Tau of 0.3 (since COSP only produces tau and re for
                    %tau>0.3).
                    dat_modis(iremove) = NaN;
                case 're CF>0.1'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>=10) )']);
                    cosp_screen_str = ' only when have COSP re value and COSP CFliq.GT.0.1';   
                    dat_modis(iremove) = NaN;
                case 're CF range' %actually, don't need the extra re screening as once CFliq screening done makes no difference
                    iremove = eval(['find( ~(  liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>thresh_CF(1)*100 & liqCF_modis_' gcm_str '<=thresh_CF(2)*100  ) )']);
                    cosp_screen_str = [' only when have COSP re value and COSP CFliq.GTE.' num2str(thresh_CF(1)) ' AND.LT.' num2str(thresh_CF(2))];   
                    dat_modis(iremove) = NaN;    
            end

% *** IMPORTANT technicality on screening :-
%   The method used here is to find set up the screen-str to represent all of the points that lie within the
%   range of the bounds. 
%   And then to do ihtot = find( ~eval(screen_str) ) to find the points to be removed.
%   I.e. all of the remaining points out of the set that satisfy
%   screen_str. 
%   Importantly, this will also include any locations where teh screening
%   data was NaN since NaN data will not satisfy the screen_str logical and
%   so will be included using the NOT. If just finding and removing the locations
%   OUTSIDE of the screening bounds then it will include locations where
%   the screening data is NaN, which we do not want to do.


            
%            dat_modis = dat_modis(time_inds_average,:,:);
            dat_modis = dat_modis + time_inds_average2;
            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [' AWNC ' thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];
            thresh_str = [''];            

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
            
        case 'Max Nd in low cloud with cloud screening'
            %             dat_modis = gcm_Nd_maxALL;
            dat_modis = eval(['gcm_Nd_max_screen_isccp_low_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            cosp_screen_str = '';

if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
            screen_using_cosp = 're'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF>0.1'; %screen based on whether have a COSP re value and COSP CFliq>10%            
            screen_using_cosp = 'none'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF range';
end

            switch screen_using_cosp
                case 're'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20) )']);
                    cosp_screen_str = ' only when have COSP re value';
                    %This is also equivalent to have a minimum in-cloud COSP
                    %Tau of 0.3 (since COSP only produces tau and re for
                    %tau>0.3).
                    dat_modis(iremove) = NaN;
                case 're CF>0.1'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>=10) )']);
                    cosp_screen_str = ' only when have COSP re value and COSP CFliq.GT.0.1';   
                    dat_modis(iremove) = NaN;
                case 're CF range' %actually, don't need the extra re screening as once CFliq screening done makes no difference
                    iremove = eval(['find( ~(  liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>thresh_CF(1)*100 & liqCF_modis_' gcm_str '<=thresh_CF(2)*100  ) )']);
                    cosp_screen_str = [' only when have COSP re value and COSP CFliq.GTE.' num2str(thresh_CF(1)) ' AND.LT.' num2str(thresh_CF(2))];   
                    dat_modis(iremove) = NaN;    
            end

% *** IMPORTANT technicality on screening :-
%   The method used here is to find set up the screen-str to represent all of the points that lie within the
%   range of the bounds. 
%   And then to do ihtot = find( ~eval(screen_str) ) to find the points to be removed.
%   I.e. all of the remaining points out of the set that satisfy
%   screen_str. 
%   Importantly, this will also include any locations where teh screening
%   data was NaN since NaN data will not satisfy the screen_str logical and
%   so will be included using the NOT. If just finding and removing the locations
%   OUTSIDE of the screening bounds then it will include locations where
%   the screening data is NaN, which we do not want to do.


%            dat_modis = dat_modis(time_inds_average,:,:);
            dat_modis = dat_modis + time_inds_average2;
            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [' AWNC ' thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];
            thresh_str = [''];            

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
            
 case 'Max Nd in low cloud without cloud screening'
            %             dat_modis = gcm_Nd_maxALL;
            dat_modis = eval(['gcm_Nd_max_noCF_isccp_low_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

%            dat_modis = dat_modis(time_inds_average,:,:);

cosp_screen_str = '';

if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
            screen_using_cosp = 're'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF>0.1'; %screen based on whether have a COSP re value and COSP CFliq>10%            
            screen_using_cosp = 'none'; %screen based on whether have a COSP re value.
            screen_using_cosp = 're CF range';
end

            switch screen_using_cosp
                case 're'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20) )']);
                    cosp_screen_str = ' only when have COSP re value';
                    %This is also equivalent to have a minimum in-cloud COSP
                    %Tau of 0.3 (since COSP only produces tau and re for
                    %tau>0.3).
                    dat_modis(iremove) = NaN;
                case 're CF>0.1'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>=10) )']);
                    cosp_screen_str = ' only when have COSP re value and COSP CFliq.GT.0.1';   
                    dat_modis(iremove) = NaN;
                case 're CF range' %actually, don't need the extra re screening as once CFliq screening done makes no difference
                    iremove = eval(['find( ~(  liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>thresh_CF(1)*100 & liqCF_modis_' gcm_str '<=thresh_CF(2)*100  ) )']);
                    cosp_screen_str = [' only when have COSP re value and COSP CFliq.GTE.' num2str(thresh_CF(1)) ' AND.LT.' num2str(thresh_CF(2))];   
                    dat_modis(iremove) = NaN;    
            end

% *** IMPORTANT technicality on screening :-
%   The method used here is to find set up the screen-str to represent all of the points that lie within the
%   range of the bounds. 
%   And then to do ihtot = find( ~eval(screen_str) ) to find the points to be removed.
%   I.e. all of the remaining points out of the set that satisfy
%   screen_str. 
%   Importantly, this will also include any locations where teh screening
%   data was NaN since NaN data will not satisfy the screen_str logical and
%   so will be included using the NOT. If just finding and removing the locations
%   OUTSIDE of the screening bounds then it will include locations where
%   the screening data is NaN, which we do not want to do.


            dat_modis = dat_modis + time_inds_average2;
            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [' AWNC ' thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];
            thresh_str = [''];            

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' ' cosp_screen_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
            
            
case 'Nd at max LWC within Cloud Layer GCM'
            dat_modis = eval(['gcm_Nd_max_lwc_cont_' gcm_str]); %Nd at position of max LWC within continuous layers

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
            
            

        case 'Max Nd in Cloud Layer GCM'
            dat_modis = eval(['gcm_Nd_maxlayer_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
        case 'Max Nd COSP MODIS CF screening GCM' 
            dat_modis = Nd;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' ' cf_type ' CF.GT.' num2str(cf_thresh(1)) '.AND.CF.LTE.' num2str(cf_thresh(2))];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;

 case 'LWP COSP MODIS CF screening GCM' 
            dat_modis = eval(['1e3*gcm_lwp_COSPCF_' gcm_str ';']);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' ' cf_type ' CF.GT.' num2str(cf_thresh(1)) '.AND.CF.LTE.' num2str(cf_thresh(2))];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;
            
            
        case 'LWP COSP tau re screening GCM'
            dat_modis = eval(['1e3*LWP_COSP_' gcm_str ';']);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' ' cf_type ' CF.GT.' num2str(cf_thresh(1)) '.AND.CF.LTE.' num2str(cf_thresh(2))];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;            

            
            
   case 'Max Nd noCF minus Max Nd in Cloud Layer GCM' 
%run gcm_process_lite to get Nd
            dat_modis = eval(['gcm_Nd_max_noCF_' gcm_str ' - gcm_Nd_maxlayer_' gcm_str ';']);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis + time_inds_average2;

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;

case 'Max Nd, no CF screening minus CF screening'
%run gcm_process_lite to get Nd, or load from load_gcm_processed_data

rel_or_abs = 'abs';
rel_or_abs = 'rel';

            dat_CF = eval(['gcm_Nd_max_screen_' gcm_str ';']);
            dat_CF(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            %            dat_modis = dat_modis(time_inds_average,:,:);
            dat_CF = dat_CF + time_inds_average2; %new time screening method
            
            daily_or_mean = 'mean';
            
switch daily_or_mean
    case 'daily'            
            %doing the difference on individual days
            dat_modis = eval(['gcm_Nd_max_noCF_' gcm_str ' - dat_CF;']);
    case 'mean'
            %difference on the time means (as done for MODIS comparison)
            [dat_CF,Ndays2] = meanNoNan(dat_CF,1);
            dat_modis = eval(['meanNoNan(gcm_Nd_max_noCF_' gcm_str ',1) - dat_CF;']);
end


            
            units_str_plot='cm^{-3}'; 
            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];


            switch rel_or_abs
                case 'rel'
                    dat_modis = 100 * dat_modis ./ dat_CF;
                    units_str_plot='%'; 
                    MODIS_varname2_plot = ['Percentage diff ' MODIS_varname2_plot];
            end
            

            switch daily_or_mean
                case 'daily'
                    [P,Ndays2] = meanNoNan(dat_modis,1);
                case 'mean'
                    P = dat_modis;
            end




      
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;
            
            cont_dat = P;
            cont_ints=[0 0];

            
            
 case 'Max Nd noCF GCM' 
%run gcm_process_lite to get Nd
            dat_modis = eval(['gcm_Nd_max_noCF_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            cosp_screen_str = '';

           if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
            screen_using_cosp = 're'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF>0.1'; %screen based on whether have a COSP re value and COSP CFliq>10%            
            screen_using_cosp = 'none'; %screen based on whether have a COSP re value.
%            screen_using_cosp = 're CF range';
end

            switch screen_using_cosp
                case 're'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20) )']);
                    cosp_screen_str = ' only when have COSP re value';
                    %This is also equivalent to have a minimum in-cloud COSP
                    %Tau of 0.3 (since COSP only produces tau and re for
                    %tau>0.3).
                    dat_modis(iremove) = NaN;
                case 're CF>0.1'
                    iremove = eval(['find(~ (liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>=10) )']);
                    cosp_screen_str = ' only when have COSP re value and COSP CFliq.GT.0.1';   
                    dat_modis(iremove) = NaN;
                case 're CF range' %actually, don't need the extra re screening as once CFliq screening done makes no difference
                    iremove = eval(['find( ~(  liqRe_modis_' gcm_str '>1e-20 & liqCF_modis_' gcm_str '>thresh_CF(1)*100 & liqCF_modis_' gcm_str '<=thresh_CF(2)*100  ) )']);
                    cosp_screen_str = [' only when have COSP re value and COSP CFliq.GTE.' num2str(thresh_CF(1)) ' AND.LT.' num2str(thresh_CF(2))];   
                    dat_modis(iremove) = NaN;    
            end

% *** IMPORTANT technicality on screening :-
%   The method used here is to find set up the screen-str to represent all of the points that lie within the
%   range of the bounds. 
%   And then to do ihtot = find( ~eval(screen_str) ) to find the points to be removed.
%   I.e. all of the remaining points out of the set that satisfy
%   screen_str. 
%   Importantly, this will also include any locations where teh screening
%   data was NaN since NaN data will not satisfy the screen_str logical and
%   so will be included using the NOT. If just finding and removing the locations
%   OUTSIDE of the screening bounds then it will include locations where
%   the screening data is NaN, which we do not want to do.


            
%            dat_modis = dat_modis(time_inds_average,:,:);
            dat_modis = dat_modis + time_inds_average2;
            
            

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=50;            

            
        case 'Mean Nd in Cloud Layer GCM'
            dat_modis = eval(['gcm_Nd_meanlayer_' gcm_str]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


%            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='cm^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;



            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;            



        case 'Max LWC in Cloud Layer GCM'
            
%            dat_modis = gcm_lwc_minthreshCF./gcm_cf; %LWC here is only if there was no ice in the
            %layer
%            dat_modis = 1e3*squeeze(max(dat_modis,[],2));
%            dat_modis(gcm_cf<0.001)=NaN;

             dat_modis = eval(['gcm_LWC_max_' gcm_str]);


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' ' gcm_str ' '];

            units_str_plot='g m^{-3}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=0.1;
            iset_max_clim=1;
            clim_max=0.6;


        case 'Cloud Fraction at max LWC GCM'
            %            dat_modis = GCM2_Nd;  %gcm_Nd_maxliq; %
            %            dat_modis = gcm_Nd_maxliq;  %./gcm_CF_maxliq;
            dat_modis = squeeze(gcm_CF_maxliq(time_inds_average,:,:));
            %N.B. this CF uses the values screeened for CF - it therefore
            %will just reflect that screening



            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Ndays2]=meanNoNan(dat_modis,1);

            MODIS_varname2_plot = [modis_data_plot time_mean_str];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;
            colormap_choose=jet; %default


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;

        case 'Number of droplets at max LWC GCM'
            %            dat_modis = GCM2_Nd;  %gcm_Nd_maxliq; %

%             if length(size(gcm_Nd_maxliq))==3
%                 dat_modis = gcm_Nd_maxliq(time_inds_average,:,:);
%                 %            dat_modis = gcm_Nd_meanlayer(time_inds_average,:,:);
%             else
%                 dat_modis = gcm_Nd_maxliq;  %./gcm_CF_maxliq;
%             end
            %            dat_modis = gcm_CF_maxliq;
            
            dat_modis = eval(['gcm_Nd_maxliq_' gcm_str]);  

            %do we need to divide by the cloud fraction to get the in-cloud CF??
            %What are units of Nd? Need to ask Chris G.


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

%            if length(size(gcm_Nd_maxliq))==3
                [P,Ndays2] = meanNoNan(dat_modis,1);
%            else
%                P=dat_modis;
%            end


            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot=' cm^{-3}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=1;
            colormap_choose=jet; %default


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;

        case 'Number of droplets vertical mean over LWC GCM'
            %            dat_modis = GCM2_Nd;  %gcm_Nd_maxliq; %

            if length(size(gcm_Nd_maxliq))==3
                %            dat_modis = gcm_Nd_maxliq(time_inds_average,:,:);
                dat_modis = gcm_Nd_meanlayer(time_inds_average,:,:);
            else
                dat_modis = gcm_Nd_maxliq;  %./gcm_CF_maxliq;
            end
            %            dat_modis = gcm_CF_maxliq;

            %do we need to divide by the cloud fraction to get the in-cloud CF??
            %What are units of Nd? Need to ask Chris G.


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            if length(size(gcm_Nd_maxliq))==3
                [P,Ndays2] = meanNoNan(dat_modis,1);
            else
                P=dat_modis;
            end


            thresh_str = [thresh_str ' for nlevs.GT.' num2str(thresh_Nlevs)];

            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot=' cm^{-3}';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=1;
            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;
            
        case 'Surface Pressure GCM'


            dat_modis = gcm_ps/100;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            dat_modis = dat_modis(time_inds_average,:,:);

            [P,Ndays2] = meanNoNan(dat_modis,1);

            %            thresh_str=['P.GTE.' num2str(thresh_P(1)/100) '.AND.LT.' num2str(thresh_P(2)/100) ' hPa'...
            %                ' for lwc.GT.' num2str(low_liq_thresh*1000) ' g kg^{-1} AND Nd.GTE.' num2str(thresh_Nd) ' cm^{-3}'];
            thresh_str='';
                                

            MODIS_varname2_plot = [modis_data_plot ' ' gcm_str ' ' time_mean_str ' ' gcm_years_loaded_str ' '];

            units_str_plot='hPa';      %no units if using normalised std dev
            %            units_str_plot='mg^{-1}';      %no units if using normalised std dev
            title_info = '';

            ifilter_ndays=0;

            icontour=0;
            inew_cticks=0;



            iset_min_clim=1;
            clim_min=600;
            iset_max_clim=1;
            clim_max=1020;



        case 'Liquid Water Path'
            dat_modis = lwp; %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P=dat_modis(row_L2_inds,col_L2_inds);


            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=500;

        case 'Brightness temp band 31'
            dat_modis = BT(:,:,2); %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P=dat_modis;

            h = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Brightness temp diff 29-31'
            dat_modis = BTD(:,:,1); %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Brightness temp diff 31-32'
            dat_modis = BTD(:,:,2); %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Various 5km flags L2 swath'

            title_info = 'Thin cirrus flag ';
            title_info = 'High cloud flag ';





            switch title_info

                case 'Thin cirrus flag '
                    dat_modis = squeeze(qapp_5km(1,:,:));
                case 'High cloud flag '
                    dat_modis = squeeze(qapp_5km(2,:,:));
            end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev



        case 'Sensor zenith angle'
            dat_modis = sensor_zenith; %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';


            %%%%% **** mock L3 plots ****    %%%%%%%

        case 'Number of rejected pixels mock L3'
            dat_modis = Nreject_mockL3; %


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Pixel rejection ratio mock L3'
            dat_modis = Nreject_mockL3./(Np_mockL3+Nreject_mockL3); %dat_modis(Np_mockL3==0)=NaN;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Number of ice or undetermined pixels mock L3'
            %                        dat_modis = Nice_mockL3./Np_mockL3; %
            dat_modis = Nice_mockL3; %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Nd distribution skewness mock L3'
            dat_modis = skewness(Nd_all,0,3); %2nd input is whether to adjust for bias
            % 0 = adjust, 1 = don't adjust

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'N datapoints mock L3'
            %dat_modis = meanNoNan(F_all,3);
            dat_modis = Np_mockL3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Mean Nd (cm^{-3}) mock L3'
            %dat_modis = meanNoNan(F_all,3);
            dat_modis = meanNoNan(Nd_all,3);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Mean Nd from mean Tau & Re (cm^{-3}) mock L3'

            Wflag='calc';
            %[dat_modis,H,W,k,Q,cw]=MODIS_N_H_func(meanTau_mockL3,meanRe_mockL3*1e-6,Wflag,NaN);
            dat_modis = N_meanTauReff;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Mean Re (\mum) mock L3'

            Re_case = '1.6 \mum';
            Re_case = '2.1 \mum';
            Re_case = '3.7 \mum';

            switch Re_case
                case '1.6 \mum'
                    dat_modis=meanRe_mockL3_16;
                case '2.1 \mum'
                    dat_modis=meanRe_mockL3;
                case '3.7 \mum'
                    dat_modis=meanRe_mockL3_37;
            end


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot ' ' Re_case];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=35;

        case 'Mean Tau mock L3'


            dat_modis=meanTau_mockL3;


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=35;


        case 'Mean Nd uncertainty (%) mock L3'
            dat_modis=meanNd_un_mockL3;
            %            dat_modis = sqrt( (0.5*meanTau_un_mockL3).^2 +
            %            (-5/2*meanRe_un_mockL3).^2 ); %very similar over 5x5 pixels -
            %            even so for over 1x1 degree, which was surprising (MPACE
            %            results)

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';



            %%%% ***  end of mock L3 plots *** %%%%
  case 'Scan Time'
            
            swath_case = 'L2';            
%            swath_case = 'L2 joint';
            switch swath_case
                case 'L2'
                    dat_modis = NaN*ones(size(re_un));
                    dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(scantime,5);
                    dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

                    dat_modis = dat_modis(row_L2_inds,col_L2_inds);
            
                case 'L2 joint'
                     dat_modis = solarZA_L2;
                     dat_modis(ihtot)=NaN;
            end
            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            %                        P=dat_modis(row_L2_inds,col_L2_inds);
            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            

        case 'Solar Zenith Angle'
            %doing this first stage to make the size the same
            %as the other arrays for use of ihtot to set values
            %to NaN
            
            swath_case = 'L2';            
%            swath_case = 'L2 joint';
            switch swath_case
                case 'L2'
            dat_modis = NaN*ones(size(re_un));
            dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(solar_zenith,5);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            dat_modis = dat_modis(row_L2_inds,col_L2_inds);
            
                case 'L2 joint'
                     dat_modis = solarZA_L2;
                     dat_modis(ihtot)=NaN;
            end
            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            %                        P=dat_modis(row_L2_inds,col_L2_inds);
            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Sensor Zenith Angle'
            %doing this first stage to make the size the same
            %as the other arrays for use of ihtot to set values
            %to NaN
            dat_modis = NaN*ones(size(re_un));
            dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(sensor_zenith,5);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            dat_modis = dat_modis(row_L2_inds,col_L2_inds);
            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            %                        P=dat_modis(row_L2_inds,col_L2_inds);
            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Surface Temperature L2 swath'
            dat_modis = NaN*ones(size(re_un));
            dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(surface_temp,5);
            %                        dat_modis = split_matrix_int(surface_temp,5);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            dat_modis = dat_modis(row_L2_inds,col_L2_inds);
            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            %                        P=dat_modis(row_L2_inds,col_L2_inds);
            P=dat_modis;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Cloud Top Temperature L2 swath'
            %             dat_modis = NaN*ones(size(re_un));
            %             dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(t_top,5);

            ctt_type = 'old?';
            ctt_type = 'full swath';
            ctt_type = 't_top3';
            
            switch ctt_type
                case 'old?'
            dat_modis = NaN*ones(size(re_un));
            temp = split_matrix(t_top,5);
            %split matrix takes the 5 km values and replicates each into
            %5x5 1 km tiles - the start of these will correspond to row and
            %column 1 in the re_un sized swath. the end values won't be
            %replicated as this would be bigger than the swath - although
            %they could be I guess. Using row_L2_inds will just ignore
            %these end values
            dat_modis(row_L2_inds,col_L2_inds) = temp(row_L2_inds,col_L2_inds);
            %need an re-un sized array for ihtot
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %
            %             dat_modis = split_matrix(t_top,5);
            %             dat_modis = dat_modis(row_L2_inds,col_L2_inds);

        

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


                case 'full swath'
            % ***     set ifull_swath=1   ***
            dat_modis = t_top3;
            dat_modis(ihtot) = NaN;
             P=dat_modis;
            
                case 't_top3'
             % ***     set ifull_swath=0   ***
             dat_modis = t_top3;
             dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
             P=dat_modis(row_L2_inds,col_L2_inds);
                    
                    
            end
            

            %                        P=dat_modis(row_L2_inds,col_L2_inds);
%            P=dat_modis(row_L2_inds,col_L2_inds);
            

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
%            iset_min_clim=1;
%            clim_min=230;
%            iset_max_clim=1;
%            clim_max=280;

       case 'Cloud Top Temperature Joint L2 swath'
           dat_modis = t_top;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='K';      %no units if using normalised std dev
            title_info = '';


        case 'Brightness Temperature L2 swath'
            %             dat_modis = NaN*ones(size(re_un));
            %             dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(squeeze(BT(:,:,2)),5);


            dat_modis = NaN*ones(size(re_un));
            temp = split_matrix(squeeze(BT(:,:,2)),5);
            dat_modis(row_L2_inds,col_L2_inds) = temp(row_L2_inds,col_L2_inds);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            %                        P=dat_modis(row_L2_inds,col_L2_inds);
            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '11 \mum ';

        case 'Bootstrap droplet number error of mean'
            err=(boot_out(:,:,12)-boot_out(:,:,2))/2;

            dat_modis=N_time3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the dataset


            [Nmean,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P = 100* err./Nmean;

            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='%';
            title_info = '';


            inew_cticks=0;



        case 'Cloud Fraction L2 swath'
            dat_modis = cf;
            dat_modis(ihtot)=NaN;
            P=dat_modis(row_L2_inds,col_L2_inds);
            
%            dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(cf,5);
%            dat_modis(ihtot)=NaN;
            
%            dat_modis = NaN*ones(size(re_un));
%            dat_modis(row_L2_inds,col_L2_inds) = split_matrix_int(cf,5);
%            dat_modis(ihtot)=NaN;

            %            dat_modis = cfL2;
            %            dat_modis(ihtot)=NaN; %make these values NaN as
            %            they will then be removed from the average
            %            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';






        case 'Various flags L2 swath'
            title_info = 'Snow-Ice flag ';
                                    title_info = 'Surface-type flag ';
%                                    title_info = 'Heavy aerosol flag ';
%            title_info = 'Thin cirrus flag ';
            %                        title_info = 'Shadow flag ';
                                    title_info = 'Multi-layer flag ';
%            title_info = 'Phase flag primary';
                                    %title_info = 'Optical Thickness out of bounds flag'; %0=in-bounds, 1=100-150, 2=>150
%                        title_info = 'Sunglint flag ';
                                    %title_info = 'Optical Thickness Usefulness flag'; %0=not useful, 1=useful
                                    title_info = 'Optical Thickness Confidence flag'; %0=no confidence, 1=marginal, 2=good, 3=v. good                                  
                                    %title_info = 'Cloud Effective Radius Usefulness flag'; %0=not useful, 1=useful
                                    %title_info = 'Cloud Effective Radius Confidence flag'; %0=no confidence, 1=marginal, 2=good, 3=v. good                                 


            switch title_info
               
                case 'Optical Thickness Usefulness flag'
                    dat_modis = squeeze(qapq_1km(1,:,:));
                case 'Optical Thickness Confidence flag'
                    dat_modis = squeeze(qapq_1km(2,:,:));
                case 'Optical Thickness out of bounds flag'
                    dat_modis = squeeze(qapq_1km(3,:,:));                    
                case 'Cloud Effective Radius Usefulness flag'
                    dat_modis = squeeze(qapq_1km(4,:,:));                    
                case 'Cloud Effective Radius Confidence flag'
                    dat_modis = squeeze(qapq_1km(5,:,:));
                case 'Phase flag primary'
                    dat_modis = squeeze(qapp_1km(1,:,:));
                case 'Snow-Ice flag '
                    dat_modis = squeeze(mask_1km(5,:,:));
                case 'Surface-type flag '
                    dat_modis = squeeze(mask_1km(6,:,:));
                case 'Heavy aerosol flag '
                    dat_modis = squeeze(mask_1km(7,:,:));
                case 'Thin cirrus flag '
                    dat_modis = squeeze(mask_1km(8,:,:));
                case 'Shadow flag '
                    dat_modis = squeeze(mask_1km(9,:,:));
                case 'Multi-layer flag '
                    dat_modis = squeeze(qapp_1km(16,:,:));
                 case 'Sunglint flag '
                    dat_modis = squeeze(mask_1km(4,:,:));     %0=sunlgint, 1= NO sunglint
            end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot ' ' title_info];

            units_str_plot='';      %no units if using normalised std dev
            



        case 'Phase flag L2 swath 1621'
            dat_modis = squeeze(qapq_1km(8,:,:));

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'N_d % Uncertainty L2 swath'

            dat_modis = percent_error_Nd;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'N_d Uncertainty (cm^{-3}) L2 swath'
            %absolute uncertainty
            dat_modis = N .* percent_error_Nd/100;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Optical Depth Uncertainty L2 swath'
            dat_modis = tau_un_abs;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Optical Depth % Uncertainty L2 swath'
            dat_modis = tau_un;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';


        case 'Reff Uncertainty (\mum) L2 swath'
            dat_modis = re_un_abs;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Reff % Uncertainty L2 swath'
            dat_modis = re_un; %re_un is in %

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Reff Difference (\mum) L2 swath'
            Re_diff_case = '1.6 minus 2.1';
%            Re_diff_case = '3.7 minus 2.1';

            switch Re_diff_case
                case '1.6 minus 2.1'
                    dat_modis = squeeze(re_diff(:,:,1));
                case '3.7 minus 2.1'
                    dat_modis = squeeze(re_diff(:,:,2));
            end

            %Cloud Particle Effective Radius two-channel retrieval using band 6 and
            %band 20 differenced from band 7 retrieval and either band 1, 2, or 5
            %(specified in Quality_Assurance_1km).
            %Plane 1 of SDS is band 6 - band 7 retrieval (Re_1.6 - Re_2.1)
            %Plane 2 is band 20 - band 7 retrieval.      (Re_3.7 - Re_2.1)
            %size(re_diff) = [1354        2030           2]
            %band 6 = 1.6um, band 7 = 2.1um (the usual one), band 20 = 3.7um
            %band 1 = 0.62, band 2 = 0.84, band 5 = 1.2 um - don't think we need to
            %worry about whether is band 1, 2 or 5?

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot Re_diff_case];

            units_str_plot='';      %no units if using normalised std dev

        case 'Reff (\mum) L2 swath'
            dat_modis = re;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            iset_min_clim=1;
            clim_min=6;
            iset_max_clim=1;
            clim_max=25;

        case 'Cloud Depth (m) L2 swath'
            dat_modis = H;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Cloud Depth (m) L2 swath mock L3'
            dat_modis = H_meanTauReff;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %
            title_info = '';

        case 'Cloud Fraction L2 swath mock L3'
            dat_modis = CF_mockL3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %
            title_info = '';

        case 'LWP homogeneity factor (Wood 2006) L2 swath mock L3'
            dat_modis = (meanW_mockL3 ./ stdW_mockL3).^2;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %
            title_info = '';

        case 'LWP homogeneity factor (Cahalan 1994) L2 swath mock L3'
            dat_modis = (logW_mockL3 ./ meanW_mockL3);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            P=dat_modis;

            %add an extra row and column so that the points are
            %plotted propery (using the cell edges)
            P(end+1,:)=NaN;
            P(:,end+1)=NaN;

            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %
            title_info = '';

        case 'Optical Depth L2 swath'
            dat_modis = tau;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
         case 'Latitude L2 swath'
            dat_modis = Plat3_L2;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
case 'Longitude L2 swath'
            dat_modis = Plon3_L2;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P=dat_modis(row_L2_inds,col_L2_inds);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            


        case 'Cloud Top Height (km) L2 swath'
   
            
%             Wflag='calc';
%             dat_modis=N;
% 
%             %                        dat_modis= N - percent_error_Nd/100.*N;
% 
%             dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
% 
%             %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)
% 
%             P=dat_modis(row_L2_inds,col_L2_inds);

            

   % Doing this in filtering_data_L2 now.
%             [sst,dates,diff_days] = sst_for_date_range_McCoy_func(SST,sst_time,modis_date_time,'nearest');
%             sst_L2 = GRIDDATA(lon_sst,lat_sst,(squeeze(sst)),Plon_L2_5km,Plat_L2_5km);
% %            ctt_L2 = GRIDDATA(Plon_L2_5km,Plat_L2_5km,(squeeze(t_top)),Plon_L2,Plat_L2);            
%             CTH_L2 = (273.15 + sst_L2 - t_top - 2.35) / 0.0069 /1e3; %CTH from Zuidema (2009) in km
%             
            P = CTH_L2_5km;
%            P(ihtot)=NaN;


            MODIS_varname2_plot = [modis_data_plot date_str];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            inew_cticks=0;
            

            


        case 'Number of droplets (cm^{-3}) L2 swath'
            Wflag='calc';
            dat_modis=N;

            %                        dat_modis= N - percent_error_Nd/100.*N;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            

            i_dpcolor=1;

            if ifull_swath==0
                P=dat_modis(row_L2_inds,col_L2_inds);
            else
                P = dat_modis;
            end



            MODIS_varname2_plot = [modis_data_plot ' ' date_str ' for ' thresh_str];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';
            
            inew_cticks=1;

            
            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=1500;
            
            lb_map = lbmap(16,'brownblue'); colormap_choose = flipdim(lb_map,1);
            
            
    

        case 'Std dev Reff (\mum) annual mean'
            dat_modis=Cloud_Effective_Radius_Liquid_Standard_Deviation.timeseries3;


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Max Nd annual'
            dat_modis=Nd_timeseries.mean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)
            P=max(dat_modis,[],3);



            MODIS_varname2_plot = [modis_data_plot];

            units_str_plot='cm^{-3}';
            %                        units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Std dev Nd annual mean'
            dat_modis=Nd_timeseries.std_dev_norm;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            MODIS_varname2_plot = [modis_data_plot  '(normalised)'];

            units_str_plot='';      %no units if using normalised std dev
            title_info = '';

        case 'Std dev LWP annual mean'

            dat_modis=W_timeseries.std_dev_norm;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            MODIS_varname2_plot = [modis_data_plot  '(normalised)'];

            units_str_plot='g m^{-2}';
            units_str_plot='';      %no units if using normalised W
            title_info = '';

        case 'LWP annual mean'

            dat_modis=W_timeseries.mean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)



            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-2}';
            title_info = '';

        case 'GOES Nd'  
            [P]=MODIS_justN_func(goes_Tau,goes_Reff*1e-6,'calc',0,goes_Teff,'N');                          


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=1;

            
        case 'Map of 2D data from outside driver script'

%            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

%            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

            P = dat_modis;
            Npoints = prod(size(P));


            MODIS_varname2_plot = titlenam_driver;
            title_info='';

%            units_str_plot='cm^{-3}';
%            title_info = '';


%            inew_cticks=1;
            
  case 'Number of droplets cell values time mean'

            dat_modis=N_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average=            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints,P_std_dev] = meanNoNan(dat_modis2,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=0;      
            
 case 'Number of droplets cell values 3.7um time mean'

            dat_modis=N_time3_37;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average=            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints,P_std_dev] = meanNoNan(dat_modis2,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=0;               

  case 'Effective radius cell values 3.7um time mean'

            dat_modis=N_time3_37;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average=            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints,P_std_dev] = meanNoNan(dat_modis2,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=0;               
                      
            
        case 'Number of droplets 1km values time mean'
            %I.e., the nuber of droplets calculated using the 1km re and
            %tau (and then averaged to 1deg) rather than 1x1 deg values

            dat_modis=Droplet_Number_Concentration_37.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average - do before the time indices are used

            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints,P_std_dev] = meanNoNan(dat_modis2,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=0;

        case 'Number of droplets multi-year monthly means OLD'
            colormap_choose=jet; %default


            P=0;
            Npoints=0;
            Ndays2=0;
            imonths=[12 1 2];
            imonths=[6 7 8];
            imonths=7;
            %             imonths=month_no; %for multiple month plotting in monthly_means_from_plot_global

            iCTPs=[5 6 7];  %CTP values
            iCTP_screen=1; %CTP, CTP-std_dev or CTP-2*std_dev
            for iCTP2=1:length(iCTPs)
                iCTP=iCTPs(iCTP2);
                for imonth2=1:length(imonths)
                    imonth=imonths(imonth2);
                    for iyearav=1:length(Nd_multi)
                        P = P + Nd_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen).*Ndata_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                        Npoints = Npoints + Ndata_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                        Ndays2 = Ndays2 + Ndays_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                    end
                end
            end
            P=P./Npoints;
            ihtot=find(Npoints==0);

            P(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

            modis_CTP_screens
            MODIS_varname2_plot = [modis_data_plot ' MODIS AQUA+TERRA ' month_str ' ' CTP_str ' ' CTP_thresh_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=1;
            inew_cticks=1;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;

        case 'Number of days multi-year monthly means'
            colormap_choose=jet; %default


            P=0;
            Npoints=0;
            Ndays2=0;
            imonths=[12 1 2];
            imonths=[6 7 8];
            imonths=7;
            iCTPs=[5 6 7];  %CTP values
            iCTP_screen=1; %CTP, CTP-std_dev or CTP-2*std_dev
            for iCTP2=1:length(iCTPs)
                iCTP=iCTPs(iCTP2);
                for imonth2=1:length(imonths)
                    imonth=imonths(imonth2);
                    for iyearav=1:length(Nd_multi)
                        P = P + Nd_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen).*Ndata_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                        Npoints = Npoints + Ndata_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                        Ndays2 = Ndays2 + Ndays_multi(iyearav).data(:,:,imonth,iCTP,iCTP_screen);
                    end
                end
            end

            PDF_temp = squeeze(sum(sum(Nd_PDF_multi(:,:,imonths,iCTPs,iCTP_screen,:),3),4));

            P=Ndays2;
            P = sum(PDF_temp,3); %actually, this will be the number of
            %datapoints rather than ndays




            %            ihtot=find(Npoints==0);

            %            P(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

            modis_CTP_screens
            MODIS_varname2_plot = [modis_data_plot ' MODIS AQUA+TERRA ' month_str ' ' CTP_str ' ' CTP_thresh_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';


            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=5;

        case 'Number of droplets multi-year monthly medians'
            colormap_choose=jet; %default


            P=0;
            Npoints=0;
            
            if ~exist('ioverride_month_select_multi_years') | ioverride_month_select_multi_years==0
            imonths=[6 7 8];
            %             imonths=[12 1 2];
            imonths=[7]; %for multiple month plotting in monthly_means_from_plot_global
            %             imonths=month_no; %for multiple month plotting in monthly_means_from_plot_global

            imonths=[1:12];
            
            end
            
            iCTP_screen=1; %(=1) CTP, (=2) CTP-std_dev or (=3) CTP-2*std_dev
            iCTPs=[5 6 7];  %CTP values - [5 6 7] is all pressures since 5 means CTP>680,
            %6 means 680-440 and 7 means <=440
%CTP-std_dev >=880,830,780,730,680 then also CTP+std_dev <= 680 & CTP-std_dev > 440  (ISCPP
%classifications) and                  CTP+std_dev <= 440 & CTP-std_dev >
%50.

            %sum the PDF to combine for several months
            PDF_temp = squeeze(sum(sum(Nd_PDF_multi(:,:,imonths,iCTPs,iCTP_screen,:),3),4));
            Nd_bins_mon = [0:25:200 300:100:1400 1500:500:5000];  %29 bins
            prcs=[5 10 20 30 40 50 60 70 80 90 95];

            [prcs_out2,imed] = percentiles_from_PDF(Nd_bins_mon,PDF_temp,prcs,3);
            P=squeeze(prcs_out2(imed,:,:));
            Ndays2 = sum(PDF_temp,3); %actually, this will be the number of
            %datapoints rather than ndays
            dat_modis=P;
            %dat_modis=shiftdim(P,-1);
            %dat_modis=permute(dat_modis,[2 3 1]);

            %              %make sure that percentiles_from_PDF_RUN has been run first
            %              %(calculates prctiles from Nd_PDF_multi
            %             P=squeeze(prcs_out(imed,:,:,imonth,4,3));
            %             Ndays2 = sum(Nd_PDF_multi(:,:,imonth,4,3,:),6);

            %            P(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            modis_CTP_screens
            MODIS_varname2_plot = [modis_data_plot ' MODIS AQUA+TERRA ' month_str ' ' CTP_str ' ' CTP_thresh_str ' '];
            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=1;

            icontour=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;

        case 'Number of droplets multi-year monthly means'
            colormap_choose=jet; %default


            P=0;
            Npoints=0;
            
            if ~exist('ioverride_month_select_multi_years') | ioverride_month_select_multi_years==0
            imonths=[6 7 8];
            %             imonths=[12 1 2];
            imonths=[7]; %for multiple month plotting in monthly_means_from_plot_global
            %             imonths=month_no; %for multiple month plotting in monthly_means_from_plot_global

            imonths=[1:12];
            
            end
            
            iCTP_screen=1; %(=1) CTP, (=2) CTP-std_dev or (=3) CTP-2*std_dev
            iCTPs=[5 6 7];  %CTP values - [5 6 7] is all pressures since 5 means CTP>680,
            %6 means 680-440 and 7 means <=440
%CTP-std_dev >=880,830,780,730,680 then also CTP+std_dev <= 680 & CTP-std_dev > 440  (ISCPP
%classifications) and                  CTP+std_dev <= 440 & CTP-std_dev >
%50.
           %for ones with no pressure screening
           iCTP_screen=1;
           iCTPs=[1]; 

            %sum the PDF to combine for several months
            PDF_temp = squeeze(sum(sum(Nd_PDF_multi(:,:,imonths,iCTPs,iCTP_screen,:),3),4));
            Nd_bins_mon = [0:25:200 300:100:1400 1500:500:5000];  %29 bins
            Nd_bins_mon_mid = 0.5 * ( Nd_bins_mon(1:end-1) + Nd_bins_mon(2:end) );

            Nd_3D = repmat(Nd_bins_mon_mid,[size(PDF_temp,1) 1 size(PDF_temp,2)]);
            Nd_3D = permute(Nd_3D,[1 3 2]);
            
            P = sum(Nd_3D.*PDF_temp,3) ./ sum(PDF_temp,3);

%            [prcs_out2,imed] = percentiles_from_PDF(Nd_bins_mon,PDF_temp,prcs,3);
%            P=squeeze(prcs_out2(imed,:,:));
            Ndays2 = sum(PDF_temp,3); %actually, this will be the number of
            %datapoints rather than ndays
            dat_modis=P;
            %dat_modis=shiftdim(P,-1);
            %dat_modis=permute(dat_modis,[2 3 1]);

            %              %make sure that percentiles_from_PDF_RUN has been run first
            %              %(calculates prctiles from Nd_PDF_multi
            %             P=squeeze(prcs_out(imed,:,:,imonth,4,3));
            %             Ndays2 = sum(Nd_PDF_multi(:,:,imonth,4,3,:),6);

            %            P(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            modis_CTP_screens
            MODIS_varname2_plot = [modis_data_plot ' MODIS AQUA+TERRA ' month_str ' ' CTP_str ' ' CTP_thresh_str ' ' gcm_years_loaded_str ' '];
            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=1;

            icontour=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;
            
            
        case 'LWP multi-year monthly means'
            colormap_choose=jet; %default


            P=0;
            Npoints=0;
            
            if ~exist('ioverride_month_select_multi_years') | ioverride_month_select_multi_years==0
            imonths=[6 7 8]; %JJA
%            imonths=[12 1 2]; %DJF
            %imonths=[7]; %for multiple month plotting in monthly_means_from_plot_global
            %             imonths=month_no; %for multiple month plotting in monthly_means_from_plot_global

            imonths=[1:12];
            
            end
            
            iCTP_screen=1; %(=1) CTP, (=2) CTP-std_dev or (=3) CTP-2*std_dev
            iCTPs=[5 6 7];  %CTP values - [5 6 7] is all pressures since 5 means CTP>680,
            %6 means 680-440 and 7 means <=440
%CTP-std_dev >=880,830,780,730,680 then also CTP+std_dev <= 680 & CTP-std_dev > 440  (ISCPP
%classifications) and                  CTP+std_dev <= 440 & CTP-std_dev >
%50.
           %for ones with no pressure screening
           iCTP_screen=1;
           iCTPs=[1]; 

            %sum the PDF to combine for several months
            PDF_temp = squeeze(sum(sum(LWP_PDF_multi(:,:,imonths,iCTPs,iCTP_screen,:),3),4));
            Nd_bins_mon = [0:25:500 500:500:2000];  
            Nd_bins_mon_mid = 0.5 * ( Nd_bins_mon(1:end-1) + Nd_bins_mon(2:end) );

            Nd_3D = repmat(Nd_bins_mon_mid,[size(PDF_temp,1) 1 size(PDF_temp,2)]);
            Nd_3D = permute(Nd_3D,[1 3 2]);
            
            P = sum(Nd_3D.*PDF_temp,3) ./ sum(PDF_temp,3);

%            [prcs_out2,imed] = percentiles_from_PDF(Nd_bins_mon,PDF_temp,prcs,3);
%            P=squeeze(prcs_out2(imed,:,:));
            Ndays2 = sum(PDF_temp,3); %actually, this will be the number of
            %datapoints rather than ndays
            dat_modis=P;
            


            modis_CTP_screens
            MODIS_varname2_plot = [modis_data_plot ' MODIS AQUA+TERRA ' month_str ' ' CTP_str ' ' CTP_thresh_str ' ' gcm_years_loaded_str ' '];
            units_str_plot='g m^{-2}';
            title_info = '';


            inew_cticks=0;

            icontour=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=150;
            

            
            
        case 'Sea ice areal coverage'
            colormap_choose=jet; %default


            P = seaice_data*100; %convert to percentage

            date_str_seaice = datestr(seaice_day+datenum(seaice_year,1,1)-1);
           

            MODIS_varname2_plot = [modis_data_plot ' NIMBUS 7 ' date_str_seaice ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='%';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=100;

        case 'Time max SZA from mockL3, selected days'
            colormap_choose=jet; %default

            dat_modis=Solar_Zenith_Mean.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = max(dat_modis(:,:,time_inds_average),[],3); %time max (over requested time range)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;
            
        case 'Max minus Min Solar Zenith Angle cell values time mean - specific days'
             colormap_choose=jet; %default
             
             %this attempts to cover two plots - the max minus min from one
             %specific day and the max over the whole timeseries minus the
             %min
             
             if ~exist('Solar_Zenith_Maximum2')
                 dat_modis = Solar_Zenith_Mean.timeseries3;
                 dat_modis(dat_modis>81.4)=NaN;
             else

                 dat_modis=Solar_Zenith_Maximum2.timeseries3 - Solar_Zenith_Minimum.timeseries3;
             end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            if ~exist('Solar_Zenith_Maximum2')
                [P] = max(dat_modis(:,:,time_inds_average),[],3) - min(dat_modis(:,:,time_inds_average),[],3);
                Npoints=1;
            else
                [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            end
            

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
            title_nice = 'Max minus Min Solar Zenith Angle (degrees)';
            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='degrees';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
           
            
        case 'Max Solar Zenith Angle cell values time mean - specific days'
            colormap_choose=jet; %default
%            colormap_choose=colormap_choose(1:55,:); %make the end of
%            colormap less dark so that 
            %it isn't as close to the black used for NaN
            
            if ~exist('Solar_Zenith_Maximum2')
                dat_modis = Solar_Zenith_Mean.timeseries3;
                dat_modis(dat_modis>81.4)=NaN;
            else

                dat_modis=Solar_Zenith_Maximum2.timeseries3;
            end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average


            
            
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            if ~exist('Solar_Zenith_Maximum2')
                P = max(dat_modis(:,:,time_inds_average),[],3);
                Npoints = 1;
            else
                [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            end

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
            %                        MODIS_varname2_plot = modis_data_plot;
            title_nice = 'Max Solar Zenith Angle (degrees)';

            units_str_plot='degrees';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=81.4;
            
        case 'Min Solar Zenith Angle cell values time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=Solar_Zenith_Minimum.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='degrees';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
        case {'MODIS 1.6\mum minus POLDER effective radius - using colocated POLDER','MODIS 2.1\mum minus POLDER effective radius - using colocated POLDER','MODIS 3.7\mum minus POLDER effective radius - using colocated POLDER'}
 %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6
  %Re_3.7' in pdf2d_plot_commands)
  
          switch modis_data_plot
              case 'MODIS 1.6\mum minus POLDER effective radius - using colocated POLDER'
                  %just directly do the difference of the means (will give the same
                  %result as the mean of all the differences - mean(xi-yi) =
                  %mean(xi) - mean(yi)
                  %But can't deal with different time selections with this
                  %approach
                  P = meanNoNan(Cloud_Effective_Radius_16_Liquid_Mean.timeseries3,3) - meanNoNan(Par2_CDR_coloc,3);                                   
              case 'MODIS 2.1\mum minus POLDER effective radius - using colocated POLDER'
                  P = meanNoNan(Cloud_Effective_Radius_Liquid_Mean.timeseries3,3) - meanNoNan(Par2_CDR_coloc,3);                                   
              case 'MODIS 3.7\mum minus POLDER effective radius - using colocated POLDER'
                  P = meanNoNan(Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,3) - meanNoNan(Par2_CDR_coloc,3);                                   
          end
          Npoints = zeros(size(P)); %don't need this
  
%            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

%            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=5;  


            
        case {'MODIS 1.6\mum minus POLDER effective radius - using mean POLDER over all years','MODIS 2.1\mum minus POLDER effective radius - using mean POLDER over all years','MODIS 3.7\mum minus POLDER effective radius - using mean POLDER over all years'}
              %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6
  %Re_3.7' in pdf2d_plot_commands)
  
          switch modis_data_plot
              case 'MODIS 1.6\mum minus POLDER effective radius - using mean POLDER over all years'
                  %just directly do the difference of the means (will give the same
                  %result as the mean of all the differences - mean(xi-yi) =
                  %mean(xi) - mean(yi)
                  %But can't deal with different time selections with this
                  %approach
                  P = meanNoNan(Cloud_Effective_Radius_16_Liquid_Mean.timeseries3,3) - meanNoNan(permute(daymean_Par2_CDR(:,ilat_loaded,ilon_loaded),[2 3 1]),3);                  
              case 'MODIS 2.1\mum minus POLDER effective radius - using mean POLDER over all years'                 
                  P = meanNoNan(Cloud_Effective_Radius_Liquid_Mean.timeseries3,3) - meanNoNan(permute(daymean_Par2_CDR(:,ilat_loaded,ilon_loaded),[2 3 1]),3);                 
              case 'MODIS 3.7\mum minus POLDER effective radius - using mean POLDER over all years'                 
                  P = meanNoNan(Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,3) - meanNoNan(permute(daymean_Par2_CDR(:,ilat_loaded,ilon_loaded),[2 3 1]),3);                                   
          end
          Npoints = zeros(size(P)); %don't need this
  
%            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

%            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=5;  


       
            
            
            
        case {'POLDER effective radius - using colocated data','POLDER effective radius - using individual days','POLDER effective radius - using mean over all years','POLDER effective radius - using individual days from VOCALS save file'}
            colormap_choose=jet; %default

            %re-arrange to put time at the end
            switch modis_data_plot
                case 'POLDER effective radius - using mean over all years'
                    dat_modis=permute(daymean_Par2_CDR,[2 3 1]);  %this is all days over all the years - usually runs out of memory
                    MODIS_varname2_plot = [modis_data_plot ' POLDER ' time_mean_str years_parasol_str ' '];
                case 'POLDER effective radius - using individual days'
                    dat_modis=permute(daymeanALL_Par2_CDR,[2 3 1]); %this is the average over all years
                    MODIS_varname2_plot = [modis_data_plot ' POLDER ' time_mean_str years_parasol_str ' '];
                case 'POLDER effective radius - using colocated data'
                    dat_modis=Par2_CDR_coloc;
                    MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
                case 'POLDER effective radius - using individual days from VOCALS save file'
                    dat_modis=daymean_Par2_CDR;
                    MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];  
            end

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)




            %                        MODIS_varname2_plot = modis_data_plot;

            titlenam =  modis_data_plot;
            units_str_plot='\mum';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=6;
            iset_max_clim=1;
            clim_max=25;


        case 'Cloud_Fraction_Day_Mean'
            %MOD35 clodu fraction
            imake_out_of_range_CTH_zero = 0; %useful for high cloud fraction
            thresh_zeroCF = 0.05;
            
            colormap_choose=jet; %default

            dat_modis=Cloud_Fraction_Day_Mean.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
            
            if imake_out_of_range_CTH_zero==1
                thresh_CTH_high = 3.2;
                inan = find(isnan(CTH.timeseries3)==1);
                dat_modis(inan)=NaN;
                icth = find(CTH.timeseries3<thresh_CTH_high);
                dat_modis(icth)=0;
                izero = find(Cloud_Fraction_Day_Mean.timeseries3 < thresh_zeroCF);
                dat_modis(izero) = Cloud_Fraction_Day_Mean.timeseries3(izero);
            end

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' (MOD35) MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            if imake_out_of_range_CTH_zero==1
                MODIS_varname2_plot = [MODIS_varname2_plot ' with CTH.LT.' num2str(thresh_CTH_high) ' made zero with no screening for CF.LT.' num2str(thresh_zeroCF)];
            end
            
            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;

            
        case 'Number of droplets differences between datasets'
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=0;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;


            end
            
            inew_cticks=0;
            iset_min_clim=1;
            clim_min=-50;
            iset_max_clim=1;
            clim_max=50;

            clim_min=-1;
            clim_max=1;

            colormap_choose=jet; %default

            %first want to match the times for the dataset that the
            %screening will NOT be based upon (e.g. G14.blah)
            
            %Nd
           [dat_out] = time_match_data({daynum_timeseries3_MODIS,modisyear_timeseries3},G14.Date_Time_Swath.timeseries3,G14.N_time3);                                    
            dat_out=flipdim(dat_out,1); %data in the files is backwards in the lat dimension
               %compared to L3 data!!
           dat_modis=N_time3 - dat_out;
            
%             %CFliq
%             [dat_out] = time_match_data({daynum_timeseries3_MODIS,modisyear_timeseries3},G14.Date_Time_Swath.timeseries3,G14.Cloud_Fraction_Liquid.timeseries3);
%             dat_out=flipdim(dat_out,1); %data in the files is backwards in the lat dimension
%                %compared to L3 data!!
%             dat_modis=Cloud_Fraction_Liquid.timeseries3 - dat_out;

            %mean SZA
%             [dat_out] = time_match_data({daynum_timeseries3_MODIS,modisyear_timeseries3},G14.Date_Time_Swath.timeseries3,G14.Solar_Zenith_Mean.timeseries3);
%             dat_out=flipdim(dat_out,1); %data in the files is backwards in the lat dimension
%                %compared to L3 data!!
%             dat_modis=Solar_Zenith_Mean.timeseries3 - dat_out;
            
            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            Ndays2 = Npoints;


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            
            if icontour==1
                switch cont_dat_choose
                    case{'calipso highCF','calipso mid+highCF','calipso lowCF'}
                        % -- Use CALIPSO mid+high CF to define a certain region
                        % But need to grid the CALIPSO data onto the MODIS grid
                        clh = meanNoNan(clhcalipso_monthly,1)/100;
                        clm = meanNoNan(clmcalipso_monthly,1)/100;
                        cll = meanNoNan(cllcalipso_monthly,1)/100;
                end
                

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');

                end
                
            end


            
        case 'SMOS soil moisture - specific days'
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=0;


                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=1;
                
                plot_num_datapoints=0;


            end

            colormap_choose=jet; %default

            dat_modis=Soil_Moisture.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
%            dat_modis2 = dat_modis(:,:,time_inds_average);

%itime_SMOS = [30];
time_inds_average = itime_SMOS;
            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

            Ndays2 = Npoints;

day_str = num2str(time_inds_average,'%03g');
time_mean_str=date_from_day_of_year_func(days_to_read(time_inds_average),years_to_read(time_inds_average));

            MODIS_varname2_plot = [modis_data_plot ' SMOS ' time_mean_str ' ' modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;
            
            if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' SMOS day ' time_mean_str ', ' modisyears_str ' '];
             end            

            units_str_plot='m^3 m^{-3}';
            title_info = '';

            
          
            
    
            

        case 'Number of droplets cell values time mean - specific days'
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=1;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;
                
                plot_num_datapoints=0;


            end

            colormap_choose=jet; %default

            dat_modis=N_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            Ndays2 = Npoints;


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;
            
            if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end            

            units_str_plot='cm^{-3}';
            title_info = '';

            
            if icontour==1
                switch cont_dat_choose
                    case{'calipso highCF','calipso mid+highCF','calipso lowCF'}
                        % -- Use CALIPSO mid+high CF to define a certain region
                        % But need to grid the CALIPSO data onto the MODIS grid
                        clh = meanNoNan(clhcalipso_monthly,1)/100;
                        clm = meanNoNan(clmcalipso_monthly,1)/100;
                        cll = meanNoNan(cllcalipso_monthly,1)/100;
                end
                

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');

                end
                
            end
            
  case 'Number of droplets QA cell values time mean - specific days'
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=1;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;
                
                plot_num_datapoints=0;


            end

            colormap_choose=jet; %default

            dat_modis=N_QA_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            [P,Npoints,Nd_std_dev] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            Ndays2 = Npoints;


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;
            
            if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end            

            units_str_plot='cm^{-3}';
            title_info = '';

            
            if icontour==1
                switch cont_dat_choose
                    case{'calipso highCF','calipso mid+highCF','calipso lowCF'}
                        % -- Use CALIPSO mid+high CF to define a certain region
                        % But need to grid the CALIPSO data onto the MODIS grid
                        clh = meanNoNan(clhcalipso_monthly,1)/100;
                        clm = meanNoNan(clmcalipso_monthly,1)/100;
                        cll = meanNoNan(cllcalipso_monthly,1)/100;
                end
                

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');

                end
                
            end
            
            
            
        case 'Number of droplets QA pixel weighted cell values time mean - specific days'
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=0;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;
                
                plot_num_datapoints=0;


            end

            colormap_choose=jet; %default

            %Apparently monthly averages (M3) prdoucts are weighted
            %according the numner of pixels in a 1x1 box rather than giving
            %equal weight to each day.
            Npix = Cloud_Fraction_Liquid_Pixel_Counts.timeseries3;
            dat_modis=N_QA_time3.*Npix;


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            Npix(ihtot)=NaN;

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
%            dat_modis2 = dat_modis(:,:,time_inds_average);
%           Npix2 = Npix(:,:,time_inds_average);

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3,'sum'); %time mean
            Npix_tot = meanNoNan(Npix(:,:,time_inds_average),3,'sum'); %sum over time (of non-NaN values)
            Npoints = Npix_tot;
            P = P ./ Npix_tot; %Divide by total no. pixels to get the mean
            
            Ndays2 = Npoints;


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;
            
            if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end            

            units_str_plot='cm^{-3}';
            title_info = '';

            
            if icontour==1
                switch cont_dat_choose
                    case{'calipso highCF','calipso mid+highCF','calipso lowCF'}
                        % -- Use CALIPSO mid+high CF to define a certain region
                        % But need to grid the CALIPSO data onto the MODIS grid
                        clh = meanNoNan(clhcalipso_monthly,1)/100;
                        clm = meanNoNan(clmcalipso_monthly,1)/100;
                        cll = meanNoNan(cllcalipso_monthly,1)/100;
                end
                

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');

                end
                
            end

            

        case 'Latest or earliest day after screening'
            %After ihtot and time screening is done what is the latest day
            %left for each location?
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0

                icontour=0;
                cont_dat_choose='';

                inew_cticks=0;


                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=150;
                
                plot_num_datapoints=0;
                
                latest_or_earliest='latest';


            end

            colormap_choose=jet; %default

            dat_modis=repmat(daynum_timeseries3,[size(N_time3,1) 1 size(N_time3,2)]);
            dat_modis=permute(dat_modis,[1 3 2]);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            
            dat_modis2 = dat_modis(:,:,time_inds_average);

            switch latest_or_earliest
                case 'latest'
                    P =max(dat_modis(:,:,time_inds_average),[],3); %latest day available
                case 'earliest'
                    P =min(dat_modis(:,:,time_inds_average),[],3); %latest day available                    
            end

            %Ndays2 = Npoints;


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' ' thresh_str];
            save_notes_filepath = write_notes_to_file([savenotes_filedir 'seasonal_cycle_fig_notes'],[MODIS_varname2_plot] );
%       write_notes_to_file( save_notes_filepath, [ 'WITH CHANGES :- ' eval(['screen_edits_' tag '{1}']) ] );

            %                        MODIS_varname2_plot = modis_data_plot;
            
            if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end            

            units_str_plot='cm^{-3}';
            title_info = '';

            
            if icontour==1
                switch cont_dat_choose
                    case{'calipso highCF','calipso mid+highCF','calipso lowCF'}
                        % -- Use CALIPSO mid+high CF to define a certain region
                        % But need to grid the CALIPSO data onto the MODIS grid
                        clh = meanNoNan(clhcalipso_monthly,1)/100;
                        clm = meanNoNan(clmcalipso_monthly,1)/100;
                        cll = meanNoNan(cllcalipso_monthly,1)/100;
                end
                

                switch cont_dat_choose
                    case 'calipso highCF'
                        cont_dat = clh;
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
                    case 'calipso mid+highCF'
                        %                    filt_dat = clm + clh;
                        %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                        cont_dat = clm + clh - (clm.*clh);
                        cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');

                end
                
            end

            

            
        case 'Number of droplets cell values time MINIMUM - specific days'
            colormap_choose=jet; %default

            dat_modis=N_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = min(dat_modis(:,:,time_inds_average),[],3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=1;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=100;
            
            
        case 'Number of droplets 1.6um cell values time mean - specific days'
            colormap_choose=jet; %default
            
            

            dat_modis=N_time3_16;


                            
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
                       
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=0;

            inew_cticks=1;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
 case 'Number of droplets 3.7um cell values time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=N_time3_37;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=0;

            inew_cticks=1;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
 case 'Number of droplets 3.7um minus 2.1um cell values time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=N_time3_37 - N_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=150;
            
            
        case 'Nd ALL or low SZA - specific days'
            colormap_choose=jet; %default
            
            band_str=''; band_str2='2.1 \mum';
%            band_str='_16'; band_str2='1.6 \mum';
            band_str='_37'; band_str2='3.7 \mum';
            
            sza_str ='_allSZA';
%            sza_str ='_lowSZA';            

            dat_modis=eval(['squeeze(N_time3' band_str ');']);
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            %remove the zero % differences for clarity?            
%            irem = find(abs(P)<3);% | isnan(P)==1 );
%            P=Npoints;
%            P(irem)=NaN;

%sza_diff_vs_LAT = abs(Solar_Zenith_Maximum_allSZA.timeseries3(:,:,time_inds_average) - Solar_Zenith_Maximum_lowSZA.timeseries3(:,:,time_inds_average));
%max_sza_diff_vs_LAT = max(sza_diff_vs_LAT,[],2);
%nvals_sza_diff_ALL = sza_diff_vs_LAT>0;
%nvals_sza_diff = mean(sum(nvals_sza_diff,2),3);

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' band_str2];
            title_nice = ['Nd for ' time_mean_str ' (' band_str2 '_' sza_str ')'];
            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

            icontour=0;

            inew_cticks=1;


            iset_min_clim=0;
            clim_min=-15;
            iset_max_clim=0;
            clim_max=15;            
            
            
            
        case 'Percentage Nd difference ALL vs low SZA - specific days'
            colormap_choose=jet; %default
            
            band_str='21'; band_str2='2.1 \mum';
%            band_str='16'; band_str2='1.6 \mum';
            band_str='37'; band_str2='3.7 \mum';

            dat_modis=eval(['squeeze(100*(Nd_' band_str '_allSZA-Nd_' band_str '_lowSZA)./Nd_' band_str '_allSZA);']);

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            %remove the zero % differences for clarity?            
%            irem = find(abs(P)<3);% | isnan(P)==1 );
%            P=Npoints;
%            P(irem)=NaN;

sza_diff_vs_LAT = abs(Solar_Zenith_Maximum_allSZA.timeseries3(:,:,time_inds_average) - Solar_Zenith_Maximum_lowSZA.timeseries3(:,:,time_inds_average));
max_sza_diff_vs_LAT = max(sza_diff_vs_LAT,[],2);
%nvals_sza_diff_ALL = sza_diff_vs_LAT>0;
%nvals_sza_diff = mean(sum(nvals_sza_diff,2),3);

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' band_str2];
            title_nice = ['Percentage Nd difference for ' time_mean_str ' (' band_str2 ')'];
            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='%';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=-15;
            iset_max_clim=1;
            clim_max=15;            
       
   case 'Nd ALL SZA, but no low SZA - specific days'
            colormap_choose=jet; %default
            
            band_str='21'; band_str2='2.1 \mum';
            band_str='16'; band_str2='1.6 \mum';
            band_str='37'; band_str2='3.7 \mum';
            
            dat_modisLOW=eval(['squeeze(Nd_' band_str '_lowSZA);']);
            dat_modisLOW(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            
            dat_modis=eval(['squeeze(Nd_' band_str '_allSZA);']);
%            dat_modis = dat_modisLOW;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
%ihtot doesn't remove anything here            
            
            irem = find(isnan(dat_modisLOW)==0);
%            dat_modis(irem)=NaN;
            
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            %remove the zero % differences for clarity?            
         
%            P=Npoints;
%            P(irem)=NaN;


            

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' band_str2];
            title_nice = ['Percentage Nd difference for ' time_mean_str '(' band_str2 ')'];
            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='%';
            title_info = '';

            icontour=0;

            inew_cticks=1;


            iset_min_clim=0;
            clim_min=-15;
            iset_max_clim=0;
            clim_max=15;            
       
            
            

        case 'Cloud depth cell values time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=H_time3*sqrt(1/1.15/0.7); modis_data_plot=[modis_data_plot ' Painemal '];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='m';
            title_info = '';

            icontour=1;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=150;
            iset_max_clim=1;
            clim_max=350;

        case 'MODIS Cloud Fraction minus CALIPSO CF'
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                MODCF_type = 'MOD06 liquid';
%                MODCF_type = 'MOD35';
                
                %what range of clouds to compare MODIS to (for bias calc)
                cal_CF_range='low';
                cal_CF_range='low+mid+high';     %Should reconsider using this as there may be cloud overlap
                  % Now have the cltcalipso field for daytime too, which I
                  % think should give the total column CF (i.e. including
                  % overlaps).
                
              %whether to filter our data where there is a the CALIPSO
              %high cloud CF is large (or mid+high - choose using cont_dat_choose below)
                ifilter_clhcalipso=1;
                thresh_clh = [-0.01 0.3];                               
                thresh_clh = [0.3 1.01]; 

                zero_CF_out_of_range = 0; %flag that causes CFs to be set to zero if the CTH is determined to be outside
                % of the CTH range (rather than being ingnored as NaN.
                % Important when looking at high clouds since if CTH is low then there is no high
                % cloud and so should count as CF=0 rather than ignoring it.
                % For low cloud is more dubious perhaps - except that this
                % is what CALIPSO would do since higher cloud would block
                % low cloud. Gives better match if use this for MOD35 vs CALIPSO comparison. 
                
                iocean_only=1; %only plot/calculate for the ocean only points.

                %What range to filter out, and also what range to plot if plotting a
                %contour
                 cont_dat_choose = 'calipso highCF';
                 cont_dat_choose = 'calipso mid+highCF';
            end
            
            %the landmask from CAMCLUBBv2.
            load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
            
             MODIS_varname2_plot = [MODCF_type ' minus CALISPO CF ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str];

             colormap_choose=jet; %default

             %            dat_modis=H_time3;

             MODIS_varname2_plot = [MODIS_varname2_plot ' for ' cal_CF_range ' CALIPSO CF'];




             %CALIPSO data is 90x180 - i.e. 2x2 deg. So interpolate the
             %MODIS data to the same grid
             % ----------------------------
             interp_to_calipso_grid  %N.B. - this just creates the grid to interpolate to. This is done later
             % ----------------------------            %using e.g. modis_cf_grid = griddata(MLON,MLAT,P2,X,Y);

             Plat2D = Y;
             Plon2D = X;

             Plon2D_edges = X2;
             Plat2D_edges = Y2;
             
             % Find the land points for removal if requried
             iland=find(gcm_landmask_cal>0.01 | isnan(gcm_landmask_cal)==1 );
             
                         %remove data where there is significant high level CALIPSO
            %cloud
            clh = meanNoNan(clhcalipso_monthly,1)/100;
            clm = meanNoNan(clmcalipso_monthly,1)/100;
            switch cont_dat_choose
                case 'calipso highCF'
                    filt_dat = clh;
                case 'calipso mid+highCF'
%                    filt_dat = clm + clh;
  %Prob better to use the randome overlap assumption  CF = CF1+CF2-CF1*CF2
                    filt_dat = clm + clh - (clm.*clh);                    
            end
                    
                    switch modis_calipso_cf_case
                        case {'mean time bias','calc and save MODIS monthly data'}

                    switch MODCF_type
                        case 'MOD06 liquid'
                            dat_modis=Cloud_Fraction_Liquid.timeseries3; modis_data_plot=[modis_data_plot];
                        case 'MOD35'
                            dat_modis=Cloud_Fraction_Day_Mean.timeseries3; modis_data_plot=[modis_data_plot];
                    end

                    dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

                    %if outside of the CTH range make zero - assume that there is none of the
                    %desired cloud
                    if zero_CF_out_of_range==1
                        dat_modis(CTH.timeseries3<thresh_CTH(1) | CTH.timeseries3>=thresh_CTH(2)) = 0;
                        MODIS_varname2_plot = [MODIS_varname2_plot ' CTHs out of range zeroed'];
                    end

                    dat_modis = dat_modis(:,:,time_inds_average);
                    
                    coarsen_method = 'interpolate'; %Current method
                    coarsen_method = 'average'; %Trial method                    
                    switch coarsen_method
                        case 'interpolate'

                    [P2,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


                    %in case are worried about the averaging to 2x2 degrees
                    %affecting things
                    %             P3 = reduce_matrix_subsample_mean(P2,2,2); %smooth out to 2x2 degree averages
                    %             MLON3 = MLON([1:2:end]);
                    %             MLAT3 = MLAT([1:2:end]);
                    %             modis_cf_grid = griddata(MLON3,MLAT3,P3,X,Y);

                    %but makes little difference, so stick with original formulation without
                    %2x2 deg averaging
                    modis_cf_grid = griddata(MLON,MLAT,P2,X,Y);



                        case 'average'
                            siz=size(dat_modis);
                            P3 = NaN*ones([siz(1)/2 siz(2)/2 siz(3)]);
                            for it_smooth=1:size(dat_modis,3)
                                P3(:,:,it_smooth) = reduce_matrix_subsample_mean(dat_modis(:,:,it_smooth),2,2);
                            end
                            
                            [P2,Npoints] = meanNoNan(P3,3); %time mean (all times)
                            
                            MLAT2 = MLAT(1)-0.5:-2:MLAT(end)-0.5;
                            MLON2 = MLON(1)+0.5:2:MLON(end)+0.5;
                            
                            %Just doing the interpolation here to make it
                            %match the global size of CALIPSO
                            modis_cf_grid = griddata(MLON2,MLAT2,P2,X,Y);
                    end
                    





                    cal_str = '_DAYTIME2'; %for the MapLowMidHigh data
                    cal_str = ''; %for the cll* etc files
                    switch cal_CF_range
                        case 'low'
                            calipsoCF = eval(['meanNoNan(cllcalipso_monthly' cal_str '(time_index_calipso,:,:),1)/100;']);
                        case 'low+mid+high'
                            calipsoCF = eval(['meanNoNan(cllcalipso_monthly' cal_str '(time_index_calipso,:,:) + clmcalipso_monthly' cal_str '(time_index_calipso,:,:) + clhcalipso_monthly' cal_str '(time_index_calipso,:,:),1)/100;']);
                    end


                    switch modis_calipso_cf_case
                        case 'calc and save MODIS monthly data'
                            modis_cf_grid_save(time_index_calipso,:,:,iCTH_set,ical_CF,iMOD_set,iscreen_set) = modis_cf_grid;
                    end


                    P = modis_cf_grid - calipsoCF;

                case {'RMSE from monthly data','Correlation coeff from monthly data'}
                    %CALIPSO data for required height range
                    cal_str = '_DAYTIME2'; %for the MapLowMidHigh data
                    cal_str = ''; %for the cll* etc files
                    switch cal_CF_range
                       % N.B. time_index_calipso is set to 1:48 in this
                       % case (i.e. all indices)
                        case 'low'
                            calipsoCF = eval(['cllcalipso_monthly' cal_str '(time_index_calipso,:,:)/100;']);
                        case 'low+mid+high'
                            calipsoCF = eval(['(cllcalipso_monthly' cal_str '(time_index_calipso,:,:) + clmcalipso_monthly' cal_str '(time_index_calipso,:,:) + clhcalipso_monthly' cal_str '(time_index_calipso,:,:) )/100;']);
                    end

%                    [P2,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)

                    modis_cf_dat = squeeze( modis_cf_grid_save(:,:,:,iCTH_set,ical_CF,iMOD_set,iscreen_set) );
                    
                     if iocean_only==1
                         gcm_landmask_cal_ND = repmat(gcm_landmask_cal,[1 1 48]);
                         gcm_landmask_cal_ND = permute(gcm_landmask_cal_ND,[3 1 2]);
                         iland_ND = find(gcm_landmask_cal>0.01 | isnan(gcm_landmask_cal)==1 );
                         calipsoCF(iland)=NaN;
                     end
                     
                     %Screening for CALIPSO cloud regions
                     filt_dat_ND = repmat(filt_dat,[1 1 48]);
                     filt_dat_ND = permute(filt_dat_ND,[3 1 2]);
                     iclh = find(filt_dat_ND<=thresh_clh(1) | filt_dat_ND>thresh_clh(2));
                     calipsoCF(iclh)=NaN;
                     
                     dat_sq = (modis_cf_dat - calipsoCF ).^2;
                     
                     switch modis_calipso_cf_case                        
                         case 'Correlation coeff from monthly data'
                             %Correlation between MODIS and CALIPSO monthly CF
                             %Find all NaNs for monthly data (either MODIS or
                             %CALIPSO, so working on the sm of the two)
                             %Thi is the overall corr - could also do a map
                             
                             P = NaN*ones([size(calipsoCF,2) size(calipsoCF,3)]);
                             
                             for ilat_corr=1:size(calipsoCF,2)
                                 for ilon_corr=1:size(calipsoCF,3)
                                     inan=find(isnan( modis_cf_dat(:,ilat_corr,ilon_corr) + calipsoCF(:,ilat_corr,ilon_corr) )==1);
                                     calipsoCF_noNan=calipsoCF(:,ilat_corr,ilon_corr); calipsoCF_noNan=calipsoCF_noNan(:);
                                     calipsoCF_noNan(inan)=[]; %remove the NaN elements
                                     modis_cf_dat_noNan = modis_cf_dat(:,ilat_corr,ilon_corr); modis_cf_dat_noNan=modis_cf_dat_noNan(:);
                                     modis_cf_dat_noNan(inan)=[];

                                     if length(calipsoCF_noNan)>10  %Only do the correlations if have some values left
                                            %And enough to be reasonably
                                            %significant.
                                         P(ilat_corr,ilon_corr) = corr(modis_cf_dat_noNan,calipsoCF_noNan);
                                     end
                                     

                                 end
                                 
                             end
                             
                           
                             

                         case 'RMSE from monthly data'
                             P = sqrt (meanNoNan(dat_sq , 1) ); %RMS difference of monthly data                             
                     end
                             
                     Psum = meanNoNan( dat_sq(:) , 1);
                     RMSE_overall = sqrt(Psum); %Overall value for the whole domain

            end
            
            
 
            
            if iocean_only==1
                P(iland)=NaN;                
            end            
            

            
             if ifilter_clhcalipso==1                          
                iclh = find(filt_dat<=thresh_clh(1) | filt_dat>thresh_clh(2));
                P(iclh) = NaN;
             end
             
               %The mean corr coeff is skewed by negative
                             %values. Know that r=1 is the correct perfect
                             %score, so can see how far from this we are.
                             %So doing 1-P - negative values will give a >1
                             %answer and r-=1 will give s score of zero
                             %(the lower the score the better).
                             score_from_one = meanNoNan(1-P(:),1);
             
            
            
             %take out the NaNs so that can calc correlation coefficients 
             inan=find(isnan(P)==1);
             calipsoCF_noNan=calipsoCF(:);
             calipsoCF_noNan(inan)=[]; %remove the NaN elements
             bias_noNan = P(:); 
             bias_noNan(inan)=[];
             modis_cf_dat_noNan = modis_cf_dat(:);
             modis_cf_dat_noNan(inan)=[];
             
             switch modis_calipso_cf_case  
                 
                 case {'mean time bias','calc and save MODIS monthly data'}

                     % correlation coefficient between the bias and CF from CALIPSO

                     [RHO,PVAL]=corr(bias_noNan,calipsoCF_noNan);
                     
                      mean_bias_MODIS = meanNoNan(P(:),1);
                      RMSE_MODIS = sqrt(meanNoNan(P(:).^2,1));


             

                 case 'RMSE from monthly data'
                     %Correlation between MODIS and CALIPSO monthly CF
                     %Find all NaNs for monthly data (either MODIS or
                     %CALIPSO, so working on the sm of the two)
                     %Thi is the overall corr - could also do a map
                      inan=find(isnan(modis_cf_dat + calipsoCF)==1);
                      calipsoCF_noNan=calipsoCF(:);
                      calipsoCF_noNan(inan)=[]; %remove the NaN elements                     
                      modis_cf_dat_noNan = modis_cf_dat(:);
                      modis_cf_dat_noNan(inan)=[];
                    
                     [RHO,PVAL]=corr(modis_cf_dat_noNan,calipsoCF_noNan);   
                     
                      mean_bias_MODIS = meanNoNan(P(:),1);
                      RMSE_MODIS = sqrt(Psum); %Overall RMSE from all monthly values and for all locations.
                     

             end
             
             
             
            icontour=1;
%            contour_label=0; %switch off contour labelling for this plot
%actually may as well leave the labelling on to make it clear to me what CF
%threshold is being used.
if icontour==1
       
            switch cont_dat_choose
                case 'calipso highCF'
%                    cont_dat = clh;
                    cont_dat = clh;                    
                    cont_ints=[0.3 0.3];  
                    %                    cont_ints=[0.2 0.2];
                    
                    MODIS_varname2_plot = [MODIS_varname2_plot ' with ' cont_dat_choose ' contours'];
                case 'calipso mid+highCF'
                    %                    cont_dat = clh;
                    cont_dat = clm+clh;
                    cont_ints=[0.3 0.3];
                    %                    cont_ints=[0.2 0.2];
                    MODIS_varname2_plot = [MODIS_varname2_plot ' with ' cont_dat_choose ' contours'];                    
            end
            
end

             

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

           

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=-0.3;
            iset_max_clim=0;
            clim_max=0.3;
            
        case 'Model minus MODIS Nd'            
            %May want to keep the option of doing the CALIPSO mid+high
            %cloud screening
            
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                MODCF_type = 'MOD06 liquid';
                MODCF_type = 'MOD35';

                error_type = 'absolute';
                error_type = 'percentage';
                
                %what range of clouds to compare MODIS to (for bias calc)
                cal_CF_range='low';
                cal_CF_range='low+mid+high';     %Should reconsider using this as there may be cloud overlap
                  % Now have the cltcalipso field for daytime too, which I
                  % think should give the total column CF (i.e. including
                  % overlaps).
                
              %whether to filter our data where there is a the CALIPSO
              %high cloud CF is large (or mid+high - choose using cont_dat_choose below)
                ifilter_clhcalipso=1;
                thresh_clh = [-0.01 0.3];                               
%                thresh_clh = [0.3 1.01]; 

                zero_CF_out_of_range = 0; %flag that causes CFs to be set to zero if the CTH is determined to be outside
                % of the CTH range (rather than being ingnored as NaN.
                % Important when looking at high clouds since if CTH is low
                % then there is no high
                % cloud and so should count as CF=0 rather than ignoring it.
                % For low cloud is more dubious perhaps - except that this
                % is what CALIPSO would do since higher cloud would block
                % low cloud. Gives better match if use this for MOD35 vs CALIPSO comparison. 
                
                iocean_only=0; %only plot/calculate for the ocean only points.

% what to plot if plotting contours
                 cont_dat_choose = 'calipso highCF';
                 cont_dat_choose = 'calipso mid+highCF';
                 cont_dat_choose = 'calipso lowCF';     
                 cont_dat_choose = 'MOD06 CF';                 
%                 cont_dat_choose = 'MOD35 CF';  
                 cont_dat_choose = 'MOD06 CF minus 1 std dev';
                 cont_dat_choose = 'MOD35 CF minus 1 std dev';                 
                 
%What to use for filtering                 
                 filt_dat_choose = 'calipso mid+highCF';
                 
                 iuse_saved_mask=0;
                 saved_mask_name = 'imask_ndays15_cal_midhigh';
                 %save the .mat files so that this works:-
%                 saved_mask_filepath = ['~/' saved_mask_name '.mat'];

            end
            
            %the landmask from CAMCLUBBv2.
            load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
            switch error_type
                case 'absolute'
                    varname_str = ' minus MODIS N_d (cm^{-3}) ';
                case 'percentage'
                    varname_str = ' vs MODIS N_d percentage difference (%) ';                    
            end
            
            MODIS_varname2_plot = [MODCF_type varname_str aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str];

            colormap_choose=jet; %default

            
% -- May need some interpolation to get MODIS and model on the same grid.  
         
            Nd_model = eval(['gcm_Nd_max_screen_' gcm_str]);
%            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
             Nd_model = Nd_model + time_inds_average2; %do the time screening for the model
            Nd_model = meanNoNan(Nd_model,1); %average over all times
            %Might want to do some model screening? E.g. for CF? Or show
            %the effect of this elsewhere and ignore here?
            %Also want to screen the MODIS data, so will need to deal with
            %two lots of screenings here.
            
            %assuming that modis screening is using ihtot here... probably
            %wrong
            modis_bias_comp = N_time3(:,:,:,1); %daytime only overpasses           
            modis_bias_comp(ihtot)=NaN; %make these values NaN as they will then be removed from the average                                                              
            modis_bias_comp = modis_bias_comp(:,:,time_inds_average);
                                    
            [Pmodis,Npoints] = meanNoNan(modis_bias_comp,3); %this tells us how much MODIS data we have in case
            %we want to screen for this
            
            P = Nd_model - Pmodis;
            
            
            switch error_type
                case 'percentage'
                    P = 100* P ./ meanNoNan(Pmodis,3);
            end
            

            
% -- Use CALIPSO mid+high CF to define a certain region
% But need to grid the CALIPSO data onto the MODIS grid
            clh = meanNoNan(clhcalipso_monthly,1)/100;
            clm = meanNoNan(clmcalipso_monthly,1)/100;
            cll = meanNoNan(cllcalipso_monthly,1)/100;            
            switch filt_dat_choose
                case 'calipso highCF'
                    filt_dat = clh;
                case 'calipso mid+highCF'
%                    filt_dat = clm + clh;
  %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                    filt_dat = clm + clh - (clm.*clh); 
                case 'calipso lowCF'
                    filt_dat = cll;                    
            end
            
            filt_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,filt_dat,MLON,MLAT');
            
             if ifilter_clhcalipso==1      
                 
                iclh = find(filt_dat<=thresh_clh(1) | filt_dat>thresh_clh(2));
                P(iclh) = NaN;
             end
             
             if iuse_saved_mask==1
                 saved_mask_filepath = ['~/' saved_mask_name '.mat'];
                 load(saved_mask_filepath);
                 eval(['isaved_mask =' saved_mask_name ';']);
                 P(isaved_mask) = NaN;
             end
             
             if iocean_only==1
                  %the landmask from CAMCLUBBv2.
                  load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
                 gcm_landmask_cal_MODIS_grid = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,gcm_landmask_cal,MLON,MLAT');
                 iland=find(gcm_landmask_cal_MODIS_grid>0.01 | isnan(gcm_landmask_cal_MODIS_grid)==1 );
                 P(iland)=NaN;
                 
                 MODIS_varname2_plot = [MODIS_varname2_plot ' ocean only'];

             end


            

             
             mean_bias_MODIS = meanNoNan(P(:),1);
             RMSE_MODIS = sqrt(meanNoNan(P(:).^2,1));
             
            icontour=1;
            
%            contour_label=0; %switch off contour labelling for this plot
%actually may as well leave the labelling on to make it clear to me what CF
%threshold is being used.
if icontour==1
    
    cont_ints=[0.3 0.3];
%    cont_ints=[0.35 0.35];    
    %N.B. - this is changed belwow some fields
    
    switch cont_dat_choose
        case 'calipso highCF'
            cont_dat = clh;
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'calipso mid+highCF'
            %                    filt_dat = clm + clh;
            %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
            cont_dat = clm + clh - (clm.*clh);
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'calipso lowCF'
            cont_dat = cll;
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'MOD35 CF'
            cont_dat = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
            cont_ints=[0.5:0.05:0.8];
        case 'MOD35 CF minus 1 std dev'
            [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
            cont_dat = mean_CF - std_CF;
            cont_ints=[0.1:0.05:0.8];        
        case 'MOD06 CF'
            cont_dat = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
            cont_ints=[0.5:0.05:0.8];
        case 'MOD06 CF minus 1 std dev'
            [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
            cont_dat = mean_CF - std_CF;
            cont_ints=[0.1:0.05:0.8];            
    end
            

            
            

    MODIS_varname2_plot = [MODIS_varname2_plot ' with ' cont_dat_choose ' contours'];
end

             

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

           

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            inew_cticks=0;

            switch error_type
                case 'absolute';
                    iset_min_clim=1;
                    clim_min=-50;
                    iset_max_clim=1;
                    clim_max=20;
                case  'percentage';
                    iset_min_clim=1;
                    clim_min=-50;
                    iset_max_clim=1;
                    clim_max=50;
            end
            
            
            
            
            
        case 'MODIS LWP minus AMSRE'
            %Here we are plotting MODIS in-cloud LWP that is multiplied by either the
            %MOD06 or MOD35 CF to give the grid-box mean LWP. Note that the
            %in-cloud LWP always comes from MOD06 because it is based upon
            %tau and reff. When multiplying by MOD35 we are assuming that
            %the LWP of the edge pixels that are detected by MOD35, but
            %ignored by MOD36 is the same as that for the non-edge pixels.
            
            %N.B. don't need to worry about selecting ocean only regions since AMSRE only
            %retrieves over the ocean.
            
            %May want to keep the option of doing the CALIPSO mid+high
            %cloud screening
            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                MODCF_type = 'MOD06 liquid';
                MODCF_type = 'MOD35';

                error_type = 'absolute';
                error_type = 'percentage';               
                
                iadd_correct_for_clear_sky_bias=1; %Subtract the clear-sky bias for AMSRE
                
                %what range of clouds to compare MODIS to (for bias calc)
                cal_CF_range='low';
                cal_CF_range='low+mid+high';     %Should reconsider using this as there may be cloud overlap
                  % Now have the cltcalipso field for daytime too, which I
                  % think should give the total column CF (i.e. including
                  % overlaps).
                
              %whether to filter our data where there is a the CALIPSO
              %high cloud CF is large (or mid+high - choose using cont_dat_choose below)
                ifilter_clhcalipso=0;
                thresh_clh = [-0.01 0.3];                               
%                thresh_clh = [0.3 1.01]; 

                zero_CF_out_of_range = 0; %flag that causes CFs to be set to zero if the CTH is determined to be outside
                % of the CTH range (rather than being ingnored as NaN.
                % Important when looking at high clouds since if CTH is low
                % then there is no high
                % cloud and so should count as CF=0 rather than ignoring it.
                % For low cloud is more dubious perhaps - except that this
                % is what CALIPSO would do since higher cloud would block
                % low cloud. Gives better match if use this for MOD35 vs CALIPSO comparison. 
                
                iocean_only=0; %only plot/calculate for the ocean only points.

% what to plot if plotting contours
                 cont_dat_choose = 'calipso highCF';
                 cont_dat_choose = 'calipso mid+highCF';
%                 cont_dat_choose = 'calipso lowCF';     
%                 cont_dat_choose = 'MOD06 CF';                 
%                 cont_dat_choose = 'MOD35 CF';  
%                 cont_dat_choose = 'MOD06 CF minus 1 std dev';
%                 cont_dat_choose = 'MOD35 CF minus 1 std dev';                 
                 
%What to use for filtering                 
                 filt_dat_choose = 'calipso mid+highCF';
                 
                 iuse_saved_mask=0;
                 saved_mask_name = 'imask_ndays15_cal_midhigh';
                 %save the .mat files so that this works:-
%                 saved_mask_filepath = ['~/' saved_mask_name '.mat'];

            end
            
            %the landmask from CAMCLUBBv2.
            load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
            switch error_type
                case 'absolute'
                    varname_str = ' minus AMSRE LWP (g m^{-2}) ';
                case 'percentage'
                    varname_str = ' vs AMSRE LWP percentage difference (%) ';                    
            end
            
            MODIS_varname2_plot = [MODCF_type varname_str aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str];

            colormap_choose=jet; %default

            %            dat_modis=H_time3;
            
% -- AMSRE and MODIS data should already be on the same grid, so no need for any interpolation       
         
switch MODCF_type
    case 'MOD06 liquid'
            CFdat_modis=Cloud_Fraction_Liquid.timeseries3; modis_data_plot=[modis_data_plot];
    case 'MOD35'
            CFdat_modis=Cloud_Fraction_Day_Mean.timeseries3; modis_data_plot=[modis_data_plot];        
end

%lwp_amsre_bias_comp = 0.5*(lwp_amsre_time3(:,:,:,1) + lwp_amsre_time3(:,:,:,2)); %an experiment
% Toniazzo Fig. 13 compares MODIS to an AMSRE daytime+nighttime average

modis_corr_factor = 0.85;


if abs(modis_corr_factor-1.0)>0.01
    MODIS_varname2_plot = [MODIS_varname2_plot ' MODIS CORRECTED BY f=' num2str(modis_corr_factor)];
end



lwp_amsre_bias_comp = lwp_amsre_time3(:,:,:,1); %daytime only overpasses
            dat_modis = 1000*(W_time3.*CFdat_modis*modis_corr_factor - lwp_amsre_bias_comp); %LWP*CF gives the grid-box average LWP
%            dat_modis = 1000*(lwp_amsre_time3(:,:,:,1)); %LWP*CF gives the grid-box average LWP
            %for AMSRE daytime is column 1, *1000 to convert to g/m2
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average                                                              
            dat_modis = dat_modis(:,:,time_inds_average);
            
            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)
            
            if iadd_correct_for_clear_sky_bias==1
                filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region'];  %7 percent refers to the threshold of MODIS CF to


                
                 filename = ['AMSRE_clear_sky_bias_9percent_VOCALS_region_20150918T065119'];
%                filename = ['AMSRE_clear_sky_bias_7percent_VOCALS_region_20150918T064749'];                         
%                filename = ['AMSRE_clear_sky_bias_5percent_VOCALS_region_20150918T065504'];        
%                filename = ['AMSRE_clear_sky_bias_3percent_VOCALS_region_20150918T065607'];
                 filename = ['AMSRE_clear_sky_bias_extrapolated_VOCALS_region_20150923T064911'];
                
                clear P_save_AMSRE_clear_sky
                
                
                filedir_savevars = '/home/disk/eos1/d.grosvenor/mat_files_various/';
                filename_savevars = [filedir_savevars filename '.mat'];
                load(filename_savevars);
                tag = '_AMSRE_clear_sky';
                
                %Try limiting to the modal value(=9 g/m2) due to
                %uncertainty about how the CF threshold affects the result
%                P_save_AMSRE_clear_sky(P_save_AMSRE_clear_sky>9) = 9;
%                P_save_AMSRE_clear_sky(P_save_AMSRE_clear_sky<0) = 0;                
                
                P = P + P_save_AMSRE_clear_sky;  %is in g/m2
                   %add teh clear-sky bias since are doing MODIS minus
                   %AMSRE and a positive AMSRE clear-sky bias means that
                   %AMSRE is too high, so we need to subtract this from
                   %AMSRE and therefore add to MODIS-AMSRE
                
                MODIS_varname2_plot = [MODIS_varname2_plot ' clear-sky correction ' filename ' '];
                
                
            end
            
            switch error_type
                case 'percentage'
                    P = 100* P ./ meanNoNan(1000*lwp_amsre_time3(:,:,:,1),3);
            end
            
            


            
% -- Use CALIPSO mid+high CF to define a certain region
% But need to grid the CALIPSO data onto the MODIS grid
            clh = meanNoNan(clhcalipso_monthly,1)/100;
            clm = meanNoNan(clmcalipso_monthly,1)/100;
            cll = meanNoNan(cllcalipso_monthly,1)/100;            
            switch filt_dat_choose
                case 'calipso highCF'
                    filt_dat = clh;
                case 'calipso mid+highCF'
%                    filt_dat = clm + clh;
  %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
                    filt_dat = clm + clh - (clm.*clh); 
                case 'calipso lowCF'
                    filt_dat = cll;                    
            end
            
            filt_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,filt_dat,MLON,MLAT');
            
             if ifilter_clhcalipso==1      
                 
                iclh = find(filt_dat<=thresh_clh(1) | filt_dat>thresh_clh(2));
                P(iclh) = NaN;
             end
             
             if iuse_saved_mask==1
                 saved_mask_filepath = ['~/' saved_mask_name '.mat'];
                 load(saved_mask_filepath);
                 eval(['isaved_mask =' saved_mask_name ';']);
                 P(isaved_mask) = NaN;
             end
             
             if iocean_only==1
                  %the landmask from CAMCLUBBv2.
                  load('~/CAMCLUBBv2_landmask_calgrid.mat','gcm_landmask_cal');
            
                 gcm_landmask_cal_MODIS_grid = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,gcm_landmask_cal,MLON,MLAT');
                 iland=find(gcm_landmask_cal_MODIS_grid>0.01 | isnan(gcm_landmask_cal_MODIS_grid)==1 );
                 P(iland)=NaN;
                 
                 MODIS_varname2_plot = [MODIS_varname2_plot ' ocean only'];

             end


            

             
             mean_bias_MODIS = meanNoNan(P(:),1);
             RMSE_MODIS = sqrt(meanNoNan(P(:).^2,1));
             
%             P=Npoints; %If want to look at how many samples we have
             
            icontour=1;
            
%            contour_label=0; %switch off contour labelling for this plot
%actually may as well leave the labelling on to make it clear to me what CF
%threshold is being used.
if icontour==1
    
    cont_ints=[0.3 0.3];
%    cont_ints=[0.35 0.35];    
    %N.B. - this is changed belwow some fields
    
    switch cont_dat_choose
        case 'calipso highCF'
            cont_dat = clh;
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'calipso mid+highCF'
            %                    filt_dat = clm + clh;
            %Prob better to use the random overlap assumption  CF = CF1+CF2-CF1*CF2
            cont_dat = clm + clh - (clm.*clh);
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'calipso lowCF'
            cont_dat = cll;
            cont_dat = griddata(gcm_Plon2D_CALIPSO_monthly,gcm_Plat2D_CALIPSO_monthly,cont_dat,MLON,MLAT');
        case 'MOD35 CF'
            cont_dat = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
            cont_ints=[0.5:0.05:0.8];
        case 'MOD35 CF minus 1 std dev'
            [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Day_Mean.timeseries3,3);
            cont_dat = mean_CF - std_CF;
            cont_ints=[0.1:0.05:0.8];        
        case 'MOD06 CF'
            cont_dat = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
            cont_ints=[0.5:0.05:0.8];
        case 'MOD06 CF minus 1 std dev'
            [mean_CF,nnums,std_CF] = meanNoNan(Cloud_Fraction_Liquid.timeseries3,3);
            cont_dat = mean_CF - std_CF;
            cont_ints=[0.1:0.05:0.8];            
    end
            

            
            

    MODIS_varname2_plot = [MODIS_varname2_plot ' with ' cont_dat_choose ' contours'];
end

             

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

           

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            inew_cticks=0;

            switch error_type
                case 'absolute';
                    iset_min_clim=1;
                    clim_min=-50;
                    iset_max_clim=1;
                    clim_max=20;
                case  'percentage';
                    iset_min_clim=1;
                    clim_min=-50;
                    iset_max_clim=1;
                    clim_max=50;
            end
            
            

            
            
            
        case 'Liquid Cloud Fraction time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=Cloud_Fraction_Liquid.timeseries3; modis_data_plot=[modis_data_plot];
            

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            dat_modis2 = dat_modis(:,:,time_inds_average);
            [P,Npoints] = meanNoNan(dat_modis2,3); %time mean (all times)

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' for ' thresh_str];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;
            
         case 'Ice Cloud Fraction time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=Cloud_Fraction_Ice.timeseries3; modis_data_plot=[modis_data_plot];
            

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            dat_modis2 = dat_modis(:,:,time_inds_average);
            [P,Npoints] = meanNoNan(dat_modis2,3); %time mean (all times)

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;    
            
  case 'Combined Cloud Fraction time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=Cloud_Fraction_Combined.timeseries3; modis_data_plot=[modis_data_plot];
            

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            dat_modis2 = dat_modis(:,:,time_inds_average);
            [P,Npoints] = meanNoNan(dat_modis2,3); %time mean (all times)

%            P=nswath;
%P = min(Solar_Zenith_Mean.timeseries3,[],3);

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;    
                        
            
        case 'Number of droplets cell values Daily from mockL3 time mean - specific days'
             colormap_choose=jet; %default

            %            dat_modis=H_time3;
            [dat_modis]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean_Daily.timeseries3,Cloud_Effective_Radius_Liquid_Mean_Daily.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean_Daily.timeseries3,'N');        
%            dat_modis=Cloud_Fraction_Liquid_L32007_Arctic.timeseries3; modis_data_plot=[modis_data_plot];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=300;         
            
        case 'Number of droplets cell values from L3 portion time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;
            [dat_modis]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean_L32007_Arctic.timeseries3,Cloud_Effective_Radius_Liquid_Mean_L32007_Arctic.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean_L32007_Arctic.timeseries3,'N');        
%            dat_modis=Cloud_Fraction_Liquid_L32007_Arctic.timeseries3; modis_data_plot=[modis_data_plot];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=300;            



        case 'Liquid Cloud Fraction Daily from L3 portion time mean - specific days'
             colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=Cloud_Fraction_Liquid_L32007_Arctic.timeseries3; modis_data_plot=[modis_data_plot];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;            


            
        case 'Liquid Cloud Fraction Daily from mockL3 time mean - specific days'
            colormap_choose=jet; %default

            %            dat_modis=H_time3;

            dat_modis=Cloud_Fraction_Liquid_Daily.timeseries3; modis_data_plot=[modis_data_plot];

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1;            


        case 'LWP cell values time mean (grid-box mean) - specific days'
            colormap_choose=jet; %default

            %are multiplying by CF to give the grid-box vlaue rather than
            %the in-cloud only value of LWP - easier for comparison to
            %model
%            dat_modis=1e3*W_time3.*Cloud_Fraction_Liquid.timeseries3*(1/1.15);  modis_data_plot=[modis_data_plot ' Painemal '];
            dat_modis=1e3*W_time3.*Cloud_Fraction_Liquid.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-2}';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=225;
            clim_max=150;            

            
            
      case 'LWP cell values time mean (grid-box mean using MOD35) - specific days'
            plot_num_datapoints=0;  %flag that plots the number of good data points
%            plot_num_datapoints=1; %instead of the actual data           
          
          
            colormap_choose=jet; %default

            %are multiplying by CF to give the grid-box vlaue rather than
            %the in-cloud only value of LWP - easier for comparison to
            %model
%            dat_modis=1e3*W_time3.*Cloud_Fraction_Day_Mean.timeseries3*(1/1.15);  modis_data_plot=[modis_data_plot ' Painemal '];
            dat_modis=1e3*W_time3.*Cloud_Fraction_Day_Mean.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
             Ndays2 = Npoints;
             
             if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-2}';
            title_info = '';

            icontour=0;
            cont_dat = Npoints;
            cont_ints=[0:5:100];

            inew_cticks=0;

            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                ifilter_ndays=1;

                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=225;
                clim_max=150;

            end
            
      case 'LWP cell values time mean (in-cloud avearge) - specific days'
            plot_num_datapoints=0;  %flag that plots the number of good data points
%            plot_num_datapoints=1; %instead of the actual data           
          
          
            colormap_choose=jet; %default

            %are multiplying by CF to give the grid-box vlaue rather than
            %the in-cloud only value of LWP - easier for comparison to
            %model
%            dat_modis=1e3*W_time3.*Cloud_Fraction_Day_Mean.timeseries3*(1/1.15);  modis_data_plot=[modis_data_plot ' Painemal '];
            dat_modis=1e3*W_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
             Ndays2 = Npoints;
             
             if plot_num_datapoints==1
                 P=Ndays2;
                 MODIS_varname2_plot = ['Number of datapoints for MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             else
                 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str 'for  ' thresh_str];
             end

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-2}';
            title_info = '';

            icontour=0;
            cont_dat = Npoints;
            cont_ints=[0:5:100];

            inew_cticks=0;

            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                ifilter_ndays=1;

                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=225;
                clim_max=150;

            end
            


            
        case 'LWC max time mean - specific days'
            colormap_choose=jet; %default

            %are multiplying by CF to give the grid-box vlaue rather than
            %the in-cloud only value of LWP - easier for comparison to
            %model
            %            dat_modis=1e3*H_time3.*cw_time3; % b/sqrt(b)=sqrt(b);
            dat_modis=1e3*H_time3.*cw_time3 * sqrt(0.7) * sqrt(1/1.15); % b/sqrt(b)=sqrt(b);
            %if cw is *0.8 then H=sqrt(2W/cw) = * 1/sqrt(cw). L=cw*H so
            %L = * 0.8/sqrt(0.8)= *sqrt(0.8). W will remain constant as is just a
            %function of tau and Re. Although if Re is too high then W will
            %be too

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-3}';
            title_info = '';

            icontour=1;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;


        case 'Cloud Top Temperature Minimum time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=Cloud_Top_Temperature_Day_Minimum.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='K';
            title_info = '';


            inew_cticks=0;


            iset_min_clim=1;
            clim_min=195;
            iset_max_clim=1;
            clim_max=295;

        case 'Cloud Top Pressure Minimum time max - specific days'
            colormap_choose=jet; %default

            dat_modis=Cloud_Top_Pressure_Day_Minimum.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %                        [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            P = max(dat_modis(:,:,time_inds_average),[],3); %time mean (all times)

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='hPa';
            title_info = '';


            inew_cticks=0;


            iset_min_clim=1;
            clim_min=300;
            iset_max_clim=1;
            clim_max=1000;

        case 'Cloud Top Pressure Minimum time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=Cloud_Top_Pressure_Day_Minimum.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            %                        P = max(dat_modis(:,:,time_inds_average),[],3); %time mean (all times)

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='hPa';
            title_info = '';


            inew_cticks=0;


            iset_min_clim=0;
            clim_min=300;
            iset_max_clim=0;
            clim_max=1000;

        case 'Cloud Top Pressure Standard Deviation time mean - specific days'
            colormap_choose=jet; %default

            dat_modis=Cloud_Top_Pressure_Day_Standard_Deviation.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            %                        P = max(dat_modis(:,:,time_inds_average),[],3); %time mean (all times)

            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='hPa';
            title_info = '';


            inew_cticks=0;


            iset_min_clim=1;
            clim_min=300;
            iset_max_clim=1;
            clim_max=1000;


        case 'BKUP - Number of Droplets Joint Histogram mean from selected days'
            icolormap_cf_grey=0;

            dat_modis=Nd_timeseries.mean;


            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later
            ihtot_other = ihtot; %store the ihtot values for the things
            %other than CF

            if icolormap_cf_grey==1
                %points where just the CF is bad
                ihtot_cf = find( ~ ( Cloud_Fraction_Liquid.timeseries3>=thresh_CF) );

                %the points where either CF or the others are bad (full screened set)
                ihtot = union(ihtot_cf,ihtot); %union combines with no repetitions
            end

            %remove points we defo don't want for Nd calc (full
            %screening)
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %do the mean over all the required times and for
            %both satellites if selected (Nd mean)
            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

            %and now do CF screening stuff - In case where one satellite
            %gives a high CF and so a valid Nd, but the
            %other does not we want to use the Nd from the good satellite
            %we only do CF if both satellites are bad

            if icolormap_cf_grey==1
                %make a NaN array - this will only include the low
                %CF points
                dat_modis2 = NaN*ones(size(dat_modis));
                %scale the cf to match the CF gray colorbar added
                %on the end
                dat_modis2(ihtot_cf) = 1e12+1e21*Cloud_Fraction_Liquid.timeseries3(ihtot_cf);
                %remove the points that are bad for other reasons
                %e.g. don't want to include the CF for a bad sensor
                %ZA
                dat_modis2(ihtot_other)=NaN;

                %new P average - only for low CF points, that are not screened
                %for other reasons. everything else is NaN & will be ignored
                [P2,Npoints2] = meanNoNan(dat_modis2(:,:,time_inds_average),3); %time mean (all times)

                %we only want to replace points that were screened
                %before - i.e. not good points where we have
                %droplet numbers from one or both satellites
                inanP = find(isnan(P)==1);
                P(inanP)=P2(inanP);

            end



            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram Nd '];

            units_str_plot='cm^{-3}';
            title_info = '';


            inew_cticks=1;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;



            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end


            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';




        case 'Number of Droplets Joint Histogram mean from selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Nd_timeseries.mean;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='cm^{-3}';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';


            MODIS_varname2_plot = [aqua_terra_str ' ' time_mean_str ', Joint Histogram Nd '];


            inew_cticks=1;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands


            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end


        case 'Mean Sensor ZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Sensor_Zenith_Mean.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [modis_data_plot ' ' time_mean_str];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end




        case 'Mean SZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Solar_Zenith_Mean.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' time_mean_str ', Mean Solar ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end
            
            
        case 'Min Sensor ZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Sensor_Zenith_Minimum.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS ' aqua_terra_str ' ' time_mean_str ', Min Sensor ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands

            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end
    
            
            
        case 'Max Sensor ZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Sensor_Zenith_Maximum.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS ' aqua_terra_str ' ' time_mean_str ', Max Sensor ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands

            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end
            

        case 'Max SZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Solar_Zenith_Maximum.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS ' aqua_terra_str ' ' time_mean_str ', Max Solar ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end
            
            
        case 'Max - Min SZA time mean, selected days'    
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Solar_Zenith_Maximum.timeseries3 - Solar_Zenith_Minimum.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS ' aqua_terra_str ' ' time_mean_str ','  modis_data_plot];


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=25;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end

        case 'Min SZA time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Solar_Zenith_Minimum.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' time_mean_str ', Min Solar ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end
            
            
        case 'Cloud Fraction Liquid Npixels time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Cloud_Fraction_Liquid_Pixel_Counts.timeseries3./Cloud_Fraction_Liquid.timeseries3;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' time_mean_str ', Max Solar ZA '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands




            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end    



        case 'Number of good data days for droplet number from cell values time mean, selected days'
            %set the screening for ihtot to NOT include CF
            %we make up the full screening set later

            icolormap_cf_grey=0;
            dat_modis=Nd_timeseries.mean;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Nd (cm^{-3})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' time_mean_str ', Joint Histogram Ngood days '];


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands

            P=Npoints; %only use the number of non-NaN days that go into the average


            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end


        case 'Number of good swaths for droplet number from mock L3, selected days'
            colormap_choose=jet; %default

            dat_modis=N_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average

            [mean_temp,P] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=90;
            
        case 'Number of different days with good Nd from mock L3, selected days'
            colormap_choose=jet; %default

            %find the number of datapoints where we had an Nd retrieval
            dat_modis=N_time3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            [mean_temp,P] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            
            %find the individual components of the date
            [Y,MO,D,H,MI,S] = datevec(Date_Time_Swath.timeseries3(:));
            %Just use the day (not hours and mins)
            all_days = datenum(Y,MO,D);   
            all_days = reshape(all_days,size(N_time3));
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
            all_days(ihtot) = NaN;
            all_days = all_days(:,:,time_inds_average);
            %remove the days for the points that are being screened
            
            
            a=zeros(size(all_days));

            %loop through all of the days
            for ii=1:size(all_days,1)
                for jj=1:size(all_days,2)
                    dat_ij = all_days(ii,jj,time_inds_average);
                    inot_nan = find(isnan(dat_ij)==0);

                    [B,I,J]=unique(dat_ij(inot_nan)); %unique doesn't work for NaNs - all are still included
                    a(ii,jj,time_inds_average(inot_nan(I)))=1;
                end
            end

            Ndays2=sum(a(:,:,time_inds_average),3);
                
                

            %Divide by the number of different days with data to get no.
            %swaths per day (with good Nd retrievals)
            P=Ndays2;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;



            

        case 'Mean number of swaths with good Nd per day from mock L3, selected days'
            colormap_choose=jet; %default

            %find the number of datapoints where we had an Nd retrieval
            dat_modis=N_time3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            [mean_temp,P] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            
            %find the individual components of the date
            [Y,MO,D,H,MI,S] = datevec(Date_Time_Swath.timeseries3(:));
            %Just use the day (not hours and mins)
            all_days = datenum(Y,MO,D);   
            all_days = reshape(all_days,size(N_time3));
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
            all_days(ihtot) = NaN;
            all_days = all_days(:,:,time_inds_average);
            %remove the days for the points that are being screened
            
            
            a=zeros(size(all_days));

            %loop through all of the days
            for ii=1:size(all_days,1)
                for jj=1:size(all_days,2)
                    dat_ij = all_days(ii,jj,time_inds_average);
                    inot_nan = find(isnan(dat_ij)==0);

                    [B,I,J]=unique(dat_ij(inot_nan)); %unique doesn't work for NaNs - all are still included
                    a(ii,jj,time_inds_average(inot_nan(I)))=1;
                end
            end

            Ndays2=sum(a(:,:,time_inds_average),3);
                
                

            %Divide by the number of different days with data to get no.
            %swaths per day (with good Nd retrievals)
            P=P./Ndays2;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;
            
            
        case 'Max number of swaths in one day from mock L3, selected days'
            colormap_choose=jet; %default

            %find the number of datapoints where we had an Nd retrieval
%            dat_modis=N_time3;
            %or use SZA for number of datapoints that had a daylight
            %overpass - if no overpass then will be NaN
            dat_modis=Solar_Zenith_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            [mean_temp,P] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)
            
            %find the individual components of the date
            [Y,MO,D,H,MI,S] = datevec(Date_Time_Swath.timeseries3(:));
            %Just use the day (not hours and mins)
            all_days = datenum(Y,MO,D);   
            all_days = reshape(all_days,size(N_time3));
            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use
            %time_inds_average
           %remove the days for the points that are being screened
%            all_days(ihtot) = NaN;
            all_days(isnan(dat_modis))=NaN; %to also remove the points that were originally NaN
            all_days = all_days(:,:,time_inds_average);

            
            
            a=zeros(size(all_days));
            max_swaths=zeros(size(all_days(:,:,1)));

            %loop through all of the days
            for ii=1:size(all_days,1)
                for jj=1:size(all_days,2)
                    dat_ij = all_days(ii,jj,time_inds_average);
                    inot_nan = find(isnan(dat_ij)==0);

                    [B,I,J]=unique(dat_ij(inot_nan)); %unique doesn't work for NaNs - all are still included
                    a(ii,jj,time_inds_average(inot_nan(I)))=1;
                    [mode_day,nmode_freq] = mode(J); %this gets the most frequent day and the number of occurences
                    %so we can know the max no. swaths per day
                    max_swaths(ii,jj) = nmode_freq;
                end
            end

           % Ndays2=sum(a(:,:,time_inds_average),3);
                
                

            
            P=max_swaths;
            P(P<1.9)=NaN;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;    
            
            
        case 'Min of the max daily SZA over period from mock L3, selected days'
            colormap_choose=jet; %default

        
            dat_modis=Solar_Zenith_Maximum_Daily.timeseries3; %made from process_mockL3_swaths_into_daily2.m
            %make NaN values for the screening (CF etc)
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            P = min(dat_modis(:,:,time_inds_average),[],3); %time minimum over all days
            %to get an idea of the best that can be hoped for for each
            %location
            
            
%            P(P>60)=NaN;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;       
            
            
 case 'Max of the max daily SZA over period from mock L3, selected days'
            colormap_choose=jet; %default

        
            dat_modis=Solar_Zenith_Maximum_Daily.timeseries3; %made from process_mockL3_swaths_into_daily2.m
            %make NaN values for the screening (CF etc)
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            P = max(dat_modis(:,:,time_inds_average),[],3); %time maximum over all days
            %to get an idea of the best that can be hoped for for each
            %location
            
            
%            P(P>60)=NaN;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;       
            
            
        case 'Max of the min daily SZA over period from mock L3, selected days'
            colormap_choose=jet; %default

        
            dat_modis=Solar_Zenith_Minimum_Daily.timeseries3; %made from process_mockL3_swaths_into_daily2.m
            %make NaN values for the screening (CF etc)
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            P = max(dat_modis(:,:,time_inds_average),[],3); %time maximum over all days
            %to get an idea of the best that can be hoped for for each
            %location
            
            
%            P(P>60)=NaN;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;           
            
        case 'Min of the min daily SZA over period from mock L3, selected days'
            colormap_choose=jet; %default

        
            dat_modis=Solar_Zenith_Minimum_Daily.timeseries3; %made from process_mockL3_swaths_into_daily2.m
            %make NaN values for the screening (CF etc)
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            P = min(dat_modis(:,:,time_inds_average),[],3); %time maximum over all days
            %to get an idea of the best that can be hoped for for each
            %location
            
            
%            P(P>60)=NaN;
            
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            %                        MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';

            icontour=0;

            inew_cticks=0;


            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=225;                       





        case 'Cloud Depth Joint Histogram mean from selected days'
            icolormap_cf_grey=0;
            dat_modis=H_timeseries.mean;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Cloud Depth (m)'
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='m';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram Cloud Depth '];


            inew_cticks=1;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands


            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end




        case 'Cloud Fraction mean from selected days'
            dat_modis=Cloud_Fraction_Liquid.timeseries3;

            %remove points we  don't want from screening
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            %                         ihtot_cf = find( ~ ( Cloud_Fraction_Liquid.timeseries3<thresh_CF) );
            %                         dat_modis(ihtot_cf)=NaN;



            %do the mean over all the required times and for
            %both satellites if selected
            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)

            %                         temp=NaN*ones(size(P));
            %                         temp(inanP)=P(inanP);
            %                         P=temp;

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS TERRA ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram Nd '];

            units_str_plot='';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;







        case 'LWP Joint Histogram mean from selected days'
            icolormap_cf_grey=0;
            dat_modis=1000*W_timeseries.mean;

            cb_lab_pos_left=0.2;
            cb_lab_str_left='Liquid Water Path (g m^{-2})';
            cb_lab_pos_right=0.85;
            cb_lab_str_right='Cloud Fraction';

            units_str_plot='g m^{-2}';
            title_info = '';

            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = [aqua_terra_str ' ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram LWP '];


            inew_cticks=1;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;

            %execute this external script
            CF_grayscale_commands


            if icolormap_cf_grey==1
                ctick_range_choose
                colormap_CF_jet_with_grey
                colormap_choose = cf_colormap;
            end



        case 'Cloud Top Temperature timeseries3 mean from selected days'
            dat_modis=Cloud_Top_Temperature_Day_Mean.timeseries3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            units_str_plot='m';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=280;
            iset_max_clim=1;
            clim_max=292;

            
        case 'Cloud Optical Depth timeseries3 mean from selected days'
            dat_modis=Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
 MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];

            units_str_plot='';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;
            
        case 'Cloud Effective Radius 1.6\mum timeseries3 mean from selected days'
%N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7' in pdf2d_plot_commands)
            match_POLDER=1;
  
            dat_modis=Cloud_Effective_Radius_16_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
            if match_POLDER==1
                ipol = find(isnan(Par2_CDR_coloc)==1);
                dat_modis(ipol)=NaN;
            end

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=6;
            iset_max_clim=1;
            clim_max=25;  
            
        case 'Cloud Effective Radius 2.1\mum timeseries3 mean from selected days'
%N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6 Re_3.7' in pdf2d_plot_commands)
  
            dat_modis=Cloud_Effective_Radius_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            match_POLDER=1;
            if match_POLDER==1
                ipol = find(isnan(Par2_CDR_coloc)==1);
                dat_modis(ipol)=NaN;
            end
            
            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=6;
            iset_max_clim=1;
            clim_max=25;  
            
            
        case 'Cloud Effective Radius 3.7\mum timeseries3 mean from selected days'
            %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6
  %Re_3.7' in pdf2d_plot_commands)
  
            dat_modis=Cloud_Effective_Radius_37_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            match_POLDER=1;
            if match_POLDER==1
                ipol = find(isnan(Par2_CDR_coloc)==1);
                dat_modis(ipol)=NaN;
            end
            
            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=6;
            iset_max_clim=1;
            clim_max=25;  
            
        case 'Cloud Effective Radius 3.7\mum minus 2.1\mum timeseries3 mean from selected days'
            %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6
  %Re_3.7' in pdf2d_plot_commands)
  
            dat_modis=Cloud_Effective_Radius_37_Liquid_Mean.timeseries3 - Cloud_Effective_Radius_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;  
            
        case 'Cloud Effective Radius 1.6\mum minus 2.1\mum timeseries3 mean from selected days'
            %N.B. - these don't have the datapoints where there are not retreivals
  %from all 3 removed (see 'R_{eff 3.7 \mum} (\mum) reduced dataset Re_1.6
  %Re_3.7' in pdf2d_plot_commands)
  
            dat_modis=Cloud_Effective_Radius_16_Liquid_Mean.timeseries3 - Cloud_Effective_Radius_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;  
            
       case 'Re3.7 minus re2.1 mockL3'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3 - Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;             
            
        case 'Re1.6 minus re2.1 mockL3'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3 - Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2; 
            
            
        case 'Mock L3 Re2.1 minus L3'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3 - Cloud_Effective_Radius_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;    
            
        case 'Mock L3 Tau minus L3'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Optical_Thickness_Liquid_Mean_mockL3.timeseries3 -  Cloud_Optical_Thickness_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=-5;
            iset_max_clim=0;
            clim_max=2;    
            
        case 'Re1.6 mockL3 minus POLDER'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3 - daymean_Par2_CDR;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;    
            
            
       case 'Re2.1 mockL3 minus POLDER'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3 - daymean_Par2_CDR;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2;  
            
         case 'Re3.7 mockL3 minus POLDER'
           %From one year of 2007 VOCALS data for mockL3
  
            dat_modis= Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3 - daymean_Par2_CDR;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=-5;
            iset_max_clim=1;
            clim_max=2; 
            
  
            
            
            
        case 'Re3.7 mockL3'
           %From one year of 2007 VOCALS data for mockL3
           
           if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
              ico_locate_POLDER=0;
           end
  
            dat_modis= Cloud_Effective_Radius_37_Liquid_Mean_mockL3.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
            if ico_locate_POLDER==1
               dat_modis(isnan( daymean_Par2_CDR ))=NaN;
            end

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=5;
            iset_max_clim=1;
            clim_max=30; 
            
           case 'Re2.1 mockL3'
           %From one year of 2007 VOCALS data for mockL3
           
           if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
              ico_locate_POLDER=0;
           end
  
  
            dat_modis = Cloud_Effective_Radius_Liquid_Mean_mockL3.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
            
            if ico_locate_POLDER==1
               dat_modis(isnan( daymean_Par2_CDR ))=NaN;
            end 
            

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=5;
            iset_max_clim=1;
            clim_max=30;    
              
        case 'Re1.6 mockL3'
           %From one year of 2007 VOCALS data for mockL3
           
           if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
               ico_locate_POLDER=0;
           end
                        
  
            dat_modis= Cloud_Effective_Radius_16_Liquid_Mean_mockL3.timeseries3;                        
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            
                        
            if ico_locate_POLDER==1
               dat_modis(isnan( daymean_Par2_CDR ))=NaN;
            end
            

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=1;
            clim_min=5;
            iset_max_clim=1;
            clim_max=30;           
            
        case 'Cloud Effective Radius QA pixel weighted timeseries3 mean from selected days'

            %Apparently monthly averages (M3) prdoucts are weighted
            %according the numner of pixels in a 1x1 box rather than giving
            %equal weight to each day.
            Npix = Cloud_Fraction_Liquid_Pixel_Counts.timeseries3;
            dat_modis=Cloud_Effective_Radius_Liquid_QA_Mean.timeseries3.*Npix;


            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average
            Npix(ihtot)=NaN;

            %need to make sure we are applying the time_inds_average AFTER
            %NaN screening since the screening does NOT use time_inds_average
            

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3,'sum'); %sum after multiplying by no. pixels
            Npix_tot = meanNoNan(Npix(:,:,time_inds_average),3,'sum'); %sum over time (of non-NaN values)
            Npoints = Npix_tot;
            P = P ./ Npix_tot; %Divide by total no. pixels to get the mean
            
            Ndays2 = Npoints;
            
            


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                inew_cticks=0;

                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=500;
            end
            
        case 'Cloud Effective Radius timeseries3 mean from selected days'
            dat_modis=Cloud_Effective_Radius_Liquid_Mean.timeseries3;
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            MODIS_varname2_plot = [modis_data_plot ' MODIS ' aqua_terra_str ' ' time_mean_str modisyears_str ' '];
 
            units_str_plot='\mum';
            title_info = '';


            if ~exist('ioverride_plotglobal_thresh') | ioverride_plotglobal_thresh==0
                inew_cticks=0;

                iset_min_clim=0;
                clim_min=0;
                iset_max_clim=0;
                clim_max=500;
            end


        case 'Cloud Depth timeseries3 mean from selected days'
            dat_modis=H_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS TERRA ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', timeseries3 Cloud Depth '];

            units_str_plot='m';
            title_info = '';


            inew_cticks=0;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=500;




        case 'OLD Cloud Depth Joint Histogram mean from selected days'
            dat_modis=H_timeseries.mean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS TERRA ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram Cloud Depth '];

            units_str_plot='m';
            title_info = '';


            inew_cticks=1;

            iset_min_clim=0;
            clim_min=0;
            iset_max_clim=0;
            clim_max=750;

        case 'Cloud Top LWC Joint Histogram mean from selected days'
            dat_modis=1000*cw_time3.*H_timeseries.mean;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis(:,:,time_inds_average),3); %time mean (all times)


            %                      MODIS_varname2_plot = modis_data_plot;
            modis_year_timeseries3='2008';
            MODIS_varname2_plot = ['MODIS TERRA ' datestr(days_required_for_mean(1)+datenum(['01-Jan-' modis_year_timeseries3])-1,'mmm-dd-yyyy') ', Joint Histogram Cloud Top LWC '];

            units_str_plot='g m^{-3}';
            title_info = '';


            inew_cticks=0; %for a non-linear colorbar

            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            clim_max=1.2;



        case 'LWC at cloud top'
            dat_modis = 1e3 * cw_time3.*H_time3; %convert to g/m3

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='g m^{-3}';
            title_info = '';


            inew_cticks=0;

        case 'Cloud Depth'
            dat_modis = H_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='m';
            title_info = '';


            inew_cticks=0;

        case 'Effective Radius'
            dat_modis = 1e6*reff_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='\mum';
            title_info = '';


            inew_cticks=0;

        case 'Optical Depth'
            dat_modis = tau_time3;

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [P,Npoints] = meanNoNan(dat_modis,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';


            inew_cticks=0;




        case 'Number of good data days for droplet number from cell values time mean'

            dat_modis=N_time3;
            %criteria is now set above
            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            [means,P] = meanNoNan(dat_modis,3); %time mean (all times)


            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='';
            title_info = '';
            
        case 'Fraction of days remaining after screening'
            %run this after doing a plot_global_map of the required data
            %with the required screening
%            dat_modis=N_time3;

            %Number of days that had a non-NaN cloud fraction
            Ndays_CF = zeros(size(Cloud_Fraction_Day_Mean.timeseries3(:,:,time_inds_average)));
            inotnan = find(isnan(Cloud_Fraction_Day_Mean.timeseries3(:,:,time_inds_average))==0);
            Ndays_CF(inotnan)=1; %allows us to count the number of possible data days for each location
            Ndays_poss = sum(Ndays_CF,3);
            Pold=P;
            P = Npoints./Ndays_poss;            

            
%             %NaN points that were excluded on a geographical basis in
%             %previous plot
%             inan=find(isnan(Pold)==1);
%             P(inan)=NaN;


            MODIS_varname2_plot = [modis_data_plot ' ' thresh_str];

            units_str_plot='';
            title_info = '';
            
            %contour of NaN points that were excluded on a geographical basis in
%             %previous plot
            icontour=1;
            cont_dat = zeros(size(Pold));
            inan=find(isnan(Pold)==1);
            cont_dat(inan)=1;
            cont_ints=[1 1];

            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            %        clim_max=225;
            clim_max=1;

case 'Fraction of days remaining after screening for AMSRE bias comp'
    %Added caveat that there has to be an AMSRE retrieval to count for the
    %possible no. days
            %run this after doing a plot_global_map of the required data
            %with the required screening
%            dat_modis=N_time3;

            %Number of days that had a non-NaN cloud fraction
            Ndays_CF = zeros(size(Cloud_Fraction_Day_Mean.timeseries3(:,:,time_inds_average)));
            inotnan = find(isnan(Cloud_Fraction_Day_Mean.timeseries3(:,:,time_inds_average))==0 & isnan(lwp_amsre_time3(:,:,time_inds_average,1))==0);
            Ndays_CF(inotnan)=1; %allows us to count the number of possible data days for each location
            Ndays_poss = sum(Ndays_CF,3);
            Pold=P;
            P = Npoints./Ndays_poss;            

            
%             %NaN points that were excluded on a geographical basis in
%             %previous plot
%             inan=find(isnan(Pold)==1);
%             P(inan)=NaN;


            MODIS_varname2_plot = [modis_data_plot ' ' thresh_str];

            units_str_plot='';
            title_info = '';
            
            %contour of NaN points that were excluded on a geographical basis in
%             %previous plot
            icontour=1;
            cont_dat = zeros(size(Pold));
            inan=find(isnan(Pold)==1);
            cont_dat(inan)=1;
            cont_ints=[1 1];

            iset_min_clim=1;
            clim_min=0;
            iset_max_clim=1;
            %        clim_max=225;
            clim_max=1;




        case 'Number of droplets histogram time mean'

            dat_modis=Nd_timeseries.mean;

            %cut some of the data out of the average
            ihtot = find( ~ ( Cloud_Fraction_Liquid_Pixel_Counts.timeseries3./Cloud_Fraction_Liquid.timeseries3>50 & Cloud_Fraction_Liquid.timeseries3>0.8 ) ); thresh_str='NP.GT.50 AND CF.GT.0.8';

            dat_modis(ihtot)=NaN; %make these values NaN as they will then be removed from the average

            P = meanNoNan(dat_modis,3); %time mean (all times)



            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            title_info = '';

        case 'Various'

            var_case ='SZA normalised std dev';
            %                        var_case ='SZA unnormalised std dev';
            %                        var_case ='SZA mean';
            %                        var_case='SZA Ndays';

            switch var_case
                case 'SZA Ndays'
                    P = Solar_Zenith_Mean.Ndata;
                case 'SZA normalised std dev'
                    P = Solar_Zenith_Mean.time_stdev./Solar_Zenith_Mean.data;


                case 'SZA unnormalised std dev'
                    P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.totNpix;
                    P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.totNpix ...
                        .* Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ndata / nMOD_av;
                    MODIS_varname2_plot = modis_data_plot; %string for the name

                    P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ndata;

                    P = Sensor_Zenith_Mean.data;

                    P = Solar_Zenith_Mean.time_stdev;

                    %                        P = Solar_Zenith_Mean.data;
                    %                        P = Sensor_Azimuth_Mean.data;

                case 'SZA mean'
                    P = Solar_Zenith_Mean.data;

            end


            MODIS_varname2_plot = var_case; %string for the name




            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Number from mode of 2D histogram - single day'
            N_H_calc_histo
            P = Nmode;

            MODIS_varname2_plot = modis_data_plot; %string for the name


            units_str_plot='cm ^{-3}';
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;


        case 'No. pixels'
            N_H_calc_histo
            %                        P = totN./cf;
            P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.totNpix;
            P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.totNpix ...
                .* Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ndata / nMOD_av;
            MODIS_varname2_plot = modis_data_plot; %string for the name

            P = Cloud_Optical_Thickness_Liquid_Joint_Histogram_vs_Effect_Radius.Ndata;

            P = Sensor_Zenith_Mean.data;
            P = Solar_Zenith_Mean.data;
            P = Sensor_Azimuth_Mean.data;

            MODIS_varname2_plot = modis_data_plot; %string for the name


            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Cloud Fraction Liq+Ice+Undetermined'
            P= Cloud_Fraction_Liquid.data + Cloud_Fraction_Ice.data + Cloud_Fraction_Undetermined.data;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Cloud Fraction Combined'
            P= Cloud_Fraction_Combined.data;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Cloud Fraction Liquid'
            P= Cloud_Fraction_Liquid.data;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Cloud Fraction Ice'
            P= Cloud_Fraction_Ice.data;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Cloud Fraction Undetermined'
            P= Cloud_Fraction_Undetermined.data;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot=Cloud_Fraction_Liquid.units_str;
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;

        case 'Number of Droplets'
            title_info = '';

            %                        MODIS_N_H_calc %runs the routine to calculate Nd
            %N calculation now done in filtering_data_get

            P=N;
            P(ihtot)=NaN;
            MODIS_varname2_plot = modis_data_plot;

            units_str_plot='cm^{-3}';
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;



            inew_cticks=0;

        case 'Number of Droplets Joint Histogram'
            N_H_calc_histo %runs the routine to calculate Nd
            P=N_histo_mean;
            MODIS_varname2_plot = 'N_d from tau,reff joint histogram';

            units_str_plot='cm^{-3}';
            modis_day_str_plot=Cloud_Fraction_Liquid.modis_day_str;
            modis_year_str_plot=Cloud_Fraction_Liquid.modis_year_str;





    end

%%  End of main switch


    %%%%   decides on the range for the colobar for new cticks
    ctick_range_choose
    %%%%  ****************************************************  %%%%%%

    figlab=[MODIS_varname2_plot ' ' proj_type];

    if inew_figure==1
        hf=figure('position',posit,'name',figlab);
    end

    %cut out points where only have a few points
    %included in the average from the display.
    if ifilter_ndays==1
        MODIS_varname2_plot = [MODIS_varname2_plot ' for ndays>' num2str(thresh_ndays) ' only '];
        
        switch mod_data_type
            case {'Monthly_data','GCM','GCM_COSP','CLOUDSAT-PRECIP','CALIPSO','Sea-ice','AMSRE'}
                %already calculated Ndays2 above
            otherwise
                %convert to "Matlab time" - days since 01-jan-0000 so that the
                %same days from different years are actually different
                daynums_matlab = datenum(modisyear_timeseries3,1,1) + daynum_timeseries3_MODIS - 1;
                all_days = repmat(repmat(daynums_matlab',[1 1 1]),[1 size(Cloud_Fraction_Liquid.timeseries3,1) size(Cloud_Fraction_Liquid.timeseries3,2)]);
                all_days = permute(all_days,[2 3 1]);
                %remove the days for the points that are being screened
%                all_days(ihtot) = NaN;
                all_days(isnan(dat_modis)==1)=NaN; %set all NaN points in the data to NaN

                a=zeros(size(all_days));

                %loop through all of the days
                for ii=1:size(all_days,1)
                    for jj=1:size(all_days,2)
                        dat_ij = all_days(ii,jj,time_inds_average);
                        inot_nan = find(isnan(dat_ij)==0);

                        [B,I,J]=unique(dat_ij(inot_nan)); %unique doesn't work for NaNs - all are still included
                        a(ii,jj,time_inds_average(inot_nan(I)))=1;
                    end
                end

                Ndays2=sum(a(:,:,time_inds_average),3);
                Ndays2_save_plot_global = Ndays2;

        end

%to filter by no. datapoints rather than Ndays2
%        Ndays2=Npoints; 
        
        P(  find( ~(Ndays2>=thresh_ndays) ) ) =NaN; thresh_str=[thresh_str ' AND ndays.GTE.' num2str(thresh_ndays)];

        %             P(find( ~(Npoints>=thresh_ndays)) )=NaN; thresh_str=[thresh_str ' AND Npoints.GTE.' num2str(thresh_ndays)];
    end

    if irestrict_domain==1
        if size(Plat2D) ~= size(P)
            error(['Size Plat2D = ' num2str(size(Plat2D)) ', size of P is ' num2str(size(P)) ' - need to be the same!']);
        end
        idomain = find( ~ ( Plat2D>=thresh_LAT(1) & Plat2D<thresh_LAT(2) & Plon2D>=thresh_LON(1) & Plon2D<thresh_LON(2) ) );
        P(idomain)=NaN;
        P_std(idomain)=NaN;
        
        if exist('dat_modis2')
        
        [domI,domJ]=ind2sub(size(Plat2D),idomain);

                %these are the lat/lon indices - need to repeat these for
                %each time and calculate a new index for the 3D matrix
                %lat/lon indices
                IJ = repmat(idomain,[1 length(time_inds_average)]);
                %K - replicate all of the time indices
                K = repmat([1:length(time_inds_average)],[length(idomain) 1]);
                iALL = sub2ind([size(dat_modis2,1)*size(dat_modis2,2) size(dat_modis2,3)] , IJ(:),K(:));
                dat_modis2(iALL)=NaN;
        end
        
        Pntotal_avail = length(P(:)) - length(idomain); %how many spatial points are within the domain

    else
        Pntotal_avail = NaN;
    end
    
    





%     switch mod_data_type
%         case {'timeseries3','daily','timeseries3 lambert','Monthly_data'}
%             dlat=1;
%             dlon=1;
% 
%             %Plat=90-0.25-[0:359];  %the example was for 0.5x0.5 degree data (360*720=lat*lon)
%             %Plon=-180+0.25+[0:719];
% 
%             %Plat=90-dlat/2-[0:180/dlat-1];  %this is generalised for any grid spacing
%             %Plon=-180+dlon/2+[0:360/dlon-1];
% 
% 
% 
%             dlat=abs(mean(diff(MLAT)));
%             dlon=abs(mean(diff(MLON)));
% 
% %            [Plon,Plat]=meshgrid(MLON-dlon/2,MLAT+dlat/2); %this actally cancels the +dlon/2 above!...?
%             
% %             [Plon,Plat]=meshgrid(MLON-dlon/2,MLAT-dlat/2);
%             
%             
%             
%     end

%    inot_cf=find(P<1e11); %was set this, but caused issues for number
%    concs where numbers do get >1e11 ! Was this set because of some cloud
%    fraction data somewhere?
    inot_cf=find(P<1e999);    
    P_linear=P(inot_cf);
    siz=size(P_linear);
    [mean_dim,idim]=max(siz);
    if siz(1)==0 %If all data is NaN
        Pmean=NaN;
        n=NaN;
        Pstd=NaN;
    else
        [Pmean,n,Pstd] = meanNoNan(P_linear,idim);
    end
    
    Pnpoints = n;
    P_RMSE = sqrt(meanNoNan(P_linear.^2,idim));

    fprintf(1,'\nOverall Mean = %f',Pmean);   

    P_save=P;
    
    if noplot==1
        return
    end


    if inew_cticks==1
        %N.B.  x_cbar_vals,cbar_vals chosen in ctick_range_choose.m (excuted earlier)
        P=map_data_onto_cticks(P,x_cbar_vals,cbar_vals);
    end



    switch mod_data_type

        case {'timeseries3','daily','Monthly_data_global','Sea-ice'}
            
 

            switch proj_type
                case 'polar'
                    %set the projection type to polar sterographic
                    %only seems to work well when set latitude to be the poles
                    %lon=98.0 makes the grid box square for the ecmwf_ml_nudgding runs
                    %since this is what stand_lon is set to in namelist.wps

                  

%the 'rad' value specified for the polar plot is the number of degrees outward from the pole to
%plot from the pole

                   if ~exist('ioverride_plotglobal_loc') | ioverride_plotglobal_loc==0    
                      lat_polar=-90; %Antarctic
                    lat_polar= 90; %Arctic

                      m_proj('stereographic','lat',lat_polar,'lon',-98,'rad',45);
                      m_proj('stereographic','lat',lat_polar,'lon',0,'rad',45);                      
                      m_proj('stereographic','lat',lat_polar,'lon',0,'rad',50);                                            
                    %  m_proj('stereographic','lat',lat_polar,'lon',-98,'rad',65);
%                    m_proj('stereographic','lat',lat_polar,'lon',0,'rad',35);
%                    m_proj('stereographic','lat',lat_polar,'lon',0,'rad',25);
                   else
                       eval(stereo_str01); %e.g. lat_polar=-90;
                       eval(stereo_str02); %e.g. m_proj('stereographic','lat',lat_polar,'lon',-98,'rad',45);                      
                   end

                    %xaxislocation property is important for south pole plots as otherwise
                    %it puts the longitude labels in the middle, which is messy
                    %setting to 'top' means they are at the northernmost part of the
                    %plot
                    %                            figure
                    %draw the coastlines - patch does it as a filled patch instead of a line
                    %                            m_coast('patch',[.6 .6 .6]);
                    %                            set(gca,'fontsize',24);


                    if i_dpcolor==1

                        m_dpcolor(Plon2D_edges-360,Plat2D_edges,P); shading flat; %colormap(map);
                        hold on
                        m_dpcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);

                    else
                        m_pcolor(Plon2D_edges-360,Plat2D_edges,P); shading flat; %colormap(map);
                        hold on
                        m_pcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);

                    end



                case 'global oval'

                    global_type = 'oval'; %oval shaped globe
%                    global_type = 'flat'; %rectangular shaped globe                    
                    switch global_type
                        case 'oval'
                            m_proj('hammer-aitoff','clongitude',0);
                            
                           

                        case 'flat'
                                    m_proj('miller','lat',82);                                   

                    end
                    
                    %                m_grid('box','fancy','tickdir','in');

                    if iseaice==1

                        iwest = find(abs(Plon2D_edges)>179);
                        Plat2D_edges(iwest)=NaN;
                        Plon2D_edges(iwest)=NaN;
                        P(iwest)=NaN;


                    end
                    
                    if i_dpcolor==1
                        m_dpcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);
                        hold on
                        m_dpcolor(Plon2D_edges-360,Plat2D_edges,P); shading flat; %colormap(map);
                       
                    else
%                        m_pcolor(Plon2D_edges-360,Plat2D_edges,P); shading flat; %colormap(map);
%                        hold on
%                        m_pcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);
                        
                        m_pcolor(Plon2D-360,Plat2D,P); shading flat; %colormap(map);
                        hold on
                        m_pcolor(Plon2D,Plat2D,P); shading flat; %colormap(map);
                        
                    end
                 

            end









            clims_get=get(gca,'clim');


            %draws some extra info in the form of contours
%            icontour=0;
            %contour_case ='Npoints';
            %contour_case ='Sea-ice';
%            contour_case ='Specified data';
            
            if icontour==1
                
                switch contour_case
                    case 'Sea-ice'
                        switch time_mean_str
                            case 'Jan'
                                seaice_days = {'0101','0201'};
                            case 'Feb'
                                seaice_days = {'0201','0301'};
                            case 'Mar'
                                seaice_days = {'0301','0401'}; 
                            case 'Apr'
                                seaice_days = {'0401','0501'}; 
                            case 'May'
                                seaice_days = {'0501','0601'}; 
                            case 'Jun'
                                seaice_days = {'0601','0701'}; 
                            case 'Jul'
                                seaice_days = {'0701','0801'}; 
                            case 'Aug'
                                seaice_days = {'0801','0901'}; 
                            case 'Sep'
                                seaice_days = {'0901','1001'}; 
                            case 'Oct'
                                seaice_days = {'1001','1101'};                                 
                            case 'Nov'
                                seaice_days = {'1101','1201'}; 
                            case 'Dec'
                                seaice_days = {'1201','1231'};     
                            case 'DJF'
                                seaice_days = {'1201','0301'};
                            case 'MAM'
                                seaice_days = {'0301','0601'};    
                            case 'JJA'
                                seaice_days = {'0601','0901'}; 
                            otherwise
                                if exist('multi_case') & strcmp(multi_case,'Fortnightly') ==1
                                    seaice_days = {seaice_start_str,seaice_end_str};
                                else
                                    error('*** This period is not specified for sea-ice!! ***');
                                    seaice_days={''};
                                end
                        end
                        
                        lwidth_seaice = 7;
                        lwidth_seaice = 3;                        
                            
                        cont_dat = eval(['seaice_data_' seaice_days{1} '*100;']); cont_ints=[0 0];                                                       
                        [cont_out,hcont]=m_contour(seaice_lon-360,seaice_lat,cont_dat,cont_ints,'r','linewidth',lwidth_seaice); shading flat; %colormap(map);
                        [cont_out,hcont]=m_contour(seaice_lon,seaice_lat,cont_dat,cont_ints,'r','linewidth',lwidth_seaice); shading flat; %colormap(map);
                        
                        cont_dat = eval(['seaice_data_' seaice_days{2} '*100;']); cont_ints=[0 0];
                        [cont_out,hcont]=m_contour(seaice_lon-360,seaice_lat,cont_dat,cont_ints,'m--','linewidth',lwidth_seaice); shading flat; %colormap(map);
                        [cont_out,hcont]=m_contour(seaice_lon,seaice_lat,cont_dat,cont_ints,'m--','linewidth',lwidth_seaice); shading flat; %colormap(map);
                                       
%                        clabel(cont_out,hcont,'fontweight','bold');
                    case 'Npoints'
                        cont_dat=Npoints;  cont_ints=[15 30 60];
                        %cont_dat=Ndays2;   cont_ints=[15 30 60];

                        [cont_out,hcont]=m_contour(Plon2D_edges-360,Plat2D_edges,cont_dat,cont_ints,'k'); shading flat; %colormap(map);
                        [cont_out,hcont]=m_contour(Plon2D_edges,Plat2D_edges,cont_dat,cont_ints,'k'); shading flat; %colormap(map);
                        clabel(cont_out,hcont,'fontweight','bold');
                        
                     case 'Specified data'
%                        cont_dat=Npoints;  cont_ints=[15 30 60];
                        %cont_dat=Ndays2;   cont_ints=[15 30 60];
                        if i_dpcolor==1
                            %N.B. - 'w' colour doesn't show up when saving for some reason?
                            %Use 'c' (cyan) instead
                            [cont_out,hcont]=m_dcontour(Plon2D_edges-360,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                            [cont_out,hcont]=m_dcontour(Plon2D_edges,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                        else
                            [cont_out,hcont]=m_contour(Plon2D_edges-360,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                            [cont_out,hcont]=m_contour(Plon2D_edges,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                        end
                        
                        clabel(cont_out,hcont,'fontweight','bold');
                        
                end
                
                set(hcont,'linewidth',3);
                set(gca,'clim',clims_get);


            end

            switch proj_type
                case 'polar'
                    if lat_polar<0
%                        m_grid('xaxislocation','top','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','c','xtick',[-180:10:180],'ytick',[-90:5:0]);   %white doesn't come out when saved for some reason so can use cyan ('c')
                        m_grid('xaxislocation','top','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','c','xtick',[-180:20:180],'ytick',[-90:10:0],'fontsize',6,'ticklen',0.005);   %white doesn't come out when saved for some reason so can use cyan ('c')                        
                    else
%                        m_grid('xaxislocation','bottom','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','k','xtick',[-180:10:180],'ytick',[0:2:90],'fontsize',6,'ticklen',0.005);
                        m_grid('xaxislocation','bottom','yaxislocation','middle','tickdir','out','linest','-','fontsize',14,'backcolor',[0 0 0],'color','k','xtick',[-180:10:180],'ytick',[0:5:90],'fontsize',6,'ticklen',0.005);                        
                    end
                    %N.B. m_ungrid removes a grid so that another
                    %can be drawn instead
                    %Also increase_font_size_map_figures may change
                    %fontsize later








                case 'global oval'                    
                    switch global_type
                        case 'oval' %rectangular shaped globe
                            m_grid('xtick',[],'ytick',[],'backcolor',[0 0 0]);
                        case 'flat'; %rectangular shaped globe    
                            m_grid('box','fancy','tickdir','in','backcolor',[0 0 0]);
                    end


            end






        case {'L2 swath','timeseries3 lambert','Monthly_data','GCM','GCM_COSP','CLOUDSAT-PRECIP','L3 time segment','CALIPSO','AMSRE','POLDER_daymean','SMOS'}


            %                m_proj('lambert','lon',[-180 -100],'lat',[50 90]);

            %                m_proj('lambert','lon',[-160 -135],'lat',[70 74]);
            if iseaice==1
                iwest = find(abs(Plon2D_edges)>179);
                Plat2D_edges(iwest)=NaN;
                Plon2D_edges(iwest)=NaN;
                P(iwest)=NaN;
            end


            if irestrict_domain==1                                
                % m_proj('lambert','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);
                
                m_proj('miller','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);
                
                
                
                %               m_proj('albers','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);
            else
                %Lon coords seem to need to be in form 0 to 360 for Miller
                %plots at least. Otherwise plot is distorted.
                %Think it may depend on where 0 deg lon is placed - think
                %they need to run continuously
                nlat=size(Plon2D,1);
                if minALL(diff(Plon2D(floor(nlat/2),:))) < -100
                    i180 = find(Plon2D<0);
                    Plon2D(i180) = Plon2D(i180)+360;
                    [Plon2D_edges,Plat2D_edges] = get_edges_lat_lon(Plon2D,Plat2D);
                end


                switch proj_type
                    case 'polar'
                        lat_polar= -90; %Antarctic
%                        lat_polar= 90; %Arctic
                        m_proj('stereographic','lat',lat_polar,'lon',0,'rad',45);
                    case 'ortho'                        
                        m_proj('ortho','lat',38,'long',-30);
                        
                    otherwise
                        %                m_proj('lambert','lon',[-180 -100],'lat',[50 90]);
                        %                    m_proj('lambert','lon',[-180 180],'lat',[50 90]);

                        %                m_proj('lambert','lon',[-120 -20],'lat',[-80 -50]);
                        %                m_proj('lambert','lon',[-80
                        %                -65],'lat',[-72 -68]);

                        %                m_proj('lambert','lon',[minALL(Plon) maxALL(Plon)],'lat',[minALL(Plat) maxALL(Plat)]);

                        %  **** NOTE! Cannot do a lambert projection if the latitude crosses the
                        %  equator! Try a Miller projection. *******
                        
                        %For some reason, if you use 0 to 360 for defining
                        %the Miller projection region to prevent the stripe at 180E. But need -180 to 180
                        %for the actual plot to prevent plotting strangeness. 
                        Plon2D_lims = Plon2D;
                        Plat2D_lims = Plat2D;
                        i180_2 = find(Plon2D_lims>180);
                        Plon2D_lims(i180_2) = Plon2D_lims(i180_2) - 360;
                        [Plon2D_lims_edges,Plat2D_lims_edges] = get_edges_lat_lon(Plon2D_lims,Plat2D_lims);
                        
                        m_proj('miller','lon',[minALL(Plon2D_lims_edges) maxALL(Plon2D_lims_edges)],'lat',[minALL(Plat2D_lims_edges) maxALL(Plat2D_lims_edges)]);
                        %m_proj('miller','lon',[0 360],'lat',[minALL(Plat2D_lims_edges) maxALL(Plat2D_lims_edges)]);

                        %                m_proj('lambert','lon',[-58 -50],'lat',[-62 -50]);
                    
                end


            end

            if ifull_swath==0
                if i_dpcolor==1
                    m_dpcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);
                    hold on
                    m_dpcolor(Plon2D_edges-360,Plat2D_edges,P); shading flat; %colormap(map);
                    m_dpcolor(Plon2D_edges+360,Plat2D_edges,P); shading flat; %colormap(map);
                else
                    m_pcolor(Plon2D,Plat2D,P); shading flat; %colormap(map);
                    hold on
                    m_pcolor(Plon2D-360,Plat2D,P); shading flat; %colormap(map);
                    m_pcolor(Plon2D+360,Plat2D,P); shading flat; %colormap(map);
                end
                
            else
                m_pcolor(Plon2_L2,Plat2_L2,P); shading flat; %colormap(map);
            end
            
            hold on

            clims_get=get(gca,'clim');

            %draws some extra info in the form of contours


            if icontour==1
%                cont_dat=Npoints;
%                cont_dat=squeeze(Ndays2);
%                cont_ints=[15 30 60 100 200 300];
%                cont_ints=[2 5 10 20 99];
                if i_dpcolor==1
                    %N.B. - 'w' colour doesn't show up when saving for some reason?
                    %Use 'c' (cyan) instead
                    [cont_out,hcont]=m_dcontour(Plon2D_edges-360,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                    [cont_out,hcont]=m_dcontour(Plon2D_edges,Plat2D_edges,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                else
                    [cont_out,hcont]=m_contour(Plon2D-360,Plat2D,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);
                    [cont_out,hcont]=m_contour(Plon2D,Plat2D,cont_dat,cont_ints,cont_col_str); shading flat; %colormap(map);

                end
                set(gca,'clim',clims_get);
                set(hcont,'linewidth',3);
                if contour_label==1
                    clabel(cont_out,hcont,'fontweight','bold');
                end
            end


            %                m_pcolor(Plon,Plat,P); shading flat; %colormap(map);
            

            
            %m_grid('box','fancy','tickdir','in','backcolor',[0 0 0],'xtick',lon_ticks,'ytick',lat_ticks);
            if iplot_mgrid_lines==1     
                
                if ioverride_ticks==1
                    % Plots lines, but with not labels
                    m_grid('xtick',[-180 -120 -60 0 60 120],'ytick',[-80:20:80],'backcolor',[0 0 0],'xticklabels',[],'yticklabels',[],'linestyle','--','color','k','linewidth',1);
                else
                    m_grid('box','fancy','tickdir','in','backcolor',[0 0 0]);
                end
                
                
               
            else
                %m_grid('box','fancy','tickdir','in','backcolor',[0 0 0],'linestyle','none');
                m_grid('xtick',[],'ytick',[],'backcolor',[0 0 0]); %this actually does lines, but with no labels
                
            end

            overpass_str='';
            if exist('file_name_h5') & strcmp(mod_data_type,'GCM')==0 & strcmp(mod_data_type,'GCM_COSP')==0 & strcmp(mod_data_type,'CLOUDSAT-PRECIP')==0 & strcmp(mod_data_type,'CALIPSO')==0 & strcmp(mod_data_type,'Sea-ice')==0 & strcmp(mod_data_type,'AMSRE')==0
                istr=strfind(file_name_h5,'MYD06');
                if length(istr)>0
                    overpass_str=['AQUA ' file_name_h5(istr+9:istr+21) ' '];
                end

                istr2=strfind(file_name_h5,'MOD06');
                if length(istr2)>0
                    istr=istr2;
                    overpass_str=['TERRA ' file_name_h5(istr+9:istr+21) ' '];
                end
                
                date_str_num = file_name_h5(istr+14:istr+16);
                date_str=datestr(datenum('01-Jan-2006')+str2num(date_str_num)-1);
                overpass_str = [overpass_str ' ' date_str];
                

            end
            
            
            title_info = [title_info overpass_str];


        case 'L2 swath 5km'
            Plat2D_edges=lat;
            Plon2D_edges=lon;
            


            if irestrict_domain==1
%                m_proj('lambert','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);
                m_proj('miller','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);                
            else
                m_proj('lambert','lon',[-180 -100],'lat',[50 90]);
            end

            m_pcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);
            hold on
            %                m_pcolor(Plon,Plat,P); shading flat; %colormap(map);
            m_grid('box','fancy','tickdir','in','backcolor',[0 0 0]);


        case 'mock L3'

            if irestrict_domain==1
                m_proj('lambert','lon',[thresh_LON(1) thresh_LON(2)],'lat',[thresh_LAT(1) thresh_LAT(2)]);
            else
                m_proj('lambert','lon',[-180 -100],'lat',[50 90]);
            end

            m_pcolor(Plon2D_edges,Plat2D_edges,P); shading flat; %colormap(map);
            hold on
            %                m_pcolor(Plon,Plat,P); shading flat; %colormap(map);
            m_grid('box','fancy','tickdir','in','backcolor',[0 0 0]);

            overpass_str='??';
            istr=strfind(file_name_h5,'MYD06');
            if length(istr)>0
                overpass_str=['AQUA ' file_name_h5(istr+9:istr+21)];
            end

            istr=strfind(file_name_h5,'MOD06');
            if length(istr)>0
                overpass_str=['TERRA ' file_name_h5(istr+9:istr+21)];
            end

            title_info = [title_info overpass_str];
            
        otherwise
            
            error('*** WARNING - no data plotted. Add to this switch. ***');

    end






    if detailed_BAS_coast==1
        for i=1:length(SHP)
            m_line(SHP(i).lon,SHP(i).lat,'color',[0.2 0.2 0.2],'linewidth',2);
        end

    else


        %m_coast('patch',[0.6 1 0.6]);
        hcoast=m_coast('line','linewidth',2,'color',[0.2 0.2 0.2]);
        %set(hcoast,'linewidth',2,'color','k');
        %m_grid('xaxis','middle','xtick',[]);

    end


    %        set(gcf,'Color',[0 0 0]); %makes the background black (so NaNs will appear black)
    
if supress_colorbar==0
    
    hc=colorbar('h');
    set(hc,'fontsize',fsize+2);
    %set(get(h,'title'),'string','MODIS plot'); %this puts a title on the
    %colorbar, so could be useful
    if ilabel_colorbar==1
        xlabel(hc,col_bar_lab_str,'fontsize',fsize+2);
    end
    %        set(gca,'Color',[0 0 0]);

    %need to set the colormap before chnage the cticks, or it goes
    %wrong
    colormap(colormap_choose);

    if inew_cticks==1
        set(gca,'clim',[min(x_cbar_vals) max(x_cbar_vals)]); %make sure we cover the full range
        set(hc,'xtick',x_cbar_vals); %set the tick marks to those requested
        ctick_text=ctick_strings_make(cbar_vals_for_text);
        set(hc,'xticklabel',ctick_text);
    end
    
    
end

    clims_get=get(gca,'clim');

    if iset_min_clim==1 & iset_max_clim==1
        set(gca,'clim',[clim_min clim_max]);
    elseif iset_min_clim==1
        set(gca,'clim',[clim_min clims_get(2)]);
    elseif iset_max_clim==1
        set(gca,'clim',[clims_get(1) clim_max]);
    end


    if strmatch(units_str_plot,'none')
        units_str_plot='';
    else
        units_str_plot=['(' units_str_plot ')'];
    end

%     if strcmp(thresh_str,'xxx')~=1
%         title_info=[title_info 'for ' thresh_str];
%     end

    title_info=[title_info ' mean=' num2str(Pmean)]

    short_plot_name=[MODIS_varname2_plot units_str_plot ' ' title_info];
    short_plot_name2=[MODIS_varname2_plot ' ' title_info];
    short_plot_name=remove_character(short_plot_name,'_','-');
    if length(short_plot_name)<71  %Ensure that goes onto two lines for consistency (unless goes on 3 lines..!)
        short_plot_name(end+1:71)=' ';
        %Doesn't seem to work...
    end
    titlenamw=textwrap({short_plot_name},71);
    title(titlenamw,'fontsize',fsize+2);

    %set(gca,'fontsize',fsize);
    %end

    if icolormap_cf_grey==1
        hmain=gca;
        axes(hc); %swith to the colobar axis
        ylims_cb=get(gca,'ylim');
        cb_lab_posy = ylims_cb(1) - diff(ylims_cb)*1.6; % a position below the cbar
        text(cb_lab_pos_left,cb_lab_posy,cb_lab_str_left,'fontsize',12)
        text(cb_lab_pos_right,cb_lab_posy,cb_lab_str_right,'fontsize',12);
        axes(hmain); %switch back to main axis
    end

    savename=[savedir short_plot_name2]

    clear ioverride_plotglobal_thresh %make sure this has been cleared if not already

    if plot_mpace_flight_path==1
        m_plot(mpace_lon,mpace_lat,'w-');
        m_plot(mpace_lon_mapped,mpace_lat_mapped,'k-');
    end

    if inew_figure==1 & i_increase_font_size_map_figures_OFF==0
    %external script
    increase_font_size_map_figures
    end
    
clear ioverride_plotglobal_thresh ioverride_month_select_multi_years ioverride_plotglobal_loc
catch plot_global_ERROR
    %clear the flag in case of a runtime error
    clear ioverride_plotglobal_thresh ioverride_month_select_multi_years ioverride_plotglobal_loc
    rethrow(plot_global_ERROR); %"re-issue" the error (also creates the links to the error etc)
end


