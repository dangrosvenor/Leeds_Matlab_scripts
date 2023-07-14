% LWP map plots for GOES data
try
    
i_mask_low_LWP=1; %Make any values below thresh_LWP equal to NaN
thresh_LWP_mask = 20;    


%What to plot
    driver_var = 'LWP';
%    driver_var = 'Visible reflectance';
    driver_var = 'Nd';
%    driver_var = 'CTH';

savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/GOES/';

if ~exist('ioverride_goes_multi') | ioverride_goes_multi==0
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file
    
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131645.nc'; %GOES file 13th Nov 18:45 UTC (14:00 LST on 13th)    
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131715.nc'; %GOES file 13th Nov 18:45 UTC (14:00 LST on 13th)        
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811131845.nc'; %GOES file 13th Nov 18:45 UTC (14:00 LST on 13th)
    
%    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811121845.nc'; %GOES file 12th Nov 18:45 UTC (14:00 LST on 12th)    
%    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200811141845.nc'; %GOES file 14th Nov 18:45 UTC (14:00 LST on 12th)        
%
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810191515.nc'; %   Hook case from Rhea's paper 
%    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810191845.nc'; %   Hook case from Rhea's paper, but at low LWP time (2pm)     
    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810172015.nc'; % Hook case from Rhea's paper, Fig. 2a (pre-hook)

%    goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810311845.nc'; %
    
    isave_plot_driver=0;
    
%% set this        
    data_source = 'load_goes_file'; %other option is 'goes_multi_array' for data stored in goes_LWP_multi, etc.
%    data_source = 'goes_multi_array'; imulti_goes=75;
end

icoarsen=0; N=6; M=6; %whether to coarsen the resolution by averaging over N*M pixels. Original data is approx 4x4km

irestrict_domain_driver=1; %flag to restrict the domain

LAT_val = [-23.5 -16.44]; LON_val = [-85.93 -78.08]; %Smaller region to but out model edges - old domain
LAT_val = [-20.5 -16.5]; LON_val = [-79 -75]; %Zooming in to top left
LAT_val = [-21 -19]; LON_val = [-76.5 -74.5]; %Zooming in to top left
LAT_val = [-20.50 -18]; LON_val =[-78.93 -76.43]; %Zooming in for high LWP, low Nd regions
%LAT_val = [-25 -20]; LON_val =[-86 -81]; %Zooming in for more open-like region on 12th Nov

LAT_val = [-20.70 -17.28]; LON_val =[-76.93 -73.08]; %Top half of UM domain for 12th Nov
LAT_val = [-22.70 -17.28]; LON_val =[-78.93 -73.08]; %FULL UM domain for 12th Nov
%LAT_val = [-25.4100  -14.5700]; LON_val =[-81.8550  -70.1550]; % 2 x FULL UM domain for 12th Nov
%LAT_val = [-30.5 -6.44]; LON_val = [-95.93 -68.08]; %Wider view
LAT_val = [-30 -10]; LON_val = [-90 -70]; %Wider view2 as used in Fig.1

iplot_markers=1; %See the end for what locations are being plotted




% -- For option setting see inside the loops



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';
clear fileUM xdat_import ydat_import
idat=1;

fileUM{idat} = 'GOES'; labs_import(idat).l = 'GOES'; pole_lat=70; pole_lon=284; idat=idat+1;

for idat=1:length(fileUM)


   
    nt_driver=1;
    


    for it_driver=1:nt_driver
        switch data_source
            case 'load_goes_file'
                %read in the GOES data for the specific time
                ioverride_goes = 1;
                goes_action = 'load a particular file';
                read_GOES_vocals_netcdf_files

                time_driver = gcm_time_matlab_GOES; %from read goes routine

            case 'goes_multi_array'
                time_driver = times_GOES_save(imulti_goes);
                goes_LWP = goes_LWP_multi{imulti_goes};
                goes_Nd = goes_Nd_multi{imulti_goes};
                goes_Vis = goes_Vis_multi{imulti_goes};                

        end
         
        %--- run the file to set up the defaults
        plot_global_maps_defaults  
        
        %--- Which data to plot
 
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

        iset_min_clim=0;
        clim_min=0;
        iset_max_clim=0;
        clim_max=300;
        

        

        savedir = savedir_driver;
        
        %Calculate the data to plot
        switch driver_var
            case 'LWP'
                dat_modis = goes_LWP; 
                
                if i_mask_low_LWP==1
                    dat_modis(dat_modis<thresh_LWP_mask)=NaN;                
                end
                
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=150;
                
            case 'Visible reflectance'
                dat_modis = goes_Vis;
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=1;
                
             case 'Nd'
                goes_Nd = MODIS_justN_func(goes_Tau,goes_Reff*1e-6,'calc',0,goes_Teff,'N');
                dat_modis = goes_Nd;
                
                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=350;
                
            case 'CTH'
%                dat_modis = goes_CTH/1e3;   %goes_Teff;  %originally in metres
                dat_modis = goes_Ztop;     %km - Minnis product?              
                iset_min_clim=1;
                clim_min=1.2;
                iset_max_clim=1;
                clim_max=2;
                
        end
        
        dat_modis(dat_modis<-8.99 & dat_modis>-9.01)=NaN;
        
        
         
         if icoarsen==1
             
             % For use later to coarsen the UM data
%              d=diff(gcm_Plat2D_GOES,[],1);
%              dlat_GOES = meanNoNan(meanNoNan(d,1),1);
%              d=diff(gcm_Plon2D_GOES,[],2);
%              dlon_GOES = meanNoNan(meanNoNan(d,1),1);
% 
%              dlat_target = dlat_GOES;
%              dlon_target = dlon_GOES;
% 
%              N = ceil(abs(dlat_target/dlat_GOES));
%              M = ceil(abs(dlon_target/dlon_GOES));

            %N=5; M=5;  %no. grid boxes to average together. GOES is at approx 4km resolution.


        %set these to zero since want to look at the effects of non-cloudy open
            %cell regions, etc.
%             dat_modis(isnan(dat_modis))=0; 
             
             dat_modis = reduce_matrix_subsample_mean(dat_modis,N,M);            
             gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
             gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
             %Work out the cell edges (as halfway between the centres)
             [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);

%          else
%              lwp_UM_n5 = lwp;
         end
        
        
        
        
        %Set various things
        time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
        time = time_driver(it_driver);
          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
        time_LST = time_driver(it_driver) + time_shift;
        time_str_LST = datestr(round(time_LST*24*60)/24/60,'HH:MM');         
        titlenam_driver = [driver_var ' for ' time_str ' UTC (' time_str_LST ' LST) ' labs_import(idat).l ' '];
        units_str_plot = '';
         
        mod_data_type='AMSRE';
        gcm_str_select='GOES';
        month_amsre = [1:length(time_driver)];
        year_amsre = [1:length(time_driver)];

        
%        i_dpcolor=1;
        ifull_swath=0;
        igcm_screen=0;
        
        

        
        %--- Apply override flags
        ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
        % iocean_only=1;
        ioverride_time_selection=0; %Override the times to include
        ioverride_plotglobal_loc=1; %Override the location of the plot window
        ioverride_years_time_screen=0; %Override years for screening?
        
        irestrict_domain = irestrict_domain_driver;
        thresh_LAT = LAT_val;
        thresh_LON = LON_val;
        
        %---  Run plot script and save
        plot_global_maps
%        colormap('gray');
        
        %Add a box for the UM simulation region
        plot_box_on_map
        
        increase_font_size_map_figures
        
        %Save the figure
        if isave_plot_driver==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
    end
   
     
end



if iplot_markers==1
        
        lat_mark01 = -20; lon_mark01 = -75;  %The ship
%        [ilat,ilon] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,lat_mark01,lon_mark01,0.1);
        m_plot(lon_mark01,lat_mark01,'wo','markersize',10,'markerfacecolor','k');
        
%        lat_mark02 = -20.8543; lon_mark02 = -76.4911; %The location where the UM matches pretty well
%        m_plot(lon_mark02,lat_mark02,'w^','markersize',15,'markerfacecolor','k');        
        
        
        
end
    

%  plot_box_on_map  - use this to plot the UM domain (set lat and lon in
%  there)

clear ioverride_goes_multi
catch goes_error
    clear ioverride_goes_multi
    rethrow(goes_error)
end





                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
