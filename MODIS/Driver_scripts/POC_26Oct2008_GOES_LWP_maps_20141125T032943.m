% LWP map plots for GOES data

goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file


% -- For option setting see inside the loops



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';
clear fileUM xdat_import ydat_import
idat=1;

fileUM{idat} = 'GOES'; labs_import(idat).l = 'GOES'; pole_lat=70; pole_lon=284; idat=idat+1;

for idat=1:length(fileUM)


    time_driver = gcm_time_matlab_GOES(1);
    nt_driver=1;
    


    for it_driver=1:nt_driver
        
        %------- Calculate the data to plot
         %read in the GOES data for the specific time
         ioverride_goes = 1;
         goes_action = 'load a particular file';        
         read_GOES_vocals_netcdf_files
        
        %--- run the file to set up the defaults
        plot_global_maps_defaults   
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=200;
        
        isave_plot=1;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/GOES/';
        
        %Calculate the data to plot
        dat_modis = goes_LWP; 
%        dat_modis = goes_Vis; 
        
        %Set various things
        time = time_driver(it_driver);
          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
        titlenam_driver = ['LWP for ' time_str ' ' labs_import(idat).l ' '];
        units_str_plot = 'g m^{-2}';
         
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
        
        %---  Run plot script and save
        plot_global_maps
        
        %Add a box for the UM simulation region
        plot_box_on_map
        
        %Save the figure
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
    end
   
     
end
%    xdat_import(idat).x =





                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
