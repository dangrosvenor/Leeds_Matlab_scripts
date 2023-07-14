% CF PDFs, where the CF is the fraction of points over a 0.25x0.25 degree
% area with an LWP bigger than a threshold. But first the UM is coarse
% grained to the resolution of GOES.

idat_driver=0;
clear fileUM xdat_import ydat_import
LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08]; %GOES regin for UM comparison xkqk 26thOct POC

time_select = datenum('26-Oct-2008 17:00'); %for UM
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file

% -- For option setting see inside the loops



%--- Load and process the data


                            


%% ------------------------------
% ------ GOES data --------
% ------------------------------
idat_driver=idat_driver+1;


time_select = datenum('26-Oct-2008 17:00');


        
%------- Calculate the data to plot
         %read in the GOES data for the specific time
         ioverride_goes = 1;
         goes_action = 'load a particular file';        
         read_GOES_vocals_netcdf_files
        
% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='GOES';
        gcm_str='GOES';
       
        month_amsre = goes_month;
        year_amsre = goes_year;

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
        
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';

        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=200;
        
        logflag=0;
        dlogflag=0;
        
        isave_plot=0;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
        y_axis_vals = 'GOES LWP';
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
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
        ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
        xdat_import(idat_driver) = xdat_norm;
        labs_import(idat_driver).l = 'GOES';
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;




%% ------------------------------
% ------ UM data --------
% ------------------------------
idat_driver=idat_driver+1;

dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';
%idat=1;
%fileUM{idat_driver} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat_driver).l = '(xkqkh) 100cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat_driver} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat_driver).l = '(xkqkj) 400cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat_driver} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat_driver).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat_driver).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat_driver} = 'xkmph_LWP_RWP_.pp.nc'; labs_import(idat_driver).l = '(xkmph)'; pole_lat=70; pole_lon=284; idat=idat+1;



for idat=1:length(fileUM)
    
    
    filename = [dirUM fileUM{idat}];
    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    


        

        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select;
        [nc,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
        %pdf2d will then use nc to get the data
        
        lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        
        icoarsen=1;
        if icoarsen==1
        
        %average to the coarser resolution of goes
        d=diff(gcm_Plat2D_GOES,[],1);
        dlat_goes = meanNoNan(meanNoNan(d,1),1);
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_goes/dlat_UM));
        
        d=diff(gcm_Plon2D_GOES,[],2);
        dlon_goes = meanNoNan(meanNoNan(d,1),1);
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_goes/dlon_UM));
        
        %Doesn't work unless M=N so set this for now
        M=N;  %N is usually bigger (=5)
        
        

        lwp_UM_n5 = reduce_matrix_subsample_mean(lwp,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        
% -----   Now calculate the CF over 0.25x0.25 degree areas  (0.25 degrees
% is the same aas AMSRE)

    %Work out the number of points to average over to get 0.25x0.25 degree
    %resolution
     %average to the coarser resolution of goes
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_goes/dlat_UM));
        
        
        
        
        
        
        
        
        
        
        
        
        else
            lwp_UM_n5 = lwp;
        end
        
        
% ----- Set various things

          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
%        titlenam_driver = ['LWP for ' time_str ' ' labs_import(idat).l];
%        units_str_plot = 'g m^{-2}';
         
%        mod_data_type='AMSRE';
        gcm_str_select='UM';
        gcm_str='UM';

       
        month_amsre = [1:length(time_out)];
        year_amsre = [1:length(time_out)];

        
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
         watervap_defaults
         pdf2D_defaults
         
        
        %--- set some options for these particular plot loops
%        set_screening = {'none'};
%        modis_data_plot = 'Map of 2D data from outside driver script';
        i577 = 'MODIS_plot_UW';

        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=200;
        
        logflag=0;
        dlogflag=0;
        
        isave_plot=0;
        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
        y_axis_vals = 'UM LWP';
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
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
        ydat_import(idat_driver) = ydat_norm; %Use the non-cumulative PDF data
        xdat_import(idat_driver) = xdat_norm;
%        labs_import(idat_driver).l = labs_import(1).l;
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;


        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'LWP PDFs';
xlab='Liquid Water Path (g m^{-2})';
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


%---  Main script to do plots and save
DRIVER_lineplot_watervap

        
        
        
        
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
   
     
end
%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
