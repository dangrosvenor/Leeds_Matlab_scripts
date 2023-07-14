% CF PDFs, where the CF is the fraction of points over a 0.25x0.25 degree
% area with an LWP bigger than a threshold. But first the UM is coarse
% grained to the resolution of GOES.
isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
        
idat_driver=0;
clear fileUM xdat_import ydat_import
LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08]; %GOES regin for UM comparison xkqk 26thOct POC
%LAT_val = [-24.5 -15.44]; LON_val = [-84.5 -77.08]; %Trying to match AMSRE and GOES domains
%LAT_val_GOES = LAT_val; LON_val_GOES = LON_val;
LAT_val = [-23.5 -16.44]; LON_val = [-85.93 -78.08]; %Smaller region to but out model edges
LAT_val = [-21 -16.44]; LON_val = [-82 -78.08]; %Much smaller region in NW corner

LAT_val = [-21 -15.44]; LON_val = [-86.93 -82]; %Revised smaller region in NW corner  

LAT_val_GOES = LAT_val; LON_val_GOES = LON_val;  %Choice of GOES lat lon
%seems not to make much difference
LAT_val_UM = LAT_val; LON_val_UM = LON_val;

%Here can choose different LAT and LON zones for the model and satellite



time_select = datenum('26-Oct-2008 17:00'); %for UM
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261645.nc'; %GOES file
goes_file_to_load = 'GOES10_cld_ret_VOCALS_200810261745.nc'; %GOES file

Ybins_driver = [0:0.1:1];
Ybins_driver = [0:0.05:1];        
        
%threshold to define cloud:-
thresh_LWP_driver = 10; %g/m2  -- ideally would do this with optical depth

target_dlat = 0.25;
target_dlon = 0.25;

%target_dlat = 0.5;
%target_dlon = 0.5;

pdf_type_driver='normal';
pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths


        
% -- For option setting see inside the loops



%--- Load and process the data


 %Quick fix to get the linestyles to be the same as without AMSRE
    idat_driver=idat_driver+1;
    labs_import(idat_driver).l = ' ';
    xdat_import(idat_driver).x=NaN;
    ydat_import(idat_driver).y=NaN;                             


%% ------------------------------
% ------ GOES data --------
% ------------------------------
idat_driver=idat_driver+1;

LAT_val = LAT_val_GOES;
LON_val = LON_val_GOES;

        
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
         pdf2D_defaults  %for pdf2D_plot_commands
         
        
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
        

        
        

                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
%        y_axis_vals = 'GOES LWP';
        y_axis_vals = 'General GCM-style';
        datatype = 'gcm_data';
        
        
        ylabelstr = ['0.25^{o} cloud fraction for LWP.GT.' num2str(thresh_LWP_driver) ' g m^{-2}'];
        Ybins = Ybins_driver; ichoose_Ybins=1;
        
        graph = 977; %new 1D PDF from 2D histo data - can choose either axis
                                %(for watervap)
                                
                axis1D = 'y';                                
                                
                logbin_norm = logbin_norm_driver;
                i_plot_norm=i_plot_norm_driver;
                i_div_bin_widths=i_div_bin_widths_driver;
                pdf_type = pdf_type_driver;
                
%        gcm_str = gcm_str_last_loaded;      


% ------- Calculate the data --------
       lwp_GOES = 5/9*1e3*1e-6*goes_Reff.*goes_Tau *1e3; %g/m2

      %Work out the number of points to average over to get target_dlat x target_dlon degree
      %resolution

        d=diff(gcm_Plat2D_GOES,[],1);
        dlat_GOES = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(target_dlat/dlat_GOES));
        

        d=diff(gcm_Plon2D_GOES,[],2);
        dlon_GOES = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(target_dlon/dlon_GOES));     
        
%        N=5; M=12;
%        N=10; M=10;
        
     %Cloud fraction will be defined as fraction of points within each N*M
     %box with an LWP greater than a threshold. Total no. points =N*M
        %Make an array of ones and make the ones that we don't want to
        %count zero
        Nlwp = zeros(size(lwp_GOES));
        Nlwp(lwp_GOES>=thresh_LWP_driver)=1;
        %Now coarse grain (avearge) the Nlwp array - this will now be our
        %cloud fraction
        cf_GOES_0pt25 = reduce_matrix_subsample_mean(Nlwp,N,M);
        Y_driver = cf_GOES_0pt25;
        gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
        gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);
        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3  %Uses pdf2D_plot_commands
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
                ylab = 'Cumulative Frequency';

        end
        
        labs_import(idat_driver).l = 'GOES';
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;




%% ------------------------------
% ------ UM data --------
% ------------------------------
LAT_val = LAT_val_UM;
LON_val = LON_val_UM;

for idat=1:99
    flag{idat}='';  %The default is netCDF
end

dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
%dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';


idat=1;
%fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = 'Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_UM(idat).l = '400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
% fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_UM(idat).l = '400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkq_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = 'No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1; %1000cm^{-3} 

%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkx_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;

%fileUM{idat} = 'xkqkm_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkn_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;

fileUM{idat} = 'xkqkv_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkw_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;

for idat_UM=1:length(fileUM)

    idat_driver=idat_driver+1;
    
    filename = [dirUM fileUM{idat_UM}];
    
    %Read in all the times in case we want to use them all
%    [nc,time_driver,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,[],pole_lat,pole_lon);
   

    


        

        
%------- Calculate the data to plot
         %read in the UM data for the specific time
        time = time_select;
        %        [nc,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);

%        lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
                
         vars_in.var = 'LWP';
         vars_in.flag = flag{idat_UM};
         vars_in.file_lwp = filename;
         vars_in.file_rho = ''; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = time_select;
    
     [lwp,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
     
        lwp = lwp*1e3;
        
        
        icoarsen=1;
        if icoarsen==1
        
        %average to the coarser resolution of goes
        %dlat_GOES and dlon_GOES calculated earlier for the original grid -
        %don't want to do again here since the Plat2D values have been
        %changed.
        
%         d=diff(gcm_Plat2D_GOES,[],1);
%         dlat_GOES = meanNoNan(meanNoNan(d,1),1);
        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(dlat_GOES/dlat_UM));
        
%         d=diff(gcm_Plon2D_GOES,[],2);
%         dlon_GOES = meanNoNan(meanNoNan(d,1),1);
        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(dlon_GOES/dlon_UM));        
        
        

        lwp_UM_n5 = reduce_matrix_subsample_mean(lwp,N,M);
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        
% -----   Now calculate the CF over 0.25x0.25 degree areas  (0.25 degrees
% is the same aas AMSRE)

    %Work out the number of points to average over to get 0.25x0.25 degree
    %resolution
            

        d=diff(gcm_Plat2D_UM,[],1);
        dlat_UM = meanNoNan(meanNoNan(d,1),1);
        N = ceil(abs(target_dlat/dlat_UM));
        

        d=diff(gcm_Plon2D_UM,[],2);
        dlon_UM = meanNoNan(meanNoNan(d,1),1);
        M = ceil(abs(target_dlon/dlon_UM));     
        

        
     %Cloud fraction will be defined as fraction of points within each N*M
     %box with an LWP greater than a threshold. Total no. points =N*M
        %Make an array of ones and make the ones that we don't want to
        %count zero
        Nlwp = zeros(size(lwp_UM_n5));
        Nlwp(lwp_UM_n5>thresh_LWP_driver)=1;
        %Now coarse grain (avearge) the Nlwp array - this will now be our
        %cloud fraction
        cf_UM_0pt25 = reduce_matrix_subsample_mean(Nlwp,N,M);
        Y_driver = cf_UM_0pt25;
        gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
        gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
        %Work out the cell edges (as halfway between the centres)
        [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
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
        
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
%        y_axis_vals = 'UM LWP';
        y_axis_vals = 'General GCM-style';
        
        ylabelstr = ['0.25^{o} cloud fraction for LWP.GT.' num2str(thresh_LWP_driver) ' g m^{-2}'];
        Ybins = Ybins_driver; ichoose_Ybins=1;
        
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
                ylab='Cumulative Frequency';
        end
        
        labs_import(idat_driver).l = labs_UM(idat_UM).l;
        xlab_import = xlab;
        ylab_import = ylab;
%        ioverride_savePDF=1;
%        save_1D_pdfs;

end


        
%% ------------------------------
% ------ plot the combined PDF using case 0 of watervap --------
% ------------------------------


%--- run the file to set up the defaults
watervap_defaults

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'CF PDFs';
%xlab='Liquid Water Path (g m^{-2})';
xlab = xlab_import;
ylab = ylab_import;
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=isave_plot_driver;
savedir = savedir_driver;

%idate_ticks_fix=1;
%iaxis_square=0; %switch to make axis square


%---  Main script to do plots and save
DRIVER_lineplot_watervap

        
        
        
        
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
   
     

%    xdat_import(idat).x =







                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
