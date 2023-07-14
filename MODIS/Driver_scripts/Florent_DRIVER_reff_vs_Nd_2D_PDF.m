%Stripped down template for running plotTime* to produce a joint 2D PDF
%Note - have not tested this - might have missed some things out, so
%compare with POC_26Oct2008_CF_0pt25deg_PDFs_20141125T032943.m , which was
%stripped down to make this template.

 LAT_val = [-1e9 1e9]; LON_val = [-1e9 1e9];
 

gcm_str='UM';

isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
        
idat_driver=0;
clear fileUM xdat_import ydat_import


Ybins_driver = [0:5:350];       
Xbins_driver = [5:1:20];     
        

pdf_type_driver='normal';
pdf_type_driver='cumulative';

logbin_norm_driver = 0;
i_plot_norm_driver=1; %Whether to normalise
i_div_bin_widths_driver=1;  %whether to divide by the bin widths


        
% -- For option setting see inside the loops



%--- Load and process the data


  
    

%% ------------------------------
% ------ some data --------
% ------------------------------
idat_driver=idat_driver+1;


        


        
        
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
%        x_axis_vals = 'Dummy data'; %dummy data
        x_axis_vals = 'General GCM-style x-axis simple'
        
%        y_axis_vals = 'GOES LWP';
        y_axis_vals = 'General GCM-style';
        
        
%        ylabelstr = ['0.25^{o} cloud fraction for LWP.GT.' num2str(thresh_LWP_driver) ' g m^{-2}'];
        ylabelstr = ['Nd (cm^{-3})'];

        xlabelstr = ['Reff (\mum)'];
        
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
%       Put the data into the variables requried by the requested action in
%       pdf2d* . Plus define the gcm_Plon2D ,etc.

        X_driver = re_UM;
        
%E.g...
        Y_driver = N_UM;
%        gcm_Plat2D_GOES = reduce_matrix_subsample_mean(gcm_Plat2D_GOES,N,M);
%        gcm_Plon2D_GOES = reduce_matrix_subsample_mean(gcm_Plon2D_GOES,N,M);
        %Work out the cell edges (as halfway between the centres)
%        [gcm_Plat2D_edges_GOES, gcm_Plon2D_edges_GOES]=get_edges_lat_lon(gcm_Plat2D_GOES,gcm_Plon2D_GOES);


        

        
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
%        close(gcf);
%        waterVapourMay2005
%        close(gcf);
        