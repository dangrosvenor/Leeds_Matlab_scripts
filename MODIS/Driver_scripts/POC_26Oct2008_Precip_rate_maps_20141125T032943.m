% LWP map plots for all times for UM
isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

idat_driver=0;

clim_driver = [0 0.4]*24; %Default for climits

icoarse = 0; %flag for whether to degrade the UM resolution to that of GOES

time_select = 'ALL';
time_select = datenum('26-Oct-2008 17:00');

% -- For option setting also see inside the loops

for idat=1:99
    flag{idat} = 'load_UM';
end


%--- Load and process the data
dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';
%dirUM='/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/';
clear fileUM xdat_import ydat_import


idat=1;
%%fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = '(xkqkf) Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
%%fileUM{idat} = 'xkqkh'; labs_UM(idat).l = '(xkqkh) 100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkj'; labs_UM(idat).l = '(xkqkj) 400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
% fileUM{idat} = 'xkqkk'; labs_UM(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqko'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkq'; labs_UM(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
%fileUM{idat} = 'xkqkr'; labs_UM(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1;
%%fileUM{idat} = 'xkqkm'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkn'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1; clim_driver = [0 100];

clear vars_UM vars_in

idat=1;
vars_UM{idat}='LWP_RWP'; idat=idat+1;
%vars_UM{idat}='Nd'; idat=idat+1;
vars_UM{idat}='rho'; idat=idat+1;
vars_UM{idat}='qR'; idat=idat+1;
vars_UM{idat}='NR'; idat=idat+1;
vars_UM{idat}='rain_rate'; idat=idat+1;

for idat_UM=1:length(fileUM)
    idat_driver = idat_driver+1;
    
%    iund = findstr(fileUM{idat_UM},'_');
%    filename_start = fileUM{idat_UM}(1:iund(1)-1);
     filename_start = fileUM{idat_UM};   

    for ivar=1:length(vars_UM)
        filename = [dirUM filename_start '_' vars_UM{ivar} '_.pp.nc'];
        switch flag{idat_UM}
            case 'load_mat'
                filename = [filename '.mat'];                
        end
        
        eval(['vars_in.file_' vars_UM{ivar} ' = filename;']);
        
    end
    

%     vars_in.var = 'Nd';
%     vars_in.flag = flag{idat_UM};
% %    vars_in.file_lwp = filename;
% %    vars_in.file_rho = filename_rho;
%     vars_in.pole_lat = pole_lat;
%     vars_in.pole_lon = pole_lon;
%     vars_in.time_in = [];

    vars_in.var = 'LWP';  %set as LWP just to get the times for the file (quicker than Nd)
    vars_in.flag = flag{idat_UM};
%    vars_in.file_lwp = file_LWP_RWP;    
%    vars_in.file_rho = ''; %filename_rho;
    vars_in.pole_lat = pole_lat;
    vars_in.pole_lon = pole_lon;
    vars_in.time_in = [];

    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
   
    if isstr(time_select)==1 & strcmp(time_select,'ALL')==1
       time_select = time_matlab; 
    end
    
%    nt_driver=length(time_driver);        
    for it_driver=1:length(time_select)   %1:nt_driver
        %--- run the file to set up the defaults
        plot_global_maps_defaults   
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

        iset_min_clim=1;
        clim_min=clim_driver(1);
        iset_max_clim=1;
        clim_max=clim_driver(2);
        

        
        %Calculate the data to plot
         %repeat the read-in to get the specific time
         time = time_select(it_driver);
%        [nc,time_out,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);

        vars_in.var = 'Precip_rate';
        vars_in.time_in = time_select(it_driver);
        [precip_rate,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);                
        
        dat_modis = 24*3600*precip_rate; %convert to mm/day
        
        if icoarse==1
            %Coarsen the UM data to approximate the resolutin of GOES
            N=5;
            M=3;
            %See POC_26Oct2008_CF_0pt25deg_PDFs_20141125T032943.m for
            %derivation of M and N
            
            dat_modis = reduce_matrix_subsample_mean(dat_modis,N,M);
            gcm_Plat2D_UM = reduce_matrix_subsample_mean(gcm_Plat2D_UM,N,M);
            gcm_Plon2D_UM = reduce_matrix_subsample_mean(gcm_Plon2D_UM,N,M);
            %Work out the cell edges (as halfway between the centres)
            [gcm_Plat2D_edges_UM, gcm_Plon2D_edges_UM]=get_edges_lat_lon(gcm_Plat2D_UM,gcm_Plon2D_UM);
            
            
            
        end

        
        %Set various things

          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
        titlenam_driver = ['Column max rain rate for ' time_str ' ' labs_UM(idat_UM).l];
        units_str_plot = 'mm day^{-1}';
         
        mod_data_type='AMSRE';
        gcm_str_select='UM';

       
        month_amsre = [1:length(time_matlab)];
        year_amsre = [1:length(time_matlab)];

        
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
        savedir = savedir_driver;
        plot_global_maps
        if isave_plot_driver==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
    end
   
     
end
%    xdat_import(idat).x =





                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
