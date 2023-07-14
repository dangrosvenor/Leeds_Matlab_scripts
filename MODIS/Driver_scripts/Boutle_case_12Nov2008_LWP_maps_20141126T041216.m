% LWP map plots for all times for UM
% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov

irestrict_domain_DRIVER=1;

tol_mins=1;
time_single = datenum('13-Nov-2008 19:00'); time_range(1) = time_single - tol_mins/24/60; time_range(2) = time_single + tol_mins/24/60;


i_mask_low_LWP=1; %Make any values below thresh_LWP equal to NaN
thresh_LWP_mask = 20;

% -- For option setting see inside the loops



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
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

fileUM{idat} = '/xmmz-u/xmmzu_LWP_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar(increased)';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
fileUM{idat} = '/xmmz-v/xmmzv_LWP_.pp.nc.mat'; labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; pole_lat=70; pole_lon=284;
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;

        

for idat_UM=1:length(fileUM)
%     filename = [dirUM fileUM{idat}];
%     nc = netcdf(filename);
%     
%     time=nc{'t'}(:);
%     nt_driver=length(time);
%     t0_str=nc{'t'}.time_origin{1};
%     t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%     time_driver = datenum(t0_str2) + time;
%     
%     lon_UM = nc{'x'}(:);
%     lat_UM = nc{'y'}(:);
%     [lon2d,lat2d]=meshgrid(lon_UM,lat_UM);
%     %Convert to normal lat lon from rotated pole coords - make sure the
%     %rotated pole lat and lon are given correctly above
%     [lats2D_driver,lons2D_driver]=em2gm(lat2d,lon2d,pole_lat,pole_lon);  %From Annette M - see email for an example of how she uses it

    
    
    %------- Calculate the data to plot
         %read in the UM data for the specific time
%        time = time_select;
        
%         [nc,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = read_UM_file(filename,time,pole_lat,pole_lon);
%         %pdf2d will then use nc to get the data
%         
%         lwp = 1e3*nc{'LWP'}(it_driver,:,:); %convert to g/m2]
        
         vars_in.var = 'LWP';
         vars_in.flag = flag{idat_UM};
         vars_in.file_lwp =  [dirUM fileUM{idat_UM}]; 
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = []; %set to this for all times
    
     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    if length(time_range)==1
        time_range(2) = time_range(1);
    end
    it_driver_find = find(time_matlab>=time_range(1) & time_matlab<=time_range(2));
     
%% select time index here    
    for it_driver=it_driver_find %22:22 %1:length(time_matlab)  %4:nt_driver
        %12=23:00 UTC on 12th, 16 = 07:00 on 13th;  22=19:00 on 13th; 
        %--- run the file to set up the defaults
        plot_global_maps_defaults  
        
        irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not
        
        thresh_LAT = LAT_val_DRIVER;
        thresh_LON = LON_val_DRIVER;
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';

        iset_min_clim=1;
        clim_min=0;
        iset_max_clim=1;
        clim_max=300;
        
        isave_plot=0;
        iplot_markers=0;
        
        %Calculate the data to plot
        dat_modis = squeeze(lwp(it_driver,:,:)); %already in g/m2
        if i_mask_low_LWP==1
           dat_modis(dat_modis<thresh_LWP_mask)=NaN; 
        end
        
        %Set various things
        time = time_matlab(it_driver);
          %Round to the nearest minute as sometimes get 18:59:59
        time_str = datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM'); 
        titlenam_driver = ['LWP for ' time_str ' for ' labs_UM(idat_UM).l ' ' ];
        units_str_plot = 'g m^{-2}';
         
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
        if isave_plot==1
            saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
            close(gcf);
        end
%         
%         if imake_anim==1
%             
%             if ==1
%                 mov = avifile('/home/disk/eos1/d.grosvenor/modis_work/plots/MOVIE_LWP_xlhgw.avi','fps',1);
%             end
%             
%             
%             animfr(jj)=getframe(hh(jj).h);
%             mov = addframe(mov,animfr(jj));
%             
%             
%         end
        
    end
    
    if iplot_markers==1
        
        lat_mark01 = -20; lon_mark01 = -75;  %The ship
%        [ilat,ilon] = getind_latlon_quick(gcm_Plat2D_UM,gcm_Plon2D_UM,lat_mark01,lon_mark01,0.1);
        m_plot(lon_mark01,lat_mark01,'wo','markersize',15,'markerfacecolor','k');
        
        lat_mark02 = -20.8543; lon_mark02 = -76.4911; %The location where the UM matches pretty well
        m_plot(lon_mark02,lat_mark02,'w^','markersize',15,'markerfacecolor','k');        
        
        
        
    end
   
     
end
%    xdat_import(idat).x =



