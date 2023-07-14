% Map plots L2 CERES data
% Run read_CERES first, and choose the swaths required to join together.

icut_for_CERES=1; %Cut out code not needed for the CERES analysis

%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];
irestrict_domain_DRIVER=1; %set to zero for the full domain
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
%LAT_val_DRIVER = [-30 -5.0]; LON_val_DRIVER =[-90 -65]; %
%LAT_val_DRIVER = [-30 -5.0]; LON_val_DRIVER =[-180 -65]; % very large region
%LAT_val_DRIVER = [-25.70 -14.28]; LON_val_DRIVER =[-81.93 -70.08]; %around 2x region around UM domain for 12th Nov

% -- For option setting see inside the loops


goes_var = 'SW TOA up';
goes_var = 'LW TOA up';

time_shift = -(4+48/60) /24; time_shift_ref=-76; %amount to shift time by for LST (from UTC) for the UM domain (76W).
clear fileUM xdat_import ydat_import
idat=1;

if icut_for_CERES==1
    load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T064655.mat';
    load(load_file_goes);  %




%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkh) 100cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkj) 400cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;idat=idat+1;
%fileUM{idat} = 'xkmph_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkmph)'; pole_lat=70; pole_lon=284; idat=idat+1;

% - 12th Nov case, as of analysis started July 2015

%fileUM{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
fileUM{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        

end


        
for idat_UM=1:1   %length(fileUM)
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
        
%          vars_in.var = 'Nd';
%          vars_in.flag = flag{idat_UM};
%          vars_in.file_lwp =  [dirUM fileUM{idat_UM}]; 
% %         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
%          vars_in.pole_lat = pole_lat;
%          vars_in.pole_lon = pole_lon;
%          vars_in.time_in = []; %set to this for all times
%     
%      [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    
    
    for it_driver=1:1 %1:length(time_matlab)  %4:nt_driver
        %36=14:00  LST on 12th - time is stored in time = times_GOES_save(it_driver) + time_shift
        %79=14:00 LST on 13th
        %--- run the file to set up the defaults
        plot_global_maps_defaults  
        
        irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not
        
        thresh_LAT = LAT_val_DRIVER;
        thresh_LON = LON_val_DRIVER;
        
        %--- set some options for these particular plot loops
        set_screening = {'none'};
        modis_data_plot = 'Map of 2D data from outside driver script';
        
        switch goes_var
            case 'SW TOA up'

                iset_min_clim=1;
                clim_min=0;
                iset_max_clim=1;
                clim_max=900;

                isave_plot=0;
                iplot_markers=0;

                dat_modis = SW_CERES_combined;

                units_str_plot = 'W m^{-2}';
                var_str='SW TOA up';  
                
            case 'LW TOA up'

                iset_min_clim=1;
                clim_min=250;
                iset_max_clim=1;
                clim_max=350;

                isave_plot=0;
                iplot_markers=0;

                dat_modis = LW_CERES_combined;

                units_str_plot = 'W m^{-2}';
                var_str='LW TOA up';                 

           

        end
        %Set various things
        
        
        %Find the indices for all the values in the region of the plot
        itime = find( gcm_Plat2D_CERES >= LAT_val_DRIVER(1) & gcm_Plat2D_CERES <= LAT_val_DRIVER(2) & gcm_Plon2D_CERES >= LON_val_DRIVER(1) & gcm_Plon2D_CERES <= LON_val_DRIVER(2) );
        
        time_local = time_CERES_combined + (gcm_Plon2D_CERES-time_shift_ref)/15/24 + time_shift; %15 degrees is one hour of time difference.
        
%        time = times_GOES_save(it_driver) + time_shift;  %time_matlab(it_driver) + time_shift;        
        time = meanNoNan(time_local(itime),1); %Find an average time - or a range might be better
        
        if ~isnan(time)        
          %Round to the nearest minute as sometimes get 18:59:59
          time_str = [datestr(round(time*24*60)/24/60,'dd-mmm-yyyy HH:MM') ' LST '];
        else
           time_str='NO TIME DATA ';
        end
        titlenam_driver = [var_str ' for ' time_str];
    
         
        mod_data_type='AMSRE';
        month_amsre = [1:length(time)]; %need to set these
        year_amsre = [1:length(time)]; 
        
        gcm_str_select='CERES';
        
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
            saveas_ps_fig_emf(gcf,[savename],'',0,1);
            close(gcf);
        end
        
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

ioverride_box_colour=1;
col_str='k-';
plot_box_on_map
