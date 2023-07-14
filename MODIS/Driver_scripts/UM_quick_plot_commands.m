%This is usuall run from :-
%   UM_quick_plot_global.m

%isave_plot=1;

clear dat_modis_save time_save



LAT_val_DRIVER = [-1e9 1e9];
LON_val_DRIVER = [-1e9 1e9]; %Can restrict the plot region with this if set irestrict_domain_DRIVER=1

LAT_val_DRIVER = [-40 0];
LON_val_DRIVER = [-100 -60];

% LAT_val_DRIVER = [-24 -16];
% LON_val_DRIVER = [-80 -72];

LAT_val_DRIVER = [-30 -10]; LON_val_DRIVER = [-86 -66]; %Close up of region surrounding the UM domain

LAT_val_DRIVER = [20 80]; LON_val_DRIVER = [-86 0]; %N Atlantic

%LAT_val_DRIVER = [-30 30]; LON_val_DRIVER = [-43 43]; %Africa

irestrict_domain_DRIVER=0;



switch file_type
    case {'nc','nc_emiss_CMIP'}
        nc = netCDF(filename);
end

gcm_str_DRIVER='UM';  %N.B. - have to set it to this and not another name - doesn't work if do (code would need changing)
%[gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(MLON,MLAT);

switch file_type
    case 'nc'


        %Slices have the lat lon converted from rotated coords for the nest
        %already.
        gcm_Plat2D_UM = nc{lat_var}(:);
        gcm_Plon2D_UM = nc{lon_var}(:);



        i180 = find(gcm_Plon2D_UM>180);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

        time_glm = nc{'Time'}(:);

    case 'nc_emiss_CMIP'

        lat = nc{'grid_latitude'}(:);
        lon = nc{'grid_longitude'}(:);
        [lon2d,lat2d] = meshgrid(lon,lat);

        [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(lat2d,lon2d,pole_lat,pole_lon);

        %        i180 = find(gcm_Plon2D_UM>180);
        %        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

        time_glm = nc{'time'}(:);

    case 'nc_emiss_CMIP_global'

        lat = nc{'latitude'}(:);
        lon = nc{'longitude'}(:);

        [gcm_Plon2D_UM,gcm_Plat2D_UM] = meshgrid(lon,lat);

        %        [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(lat2d,lon2d,pole_lat,pole_lon);

        i180 = find(gcm_Plon2D_UM>180);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

        time_glm = nc{'time'}(:);


    case 'mat'
        load(filename);
        [lon2d,lat2d] = meshgrid(lon,lat);

        [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(lat2d,lon2d,pole_lat,pole_lon);

        %         i180 = find(gcm_Plon2D_UM>180);
        %        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

        time_glm=(datenum('01-Aug-1970') - datenum('01-Jan-1970') )*24;
        it=1;

end

[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);

%Back up the original lat lons
gcm_Plat2D_UM_read = gcm_Plat2D_UM;
gcm_Plon2D_UM_read = gcm_Plon2D_UM;
gcm_Plat2D_edges_UM_read = gcm_Plat2D_edges_UM;
gcm_Plon2D_edges_UM_read = gcm_Plon2D_edges_UM;


time_gml_matlab = datenum('01-Jan-1970') + time_glm/24;

daynum_timeseries3_UM = 1; %[1:length(time)];
gcm_time_UTC_UM = 1; %[1:length(time)];

isquare=0;

if isstr(itimes)
    switch itimes
        case 'all'
            itimes=1:length(time_glm);
    end
end



if iUTC==1
    time_shift=0;
    time_format_str=' UTC';
else
    time_format_str=' LST';
end


for it2=1:length(itimes)
%     %Revert back to originals (since may be changed by icoarse)
%     gcm_Plat2D_UM_read = gcm_Plat2D_UM;
%     gcm_Plon2D_UM_read = gcm_Plon2D_UM;
%     gcm_Plat2D_edges_UM_read = gcm_Plat2D_edges_UM;
%     gcm_Plon2D_edges_UM_read = gcm_Plon2D_edges_UM;
    

    it=itimes(it2);

    time_LST = time_gml_matlab(it) + time_shift;

    [Y,M,D,HH,MM,SS]=datevec(time_LST);
    if SS>59
        SS=0;
        MM=MM+1;
    end
    time_round = datenum(Y,M,D,HH,MM,SS);

    clear clim_min clim_max


    units_str_plot = '';



    isplit=0;
    if isplit==1

        file_grid = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-at459/umglaa_LWP_ALL.mat';
        dat_grid = load(file_grid);

        lat_target = dat_grid.gcm_Plat2D_UM(:,1);
        lon_target = dat_grid.gcm_Plon2D_UM(1,:);

        ilat=find(lat_target>=LAT_val_DRIVER(1) & lat_target<=LAT_val_DRIVER(2));
        ilon=find(lon_target>=LON_val_DRIVER(1) & lon_target<=LON_val_DRIVER(2));

        %N216 is around 1 degree, so can fit around 25 4km boxes in 1x1deg
        lat_fine = split_matrix_int(dat_grid.gcm_Plat2D_UM(ilat,ilon),25);

        [lat_coarse] = reduce_matrix_subsample_mean(lat_fine,25,25);

    end



    %%
    switch var_UM
        case 'wind speed';

            U_glm = nc{'U'}(it,:,:);
            V_glm = ncV{'V'}(it,:,:);
            V_glm = ncV{'V'}(it,1:end-1,:);

            dat_modis = sqrt(U_glm.^2 + V_glm.^2);

            titlenam_driver = ['Wind speed (m s^{-1}) at ' datestr(time_round) time_format_str];
            units_str_plot = '';

        case 'accum_number_ukca'
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['Soluble accumulation mode number (# per #molecules air) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 3e-17;

        case 'accum_number_ukca diff'
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['Difference'];

            clim_min = -0.5e-170;
            clim_max = 0.5e-17;

        case 'accum_mass_H2SO4_ukca'
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['Soluble accumulation mode H2SO4 MMR (kg kg^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-9;

        case 'accum_mass_OC_ukca'
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['Soluble accumulation mode OC MMR (kg kg^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 2e-9;

        case 'accum_mass_BC_ukca'

            dat_modis = nc{var_UM}(it,:,:);
            clim_min = 0;
            clim_max = 5e-10;

            titlenam_driver = ['Soluble accumulation mode BC MMR (kg kg^{-1}) at ' datestr(time_round) time_format_str];

        case 'emissions_BC_biomass'
            dat_modis = BC_model_level01_time08;
            titlenam_driver = ['BC biomass emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];
            clim_min = 0;
            clim_max = 5e-11;


        case 'emissions_BC_biofuel'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-12;

        case 'emissions_BC_fossil'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['BC fossil fuel emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-12;

        case 'emissions_DMS'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['DMS emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-12;

        case 'emissions_Monoterp'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['Monoterpene emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-12;

        case 'emissions_OC_biofuel'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 5e-12;

            %     case 'emissions_BC_biofuel'
            %
            %         dat_modis = nc{var_UM}(it,:,:);
            %         titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];
            %
            %         clim_min = 0;
            %         clim_max = 5e-12;
            %
            %     case 'emissions_BC_biofuel'
            %
            %         dat_modis = nc{var_UM}(it,:,:);
            %         titlenam_driver = ['BC biofuel emissions (kg m^{-2} s^{-1}) at ' datestr(time_round) time_format_str];
            %
            %         clim_min = 0;
            %         clim_max = 5e-12;


        case 'LWP'

            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['LWP (g m^{-2}) at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 300;
            
        case 'Nd_lwc_weighted_UKCA'
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = ['UKCA N_d weighted by LWC at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 300;            
            clim_max = 200;

        otherwise
            dat_modis = nc{var_UM}(it,:,:);
            titlenam_driver = [var_UM ' at ' datestr(time_round) time_format_str];

            clim_min = 0;
            clim_max = 300;



    end


    if icoarse==1
        switch coarse_grid
            case 'global grid'
                file_grid = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-at459/umglaa_LWP_ALL.mat';
                dat_grid = load(file_grid);


        
            case '1 deg'
                MLAT=[-89.5:89.5];
                MLON=[-179.5:179.5];
                
                [dat_grid.gcm_Plon2D_UM,dat_grid.gcm_Plat2D_UM] = meshgrid(MLON,MLAT);     
                
                [dat_grid.gcm_Plon2D_edges_UM,dat_grid.gcm_Plat2D_edges_UM] = get_edges_lat_lon(dat_grid.gcm_Plon2D_UM,dat_grid.gcm_Plat2D_UM);
                
                
        end
        
        lat_target = dat_grid.gcm_Plat2D_UM(:,1);
        lon_target = (dat_grid.gcm_Plon2D_UM(1,:))';

        dlat_target = mean(abs(diff(lat_target)));
        dlon_target = mean(abs(diff(lon_target)));

        %    ilat=find(lat_target>=LAT_val_DRIVER(1) & lat_target<=LAT_val_DRIVER(2));
        %    ilon=find(lon_target>=LON_val_DRIVER(1) & lon_target<=LON_val_DRIVER(2));

        %lat2d = dat_grid.gcm_Plat2D_UM(ilat,ilon);
        %lon2d = dat_grid.gcm_Plon2D_UM(ilat,ilon);

        %    coarsen_method = 'lat lon blocks';
        [dat_modis,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM] = coarse_grain(dat_modis,gcm_Plat2D_UM_read,gcm_Plon2D_UM_read,dlat_target,dlon_target);

        dat_modis = griddata(gcm_Plat2D_UM(:),gcm_Plon2D_UM(:),dat_modis(:),dat_grid.gcm_Plat2D_UM,dat_grid.gcm_Plon2D_UM);

        gcm_Plat2D_UM = dat_grid.gcm_Plat2D_UM;
        gcm_Plon2D_UM = dat_grid.gcm_Plon2D_UM;
        gcm_Plat2D_edges_UM = dat_grid.gcm_Plat2D_edges_UM;
        gcm_Plon2D_edges_UM = dat_grid.gcm_Plon2D_edges_UM;

    end



    %%


    %--- run the file to set up the defaults
    plot_global_maps_defaults

    irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not

    thresh_LAT = LAT_val_DRIVER;
    thresh_LON = LON_val_DRIVER;

    %--- set some options for these particular plot loops
    set_screening = {'none'};
    modis_data_plot = 'Map of 2D data from outside driver script';

    if exist('clim_min')
        iset_min_clim=1;
        %clim_min=0;
        iset_max_clim=1;
        %clim_max=64800;
    end


    iplot_markers=0;



    mod_data_type='AMSRE';
    gcm_str_select = gcm_str_DRIVER;
    %        daynum_timeseries3_UM = [1:length(time)];
    %        gcm_time_UTC_UM = [1:length(time)];


    ifull_swath=0;
    igcm_screen=0;




    %--- Apply override flags
    ioverride_plotglobal_thresh=1; %Override most of the options (what to plot, etc.)
    % iocean_only=1;
    ioverride_time_selection=0; %Override the times to include
    ioverride_plotglobal_loc=1; %Override the location of the plot window
    ioverride_years_time_screen=0; %Override years for screening?
    iover_ride_plot_global=1; %overrides inew_figure=1; supress_colorbar=0; i_increase_font_size_map_figures_OFF = 0;
    %(all set in plot_global_maps_defaults)

    %---  Run plot script and save
    plot_global_maps
    %-------------------------------

    title(titlenam_driver);


    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/ACSIS_plots_general/' num2str(ifile,'%03d') '_' titlenam_driver];

    if isave_plot==1
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
%        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0);
    end
    if iclose_plot==1
        close(gcf);
    end



    % General instructions for case 'Generic plot specified outside of script'
    %Run plot_global_maps_defaults (might want to change these)
    %Then need to set e.g. P, gcm_Plat2D_CERES, gcm_Plon2D_CERES,
    % gcm_Plat2D_edges_CERES, gcm_Plon2D_edges_CERES],
    % gcm_str, daynum_timeseries3_CERES = 1; gcm_time_UTC_CERES = 1; month_amsre=1; year_amsre=1;
    % mod_data_type='AMSRE';
    % Can set these times to one, or do properly depending on what
    % time screening you want, etc.



    if iplot_wind_arrows==1


        nx_quiver=25; %number of arrows to draw for x
        ny_quiver=25; %number for y
        scale_speed_quiver = [15 15]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
        %scale_speed_quiver = [7 7]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
        %scale_speed_quiver = [10 10]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
        scale_speed_quiver = [1 1]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed
        scale_speed_quiver = [0.5 0.5]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed



        nx_quiver=100; %number of arrows to draw for x
        ny_quiver=100; %number for y

        nx_quiver = size(gcm_Plon2D_UM,2)/4; %All of th arrows divided by a factor
        ny_quiver = size(gcm_Plon2D_UM,1)/4;

        %nx_quiver=500; %number of arrows to draw for x
        %ny_quiver=500; %number for y

        %nx_quiver=25; %number of arrows to draw for x
        %ny_quiver=25; %number for y


        %[ax_orig]=plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

        %m_quiver(280,-20,U_glm(300,300),V_glm(300,300),0,'w');

        %m_quiver(gcm_Plon2D_UM,gcm_Plat2D_UM,U_glm,V_glm(1:end-1,:));



        %convert the lat and lon values to map coords
        [x_quiver,y_quiver]=m_ll2xy(gcm_Plon2D_UM,gcm_Plat2D_UM);
        %[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

        [ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');
        %[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,'k');


        %Plot the UM nest domain
        ioverride_box_colour=1;
        box_lwidth = 3;
        LAT_val = [-22.70 -17.28]; LON_val =[-78.93 -73.08]; %12th Nov UM FULL domain
        col_str='k-';
        itext_in_box=0;

    end

    plot_box_on_map


    dat_modis_save{it} = dat_modis;
    time_save(it) = time_round;

end %it loop