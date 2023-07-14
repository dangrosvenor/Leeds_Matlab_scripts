%Stripped down template for running plot_global_maps to produce a global
%map

% Set some stuff

var_UM = 'wind speed';
var_UM = 'RH';
var_UM = 'potemp';
%var_UM = 'W';

%Data to plot
it=10; %04:30 UTC 12Nov2008
it=87; %19:00 UTC 13 Nov
it=79; %15:00 UTC 13 Nov
%it=81; %16:00 UTC 13 Nov
%it=83; %17:00 UTC 13 Nov
%it=85; %18:00 UTC 13 Nov

iz=6; %133m
%iz=9; %300m
iz=12; %513m

iplot_wind_arrows


%Define lat and lon grids
filename = '/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_UVW_component_3D_.pp.nc';
nc = netCDF(filename);
%filenameV = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-ar365/global/u-ar365_V_200811131500___glm_saved_slice_.pp.nc';
%ncV = netCDF(filenameV);

gcm_str_DRIVER='UM';  %N.B. - have to set it to this and not another name - doesn't work if do (code would need changing)


global_or_nest = 'nest';
switch global_or_nest
    case 'nest'
        lat = nc{'y'}(:);
        lon = nc{'x'}(:);
        [lon2d,lat2d] = meshgrid(lon,lat);

        pole_lat=70; pole_lon=284;

        [gcm_Plat2D_UM,gcm_Plon2D_UM]=em2gm(lat2d,lon2d,pole_lat,pole_lon);

    case 'global'
        gcm_Plat2D_UM = nc{'Latitude'}(:);
        gcm_Plon2D_UM = nc{'Longitude'}(:);

        i180 = find(gcm_Plon2D_UM>180);
        gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;
end

% i180 = find(gcm_Plon2D_UM>180);
% gcm_Plon2D_UM(i180) = gcm_Plon2D_UM(i180) - 360;

[gcm_Plon2D_edges_UM,gcm_Plat2D_edges_UM] = get_edges_lat_lon(gcm_Plon2D_UM,gcm_Plat2D_UM);

t0 = datenum('12-Nov-2008 00:00');
time_glm = nc{'t'}(:);
% time_gml_matlab = datenum('01-Jan-1970') + time_glm/24;
time_gml_matlab = t0 + time_glm;

time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
time_LST = time_gml_matlab(it) + time_shift;

[Y,M,D,HH,MM,SS]=datevec(time_LST);
if SS>59
    SS=0;
    MM=MM+1;
end
time_round = datenum(Y,M,D,HH,MM,SS);


 
daynum_timeseries3_UM = 1; %[1:length(time)];
gcm_time_UTC_UM = 1; %[1:length(time)];



%get data for wind vectors
U_glm = nc{'U'}(it,iz,:,:);
V_glm = nc{'V'}(it,iz,:,:);
V_glm(end+1,:) = V_glm(end,:);

titlenam_driver = [var_UM ' at ' datestr(time_round) ' LST'];
units_str_plot = '';

switch var_UM
    case 'wind speed'

        %dat_modis = V_glm;
        %dat_modis = U_glm;
        dat_modis = sqrt(U_glm.^2 + V_glm.^2);        
        
        titlenam_driver = ['Wind speed (m s^{-1}) at ' datestr(time_round) ' LST'];
        units_str_plot = '';
        
    case 'W'
        dat_modis = nc{'W'}(it,iz,:,:);
    case 'qL'
        filename2 = '/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_qL_.pp.nc';
        nc2=netcdf(filename2)
        dat_modis = nc2{'qL'}(it,iz,:,:);
        
    case 'potemp'
        filename2 = '/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_th_qv_.pp.nc';
        nc2=netcdf(filename2)
        dat_modis = nc2{'th'}(it,iz,:,:);     
        
    case 'RH'
        filename2 = '/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_th_qv_.pp.nc';
        nc2=netcdf(filename2);
        th = nc2{'th'}(it,iz,:,:);
        qv = nc2{'qV'}(it,iz,:,:);
        
        filename3 = '/home/disk/eos15/d.grosvenor/UM/12Nov2008_Boutle/xmmz-u/xmmzu_rho_.pp.nc';
        nc3=netcdf(filename3);
        rho = nc3{'rho'}(it,iz,:,:);
        
        R = 6371000; %Radius of the Earth!
        rho = rho/R/R;
        
        [T,P]=T_P_from_th_rho(th,rho);
        
        dat_modis = calc_RH(qv,P,T);
        
end


%%


LAT_val_DRIVER = [-1e9 1e9];
LON_val_DRIVER = [-1e9 1e9]; %Can restrict the plot region with this if set irestrict_domain_DRIVER=1

LAT_val_DRIVER = [-40 0];
LON_val_DRIVER = [-100 -60]; 

LAT_val_DRIVER = [-24 -16];
LON_val_DRIVER = [-80 -72]; 


irestrict_domain_DRIVER=0;

%--- run the file to set up the defaults
        plot_global_maps_defaults
        
irestrict_domain=irestrict_domain_DRIVER; %whether to restrict the domain or not

thresh_LAT = LAT_val_DRIVER;
thresh_LON = LON_val_DRIVER;

%--- set some options for these particular plot loops
set_screening = {'none'};
modis_data_plot = 'Map of 2D data from outside driver script';

iset_min_clim=0;
clim_min=0;
iset_max_clim=0;
clim_max=64800;

isave_plot=0;
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

savename=['/home/disk/eos1/d.grosvenor/modis_work/lwp_vs_Nd_plots/nested_' titlenam_driver];

if isave_plot==1            
    saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],1); %last option is whether to do eps+pdf plot - can take a while
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
    scale_speed_quiver = [5 5]; %max speed expected - arrows are scaled according to this speed - i.e. 15 m/s would produce the biggest arrow that looks decent on the plot - so that plots with different wind speeds will produce arrows that consistently proportional to the wind speed

    %scale_speed_quiver = [0.003 0.003];

    nx_quiver=100; %number of arrows to draw for x
    ny_quiver=100; %number for y

    nx_quiver=25; %number of arrows to draw for x
    ny_quiver=25; %number for y


    %[ax_orig]=plot_wind_quiver_arrows(u_quiver,v_quiver,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

    %m_quiver(280,-20,U_glm(300,300),V_glm(300,300),0,'w');

    %m_quiver(gcm_Plon2D_UM,gcm_Plat2D_UM,U_glm,V_glm(1:end-1,:));



    %convert the lat and lon values to map coords
    [x_quiver,y_quiver]=m_ll2xy(gcm_Plon2D_UM,gcm_Plat2D_UM);
    %[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');

    [ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,scale_speed_quiver,1,1,isquare,'k');
    %[ax_orig]=plot_wind_quiver_arrows(U_glm,V_glm,x_quiver,y_quiver,nx_quiver,ny_quiver,'k');

end