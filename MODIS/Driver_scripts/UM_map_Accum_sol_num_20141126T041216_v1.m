% LWP map plots for all times for UM
% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
%LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

clear UM_map

UM_map.isave_plot=0;


UM_map.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
UM_map.thresh_LWP_mask = 20;

UM_map.iset_min_clim=0;
UM_map.clim_min=0;
UM_map.iset_max_clim=0;
UM_map.clim_max=300;

UM_map.varname = 'Accum_sol_num'; %The name of the variable in the nc file
UM_map.VAR_NAME_STR=''; %The name part in the filename to replace VAR_NAME with
            %OR, if the VAR_NAME is not in the supplied filename then will
            %just use the supplied one - make sure it is correct!
UM_map.flag2=''; %replacing the old flag for whether is a mat file or not - set to '' if not.
UM_map.var_units_str = 'kg^{-1}';

%UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_map.LAT_val_DRIVER = [-90 90]; UM_map.LON_val_DRIVER =[-180 180]; %Not needed unless restricting domain below
UM_map.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

%Choose a height index for 3D data if needed
UM_map.iz = 30; %=3km for Iceland case (can do z_um = nc{'hybrid_ht'}(:); to check (using previously opened nc file)

tol_mins=1;
time_single = datenum('01-Sep-2014 12:00'); 
%time_single = datenum('31-Aug-2014 03:00');
%time_single = datenum('01-Sep-2014 21:00');


%UM_map.time_range(1) = time_single - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;

%UM_map.time_range=[]; %set to thsi for all the tiems
UM_map.time_range(1) = time_single;

% -- For option setting see inside the loops



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/Iceland_Anja/';

clear fileUM xdat_import ydat_import
idat=1;
%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkh) 100cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkj) 400cm^{-3}'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;idat=idat+1;
%fileUM{idat} = 'xkmph_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkmph)'; pole_lat=70; pole_lon=284; idat=idat+1;

% - 12th Nov case, as of analysis started July 2015

%fileUM{idat} = '/xlhg-u/xlhgu_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlhg-v/xlhgv_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlhg-w/xlhgw_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
        
%fileUM{idat} = '/xlyd-x/xlydx_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Nd_fixed_act-10';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlyd-y/xlydy_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-new';  flag{idat} = 'load_mat'; fileUM_Nd{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlyd-p/xlydp_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-p/xlydp_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-z/xlydz_VAR_NAME_.pp.nc'; labs_UM(idat).l = 'CASIM-Ndvar-10-SWNd';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-z/xlydz_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0];
%marker_styleUM(idat).m='^'; idat=idat+1;

%UM_map.fileUM{idat} = '/xmmz-u/xmmzu_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'CASIM-Ndvar(increased)';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-u/xmmzu_rho_.pp.nc';UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-v/xmmzv_VAR_NAME_.pp.nc';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  UM_map.flag{idat} = 'load_mat'; UM_map.fileUM_rho{idat} = '/xmmz-v/xmmzv_rho_.pp.nc'; UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;
%UM_map.fileUM{idat} = '/xmmz-w/xmmzw_VAR_NAME_.pp.nc';  UM_map.labs_UM(idat).l ='CASIM-Ndvar(increased)-0.1';  UM_map.flag{idat} = 'load_mat';  UM_map.pole_lat=70; UM_map.pole_lon=284;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0 0 0]; UM_map.marker_styleUM(idat).m='v'; idat=idat+1;

%UM_map.fileUM{idat} = '/u-ad234/u-ad234_Accum_sol_mass_num_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad234'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
UM_map.fileUM{idat} = '/u-ad571/u-ad571_Accum_sol_mass_num_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad571, low aerosol, volcano ON'; UM_map.pole_lat=25; UM_map.pole_lon=165;
    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
UM_map.fileUM{idat} = '/u-ad572/u-ad572_Accum_sol_mass_num_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad572, low aerosol, volcano OFF'; UM_map.pole_lat=25; UM_map.pole_lon=165;
    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
  
    



UM_map.iplot_markers=0;
    
    
    
 UM_maps_20141126T041216_FUNC(UM_map)   
