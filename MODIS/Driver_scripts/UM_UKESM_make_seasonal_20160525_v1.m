% Loops through all 4 seasons for all the years and makes a climatology for
% multiple models. Saves in Nd_UKESM_savefile in the form of 
%  Name              Size              Bytes  Class     Attributes
% 
%   Nd_UKESM          3x4             2655552  cell                
%   Ndatap_UKESM      3x4             2655552  cell                
%   UM_map            1x1                5346  struct 
%
% Where Nd_UKESM{imodel,iseason} is a cell array for the different models
% selected.

Nd_UKESM_savefile = '/home/disk/eos8/d.grosvenor/UM/UKESM/Nd_UKESM_seasonal_clim_2000-2007.mat';

%Nd_file_UKESM = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/MODIS_TERRA_Nd_QA_from_L3_daily_UKESM_Jane_2000_2014_allCF_CTT273_region_mask_20150416T063530.mat';



clear UM_map

UM_map.isave_plot=0;

iplot_maps=1;
plot_diff=0;
plot_pdf_diff=0;

UM_map.varname = 'Nd_zweight'; %The name of the variable in the nc file
UM_map.VAR_NAME_STR='Nd_GCM'; %The name part in the filename to replace VAR_NAME with
            %OR, if the VAR_NAME is not in the supplied filename then will
            %just use the supplied one - make sure it is correct!
            
%UM_map.flag2=''; %replacing the old flag for whether is a mat file or not - set to '' if not.
UM_map.flag2='load_mat'; %replacing the old flag for whether is a mat file or not - set to '' if not.

UM_map.var_units_str = 'cm^{-3}';

UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov
UM_map.LAT_val_DRIVER = [-22.70 -17.28]; UM_map.LON_val_DRIVER =[-78.93 -73.08]; %FULL UM domain for 12th Nov

UM_map.irestrict_domain_DRIVER=0; %set to zero to plot the full domain

UM_map.time_tol='same month'; %set to this to match by month and year

%Setting time in the time loop below
%tol_mins=1;
%time_single = datenum('01-May-2007 00:00'); 
%UM_map.time_range(1) = time_single;  % - tol_mins/24/60; UM_map.time_range(2) = time_single + tol_mins/24/60;


UM_map.i_mask_low_LWP=0; %Make any values below thresh_LWP equal to NaN
UM_map.thresh_LWP_mask = 20;


%% Options for PDF plotting

Ybins_DRIVER = [0:20:800];
logbin_norm_driver = 0; %Whether hav lognormal bins and normalise by these
i_plot_norm_driver = 1; %whether to normalise
i_div_bin_widths_driver = 1; %Whether to divide by the bin widths (also normalises)

pdf_type_driver='normal';  %normal or cumulative PDF
%pdf_type_driver='cumulative';



%--- Load and process the data
%dirUM='/home/disk/eos1/d.grosvenor/UM/26thOct_POC/';
%UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';
%UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/Iceland_Anja/';
UM_map.dirUM='/home/disk/eos8/d.grosvenor/UM/UKESM/';

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

%UM_map.fileUM{idat} = '/u-ad234/u-ad234_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad234'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
%UM_map.fileUM{idat} = '/u-ad571/u-ad571_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad571, low aerosol, volcano ON'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
%UM_map.fileUM{idat} = '/u-ad572/u-ad572_SWLW_TOA_outgoing_.pp.nc'; UM_map.labs_UM(idat).l = 'u-ad572, low aerosol, volcano OFF'; UM_map.pole_lat=25; UM_map.pole_lon=165;
%    UM_map.line_patternUM(idat).p= '--';  UM_map.line_colourUM(idat).c=[0.4 0.4 0.4]; UM_map.marker_styleUM(idat).m='d'; idat=idat+1;
    
UM_map.fileUM{idat} = '/u-ab642/u-ab642_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA7-N96 AMIP'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;
UM_map.fileUM{idat} = '/u-ab754/u-ab754_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA6-N96 AMIP'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;
UM_map.fileUM{idat} = '/u-ac043/u-ac043_VAR_NAME_.pp.nc'; UM_map.labs_UM(idat).l = 'GA7-N96 ORCA1.0'; UM_map.pole_lat=70; UM_map.pole_lon=284; idat=idat+1;

    
    
UM_map.iset_min_clim=1;
UM_map.clim_min=0;
UM_map.iset_max_clim=1;
UM_map.clim_max=650;


UM_map.iplot_markers=0;

UM_map.noplot=1; %suppress the plotting part - just read in the data.

% time_single = datenum(['01-Jul-2004 00:00']); datestr(time_single)
% UM_map.time_range(1) = time_single;
% vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);
% 


%% Loop over time for creating seasonal data

% Loop over time - create seasonal climatologies
% Requeseted the files from 2000-2008,.
% At first hadn't used the 360 day caldendar, so the dates looked weird.
% Now that I have used 360 day calendar, we have data for 16th of every month
% from 2000-2007 (Jan to Dec). Not sure what happened to 2008 - think the
% model stopped at start of 2008 since said ran from 1982-2008.
% Need to decide what to do RE DJF ,etc. Prob is ok to use Jan and Feb 2000
% and Dec 2007 since am just doing a climatology
% (7 years). Should be enough for Nd, which doesn't very that much I think
% - can check this.
% So, just using DJF of each year, e.g. Dec, Jan and Feb of 2000.

years_clim = [2000:2007];
im=1; clear months
months{im}='Dec'; im=im+1;
months{im}='Jan'; im=im+1;
months{im}='Feb'; im=im+1;
months{im}='Mar'; im=im+1;
months{im}='Apr'; im=im+1;
months{im}='May'; im=im+1;
months{im}='Jun'; im=im+1;
months{im}='Jul'; im=im+1;
months{im}='Aug'; im=im+1;
months{im}='Sep'; im=im+1;
months{im}='Oct'; im=im+1;
months{im}='Nov'; im=im+1;



Ndatap = length(years_clim)*3; %3 months for each year

for iseason_clim=1:4
    ifirst=1; %whether we are starting a new season for averaging
    for iyear_clim=1:length(years_clim)
        for imonth_clim=1:3
%             %Special case for Dec of each year - want to use the current
%             %year, whereas for all other months will use the next year
%             if iseason_clim==1 & imonth_clim==1
%                 year = years_clim(iyear_clim);
%             else
%                 year = years_clim(iyear_clim+1);                 
%             end

            year = years_clim(iyear_clim);
            im = (iseason_clim-1)*3 + imonth_clim;
            time_single = datenum(['01-' months{im} '-' num2str(year) ' 00:00']); datestr(time_single)
            UM_map.time_range(1) = time_single;        
            vars_map_out = UM_maps_20141126T041216_FUNC(UM_map);            
            
            for imodel=1:3
                if ifirst==1
%                    dat_season{imodel}=zeros(size(vars_map_out.P_save{imodel}));
                    dat_season{imodel}=zeros([length(years_clim)*3 size(vars_map_out.P_save{imodel})]);                    
                end
%                dat_season{imodel} = dat_season{imodel} + vars_map_out.P_save{imodel} / Ndatap;
                i_ind = (iyear_clim-1)*3 + imonth_clim;
%                dat_season{imodel}(iyear_clim,imonth_clim,:,:) = vars_map_out.P_save{imodel};                
                dat_season{imodel}(i_ind,:,:) = vars_map_out.P_save{imodel};                  
            end
            ifirst=0;
        end 
        
        
    end
    
    for imodel=1:3
        Nd_UKESM{imodel,iseason_clim} = dat_season{imodel};
        [Nd_UKESM{imodel,iseason_clim}, Ndatap_UKESM{imodel,iseason_clim}] = meanNoNan(dat_season{imodel},1);        
        %cut out data if there were any months with NaN for that location
        %(essentially the same as doing a straight mean)
        i_cut = find(Ndatap_UKESM{imodel,iseason_clim} < length(years_clim)*3);
        Nd_UKESM{imodel,iseason_clim}(i_cut) = NaN;
    end
            
end
  
clear UM_grid
UM_grid.gcm_Plat2D_edges_UM = vars_map_out.gcm_Plat2D_edges_UM;
UM_grid.gcm_Plon2D_edges_UM = vars_map_out.gcm_Plon2D_edges_UM;
UM_grid.gcm_Plat2D_UM = vars_map_out.gcm_Plat2D_UM;
UM_grid.gcm_Plon2D_UM = vars_map_out.gcm_Plon2D_UM;

save(Nd_UKESM_savefile,'Nd_UKESM','Ndatap_UKESM','UM_map','UM_grid');






