% Timeseries plot for UM
isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

iplot_goes=0; %Whether to plot the GOES timeseries
iplot_RHB=0;  %Whether to plot the RHB ship timeseries
iplot_amsre=0;
ismooth_RHB=0;

%All based on the max in columns
dBZ_timeser_type = 'Median';
dBZ_timeser_type = 'Max';

%% N.B. - the REMSS and GOES files are based on a specific lat lon domain. So make sure that are using the same one for the
%% UM, and GOES (assuming that is what is required)

LAT_val_DRIVER = [-23.5 -16.44]; LON_val_DRIVER = [-85.93 -78.08]; %Smaller region to but out model edges
LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-86.93 -77.08]; %GOES region for UM comparison xkqk 26thOct POC
LAT_val_DRIVER = [-21 -16.44]; LON_val_DRIVER = [-82 -78.08]; %Much smaller region in NW corner
LAT_val_DRIVER = [-24.5 -15.44]; LON_val_DRIVER = [-84 -77.08]; %Same region as for the AGU PDFs - Trying to match AMSRE and GOES domains
LAT_val_DRIVER = [-21 -15.44]; LON_val_DRIVER = [-86.93 -82]; %Revised PDFs post-AGU. Revised smaller region in NW corner  

% 12th Nov case
%Need to come up with a new corner for the 12th Nov case since this is
%centred at 76W, 20S (not at 82W).
% Domain runs from -78.93 to -73.08 W and -22.70 to -17.28 (based on the edges)
LAT_val_DRIVER = [-22.70 -17.28]; LON_val_DRIVER = [-78.93 -73.08];

thresh_LWP_DRIVER = 5; %threshold for UM Nd - sets values with lower LWP to NaN, so that they aren't included in the mean

% Smaller region to account for boundary inflow of LWP and
% spin-up during advection.
%Looks like this mainly affects the south of the domain and to
%the east (for 26th Oct POC case the east was also affected).
%Also remove a bit for the boundary itself (around 0.25 deg
%should be enough).
LAT_val_DRIVER = [-20.5 -17.5]; LON_val_DRIVER = [-78.75 -73.25];

%The location of the ship = 20S, 75W
%One degree region around it
%LAT_val_DRIVER = [-20.5 -19.5]; LON_val_DRIVER = [-75.5 -74.5];


%LAT_val_DRIVER = [-20.5 -20.4]; LON_val_DRIVER = [-75.1 -75.0];

% See --- read_GOES_vocals_netcdf_files_MULTI.m ---
%           for reading in muliple files and saving the data in a .mat file
%           26-27th Oct
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150509T124841.mat';
%  12-14th Nov
%   LAT_val = [-22.70 -17.28]; LON_val = [-78.93 -73.08];
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150730T092051.mat';

%Using partial domain to avoid boundary issues
%   LAT_val= [-20.5 -17.5]; LON_val = [-78.75 -73.25];
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150731T024411.mat';
load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T064655.mat';
  %Special test case for 1 deg region around the ship to see if GOES
  %compares well with the aircraft :-
%load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20151112T074537.mat';  

thresh_SZA = 64;

%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150727T071106.mat';
%Using full domain
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T083601.mat';
%Using partial domain to avoid boundary issues
%load_file_remss = '/home/disk/eos1/d.grosvenor/UM/12Nov2008_Boutle/remss_lwp_saved_20150729T101619.mat';


%runt his script to get cdan(1:9) and markers(1:3) and pdan(1:3)
LInestyles_etc

%--- run the file to set up the defaults
watervap_defaults
for idat=1:99
    ismooth_x_import(idat)=0;
    ismooth_y_import(idat)=0;
end
idat_driver=0;

clear fileUM xdat_import ydat_import flag 

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = [dBZ_timeser_type 'of column max dBZ timeseries'];
xlab='Time (UTC)';
xlab='Time (Local Solar Time)';
ylab='Radar Reflectivity (dBZ)';

% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 12th Nov this was at 16:48 for
% -20, -76 lat lon (centre of the domain). I.e. they are 4hrs 48 mins behind UTC
time_shift = -(4+48/60) /24; %amount to shift time by for LST (from UTC)
xlims=1;
xlimits=[datenum('00:00 12-Nov-2008') datenum('00:00 16-Nov-2008')] + time_shift; %shift to LST since final plot will be in LST
xlimits=[datenum('00:00 12-Nov-2008') datenum('00:00 14-Nov-2008')] + time_shift; %shift to LST since final plot will be in LST

izlim=0;
zmin=1500;
zmax=3000;

lor=-1; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.



idate_ticks_fix=1;
iaxis_square=0; %switch to make axis square

idat_micro=[];
if iplot_amsre==1
%% ------------------------------
% ------ AMSRE, SSMI and TMI data --------
% ------------------------------


%Load the data file with the remss data for the specific model domain
%region - made using AMSRE_LWP_maps_20150717T124800.m
load(load_file_remss);

        
%------- Calculate the data to plot         
         
%          % For use later to coarsen the UM data
%          d=abs(diff(gcm_Plat2D_AMSRE,[],1));
%          dlat_AMSRE = meanNoNan(meanNoNan(d,1),1);
%          d=abs(diff(gcm_Plon2D_AMSRE,[],2));
%          dlon_AMSRE = meanNoNan(meanNoNan(d,1),1);
%          
%          %Making a "new" gcm_str to avoid using timeseries3 type data in
%          %pdf2D
%          gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
%          gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
%          gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
%          gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
%          %Actually can't really define this it varies over the globe for
%          %ASMRE
% %         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre); 
%             gcm_time_matlab_AMSRE2 = 0;
%             gcm_time_UTC_AMSRE2 = 0;
%             daynum_timeseries3_AMSRE2 = 1;
%             modisyear_timeseries3_AMSRE2 = 1;
%          
%         
% 
% % ----- Set various things
% 
%           
%          
% %        mod_data_type='AMSRE';
%         gcm_str_select='AMSRE2';
%         gcm_str='AMSRE2';
       
%        month_amsre = goes_month;
%        year_amsre = goes_year;


for isat=1:length(sats)
    idat_driver=idat_driver+1;
    idat_micro(isat)=idat_driver; %store the indices for the microwave lines for joining together later
    
    itimser2=0;
    
    [lc,lp,mark]=choose_linestyles_func(idat_driver,cdan,pdan,markers);
    line_pattern(idat_driver) = lp;  line_colour(idat_driver) = lc; marker_style(idat_driver)=mark;    
    
    
for itimser=1:size(lwp_remss,3)
    for idaynight=2:-1:1
        
        itimser2=itimser2+1;
    
    
%     Y_driver = 1e3*squeeze(lwp_amsre(:,:,itimser,idaynight));
%     Y_driver = 1e3*squeeze(lwp_remss(:,:,itimser,idaynight,isat));     
        
%         %--- run the file to set up the defaults
% %        plot_global_maps_defaults   
%          pdf2D_defaults
%          
%         
%         %--- set some options for these particular plot loops
% %        set_screening = {'none'};
% %        modis_data_plot = 'Map of 2D data from outside driver script';
%         i577 = 'MODIS_plot_UW';
% 
%         iset_min_clim=1;
%         clim_min=0;
%         iset_max_clim=1;
%         clim_max=200;
%         
%         logflag=0;
%         dlogflag=0;
%         
%         isave_plot=0;
% %        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
%         
%         
%                         
% 
% %        screen_type = 'gcm_screening';
% 
%         %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
%         x_axis_vals = 'Dummy data'; %dummy data
%         y_axis_vals = 'General GCM-style';
%         y_axis_vals =  'General y-axis no ilat simple';
%         
%         ylabelstr='LWP (g m^{-2})';
%         
% %        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;
%         Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;
%         
%                                
%                                 
% %          logbin_norm = logbin_norm_driver;
% %          i_plot_norm=i_plot_norm_driver;
% %          i_div_bin_widths=i_div_bin_widths_driver;
% %          pdf_type = pdf_type_driver;
%                                 
% %        gcm_str = gcm_str_last_loaded;        
% 
%         
%  % --------- Override flags for 2D PDF --------
%         ioverride_pdf=1;
%         %iocean_only=1;
%         man_choose_plotTimeHeight_graph=1;
%         ioverride_location_selection=1;
%         ioverride_pdf_varchoose = 1;
%         datatype = 'gcm_data';        
% 
%         % --------- Override flags for watervap --------
%         man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
%        plotTimeHeightVap3
%        close(gcf);
%        waterVapourMay2005
%        close(gcf);

% ydat_import(idat_driver).y(itimser2) = Y_mean_overall; lwp_mean_remss
 ydat_import(idat_driver).y(itimser2) = lwp_mean_remss(itimser,idaynight,isat);  

%  year = str2num(chosen_files(itimser).name(7:10));
%  month = str2num(chosen_files(itimser).name(11:12));
%  day = str2num(chosen_files(itimser).name(13:14)); 
%  
%  %Aqua overpass is 13:30 or 01:30 local time
%  if idaynight==1
%      hour = mod( 13.5 - time_shift*24 , 24); %time_shift is in days
%  else
%      hour = mod( 1.5 - time_shift*24 , 24 );          
%  end
 
 %time = datenum(year,month,day,hour,0,0);
 
 time = date_remss(itimser,idaynight,isat) + timeUTC_remss(itimser,idaynight,isat)/24;
 
 xdat_import(idat_driver).x(itimser2) = time + time_shift;

 
 
    end


end

labs_import(idat_driver).l = sats{isat};
inan = find(isnan(xdat_import(idat_driver).x)==1);
xdat_import(idat_driver).x(inan)=[];
ydat_import(idat_driver).y(inan)=[];

end

   

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;   


end
    
    



%% ------------------------------
% ------ GOES data --------
% ------------------------------
if iplot_goes==1

    idat_driver=idat_driver+1;

    load(load_file_goes);  %Values stored in a .mat file previously - used read_GOES_vocals_netcdf_files_MULTI.m 
    

    xdat_import(idat_driver).x = times_GOES_save + time_shift;
    ydat_import(idat_driver).y = goes_Nd_mean;

    labs_import(idat_driver).l = 'GOES';

    sza = sun_pos(times_GOES_save,mean(LAT_val_DRIVER),mean(LON_val_DRIVER));
    xdat_import(idat_driver).x(sza>thresh_SZA) = NaN;
    

line_pattern(idat_driver).p= '-';  line_colour(idat_driver).c=[0.6 0.6 0.8]; marker_style(idat_driver).m='o';


end





%% -------- UM data -----------------------------------------------------

%--- Load and process the data
dirUM='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/';

for idat=1:99
   flag{idat} = 'nc';
   fileUM_rho{idat} = ''; 
end

% idat=1;
% fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkh) 100cm^{-3}';  idat=idat+1;
% fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkj) 400cm^{-3}';  idat=idat+1;
% fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; ; idat=idat+1;
% fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; idat=idat+1;
% fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7';  idat=idat+1;
% fileUM{idat} = 'xkqkq_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; idat=idat+1;
% fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_import(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme'; idat=idat+1;
% fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_import(idat).l = '(xkqkf) 1000cm^{-3} No cloud-scheme'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; idat=idat+1;

idat=1;
%fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = '(xkqkf) Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkh) 100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkj) 400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
% fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkk)
% 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkq_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
%fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkm_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkn_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;

% 
% %New runs as done for AGU
% fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = 'No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1; %1000cm^{-3} 
% 
% %fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkx_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
% 
% %fileUM{idat} = 'xkqkm_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
% %fileUM{idat} = 'xkqkn_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
% 
% fileUM{idat} = 'xkqkv_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkw_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;

% - 12th Nov case, as of analysis started July 2015
fileUM{idat} = '/xlhg-u/xlhgu_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
    line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
% High and low aerosol cases.
fileUM{idat} = '/xlhg-v/xlhgv_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-0.1';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
fileUM{idat} = '/xlhg-w/xlhgw_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284;
        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.8 0.4 0]; marker_styleUM(idat).m='^'; idat=idat+1;
%fileUM{idat} = '/xlyd-m/xlydm_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-aproc';  flag{idat} = 'load_mat'; pole_lat=70; pole_lon=284; 
%        line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0]; marker_styleUM(idat).m='^'; idat=idat+1;

%fileUM{idat} = '/xlyd-x/xlydx_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Nd-fixed-act-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-x/xlydx_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.4 0.4]; marker_styleUM(idat).m='d'; idat=idat+1;
%fileUM{idat} = '/xlyd-y/xlydy_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-y/xlydy_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 0]; marker_styleUM(idat).m='v'; idat=idat+1;
%fileUM{idat} = '/xlyd-q/xlydq_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-fixed-act-NoRain';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-q/xlydq_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[1 0.0 0.0]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xlyd-r/xlydr_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-fixed-act-NoRain-Less-Mixing';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-r/xlydr_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0.4 0.8 0.0]; marker_styleUM(idat).m='<'; idat=idat+1;
%fileUM{idat} = '/xlyd-s/xlyds_Nd_.pp.nc.mat'; labs_UM(idat).l = 'CASIM-Ndvar-10-NoRain-Less-Mixing';  flag{idat} = 'load_mat'; fileUM_rho{idat} = '/xlyd-s/xlyds_rho_.pp.nc'; pole_lat=70; pole_lon=284;
%line_patternUM(idat).p= '--';  line_colourUM(idat).c=[0 0 1]; marker_styleUM(idat).m='<'; idat=idat+1;

        

%Bypass this for now
for idat_UM=1:length(fileUM)
%for idat_UM=1:length(fileUM)
    idat_driver = idat_driver+1;
    
    
    
    filename = [dirUM fileUM{idat_UM}];
    filename_rho = [dirUM fileUM_rho{idat_UM}];
    

    line_pattern(idat_driver)=line_patternUM(idat_UM);  line_colour(idat_driver)=line_colourUM(idat_UM); marker_style(idat_driver)=marker_styleUM(idat_UM);
    

%    lwp = get_LWP_RWP_UM('LWP',flag{idat_UM},filename,filename_rho);
%     clear vars_in
%     vars_in.var = 'LWP';
%     vars_in.flag = flag{idat_UM};
%     vars_in.file_lwp = filename;
%     vars_in.file_rho = filename_rho;
%     vars_in.pole_lat = pole_lat;
%     vars_in.pole_lon = pole_lon;
%     vars_in.time_in = [];
% 
%     [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
%     
    
         vars_in.var = 'dBZ';
         vars_in.flag = flag{idat_UM};
         vars_in.file_lwp =  [dirUM remove_character(fileUM{idat_UM},'Nd','dBZ')];   %Keep as file_lwp for loading .mat files
%         vars_in.file_rho = [dirUM fileUM_rho{idat_UM}]; %filename_rho;
         vars_in.pole_lat = pole_lat;
         vars_in.pole_lon = pole_lon;
         vars_in.time_in = []; %set to this for all times
    
     [dBZ_3D,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it_driver,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

     dBZ_2D = max(dBZ_3D,[],2);



    stime=size(dBZ_2D,1);
    [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val_DRIVER,LON_val_DRIVER,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
    dBZ_2D=permute(dBZ_2D,[2 3 1]);
    a=NaN*ones(size(dBZ_2D));
    a(iregion_lin)=0;
    dBZ_2D=dBZ_2D+a;
    dBZ_2D=permute(dBZ_2D,[3 1 2]);
    
% vars_in.var = 'RWP';
%     vars_in.flag = flag{idat_UM};
%     vars_in.file_lwp = filename;
%     vars_in.file_rho = filename_rho;
%     vars_in.pole_lat = pole_lat;
%     vars_in.pole_lon = pole_lon;
%     vars_in.time_in = [];

%     [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
% 
%     stime=size(rwp,1);
%     [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val_DRIVER,LON_val_DRIVER,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
%     rwp=permute(rwp,[2 3 1]);
%     a=NaN*ones(size(rwp));
%     a(iregion_lin)=0;
%     rwp=rwp+a;
%     rwp=permute(rwp,[3 1 2]); 
    
%    rwp=0;
    
%    tlwp = lwp+rwp;
    
switch dBZ_timeser_type
    case 'Median'
        %    ydat_import(idat_driver).y = meanNoNan(Nd(:,:),2);
        ydat_import(idat_driver).y = prctile(dBZ_2D(:,:),50,2);
    case 'Max'
           ydat_import(idat_driver).y = max(dBZ_2D(:,:),[],2); 
end

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;   


%    xdat_import(idat_driver).x = time_matlab; %
    xdat_import(idat_driver).x = time_matlab + time_shift; %    
    
    labs_import(idat_driver).l = labs_UM(idat_UM).l;
end


%% -------- RH Brown (ship) LWP data -----------------------------------------------------
if iplot_RHB==1
    
idat_driver = idat_driver+1;

%--- Load and process the data
fileRHB='/home/disk/eos8/d.grosvenor/UM/12Nov2008_Boutle/VOCALS_obs_from_IanB/vocals_rhb_lwp10min_v1.nc';
ncRHB=netcdf(fileRHB);
 

time=ncRHB{'jd_jdfrac'}(:); %"julian day + day fraction, 10-min res,time marker at beginning"
time(time>1e30)=NaN;
%convert into matlab time
time_matlab_RHB = datenum(2008,1,1) -1 + time;

lwp_RHB = ncRHB{'lwp'}(:);
lwp_RHB(lwp_RHB<-98.9)=NaN;
lwp_RHB(lwp_RHB>1e30)=NaN;

ydat_import(idat_driver).y = lwp_RHB;
xdat_import(idat_driver).x = time_matlab_RHB + time_shift; %
labs_import(idat_driver).l = 'RHB ship';

inan = find(isnan(xdat_import(idat_driver).x)==1);
xdat_import(idat_driver).x(inan)=[];
ydat_import(idat_driver).y(inan)=[];

if ismooth_RHB==1
    ismooth_y_import(idat_driver)=1;
    smooth_mode='mean';
    Nsmooth_window = 6;
end


end










%% ---  Main script to do plots and save
savedir = savedir_driver;
%
ichoose_styles=1;
DRIVER_lineplot_watervap

for i=idat_micro   %length(h)-2
set(h(i).h,'linestyle','none');
end

% set(h(1).h,'marker','*');
% set(h(1).h,'markersize',20);
% uistack(h(1).h,'top');

%Joint together the REMSS satellite values with a line.
x_all=[]; y_all=[];
for i=idat_micro  %1:length(xdat)-2
   x_all = cat(2,x_all,xdat(i).x);
   y_all = cat(2,y_all,ydat(i).y);   
end
[x_all,I]=sort(x_all);
y_all=y_all(I);

plot(x_all,y_all,'b','linewidth',3);

% for i=idat_micro(end)+1:length(xdat)
%     set(h(i).h,'color',line_colour(i).c);             
%     set(h(i).h,'marker',marker_style(i).m,'markerEdgeColor',line_colour(i).c,'markerFaceColor',line_colour(i).c);
%     set(h(i).h,'linestyle',line_pattern(i).p);
% 
% end


isave_plot=isave_plot_driver;
if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1);
end
