% Timeseries plot for UM
isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

LAT_val = [-23.5 -16.44]; LON_val = [-85.93 -78.08]; %Smaller region to but out model edges
LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08]; %GOES region for UM comparison xkqk 26thOct POC
LAT_val = [-21 -16.44]; LON_val = [-82 -78.08]; %Much smaller region in NW corner
LAT_val = [-24.5 -15.44]; LON_val = [-84 -77.08]; %Same region as for the AGU PDFs - Trying to match AMSRE and GOES domains
LAT_val = [-21 -15.44]; LON_val = [-86.93 -82]; %Revised PDFs post-AGU. Revised smaller region in NW corner  

load_file_goes = '/home/disk/eos8/d.grosvenor/VOCALS/GOES_cloud/cloud-products/saved_multiple_days_20150509T124841.mat';

%--- run the file to set up the defaults
watervap_defaults
idat_driver=0;

clear fileUM xdat_import ydat_import flag

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'LWP timeseries';
xlab='Time (UTC)';
xlab='Time (Local Solar Time)';
ylab='Liquid Water Path (g m^{-2})';

% Shift to local time (Local Solar Time - so will base this on the time at
% which the Sun is highest in the sky. On 26th Oct this was at 17:32 for
% -15, -87 lat lon and 17:12 for at -82 lon, or 17:22 for the centre of
% these two lons. I.e. they are 5hrs 12 mins behind UTC
time_shift = -(5+12/60) /24; %amount to shift time by for LST (from UTC)
xlims=1;
xlimits=[datenum('12:00 25-Oct-2008') datenum('21:00 26-Oct-2008')] + time_shift; %shift to LST since final plot will be in LST

izlim=0;
zmin=1500;
zmax=3000;

lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.



idate_ticks_fix=1;
iaxis_square=0; %switch to make axis square


%% ------------------------------
% ------ AMSRE data --------
% ------------------------------
idat_driver=idat_driver+1;

        
%------- Calculate the data to plot         
         
         % For use later to coarsen the UM data
         d=abs(diff(gcm_Plat2D_AMSRE,[],1));
         dlat_AMSRE = meanNoNan(meanNoNan(d,1),1);
         d=abs(diff(gcm_Plon2D_AMSRE,[],2));
         dlon_AMSRE = meanNoNan(meanNoNan(d,1),1);
         
         %Making a "new" gcm_str to avoid using timeseries3 type data in
         %pdf2D
         gcm_Plat2D_AMSRE2 = gcm_Plat2D_AMSRE;
         gcm_Plat2D_edges_AMSRE2 = gcm_Plat2D_edges_AMSRE;
         gcm_Plon2D_AMSRE2 = gcm_Plon2D_AMSRE;
         gcm_Plon2D_edges_AMSRE2 = gcm_Plon2D_edges_AMSRE;
         %Actually can't really define this it varies over the globe for
         %ASMRE
%         gcm_time_matlab_AMSRE2 = datenum(year_amsre,month_amsre,day_amsre); 
            gcm_time_matlab_AMSRE2 = 0;
            gcm_time_UTC_AMSRE2 = 0;
            daynum_timeseries3_AMSRE2 = 1;
            modisyear_timeseries3_AMSRE2 = 1;
         
        

% ----- Set various things

          
         
%        mod_data_type='AMSRE';
        gcm_str_select='AMSRE2';
        gcm_str='AMSRE2';
       
%        month_amsre = goes_month;
%        year_amsre = goes_year;

itimser2=0;
for itimser=1:size(lwp_amsre,3)
    for idaynight=2:-1:1
        
        itimser2=itimser2+1;
    
    
     Y_driver = 1e3*squeeze(lwp_amsre(:,:,itimser,idaynight));
        
        %--- run the file to set up the defaults
%        plot_global_maps_defaults   
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
        
        isave_plot=0;
%        savedir='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
        
        
                        

%        screen_type = 'gcm_screening';

        %                            x_axis_vals = 'LWP+RWP GCM grid-box mean'; %dummy data
        x_axis_vals = 'Dummy data'; %dummy data
        y_axis_vals = 'General GCM-style';
        
        ylabelstr='LWP (g m^{-2})';
        
%        Ybins = [-0.01 30:10:2500]; ichoose_Ybins=1;
        Ybins = [-0.01 10.^[log10(30):0.1:log10(2500)]]; ichoose_Ybins=1;
        
                               
                                
%          logbin_norm = logbin_norm_driver;
%          i_plot_norm=i_plot_norm_driver;
%          i_div_bin_widths=i_div_bin_widths_driver;
%          pdf_type = pdf_type_driver;
                                
%        gcm_str = gcm_str_last_loaded;        

        
 % --------- Override flags for 2D PDF --------
        ioverride_pdf=1;
        %iocean_only=1;
        man_choose_plotTimeHeight_graph=1;
        ioverride_location_selection=1;
        ioverride_pdf_varchoose = 1;
        datatype = 'gcm_data';        

        % --------- Override flags for watervap --------
        man_choose_water_graph=1;    %for watervap 
        
        %---  Run plot script and save
        plotTimeHeightVap3
        close(gcf);
%        waterVapourMay2005
%        close(gcf);

 ydat_import(idat_driver).y(itimser2) = Y_mean_overall;

 year = str2num(chosen_files(itimser).name(7:10));
 month = str2num(chosen_files(itimser).name(11:12));
 day = str2num(chosen_files(itimser).name(13:14)); 
 
 %Aqua overpass is 13:30 or 01:30 local time
 if idaynight==1
     hour = mod( 13.5 - time_shift*24 , 24); %time_shift is in days
 else
     hour = mod( 1.5 - time_shift*24 , 24 );          
 end
 
 time = datenum(year,month,day,hour,0,0);
 xdat_import(idat_driver).x(itimser2) = time + time_shift;
 
    end


end

   

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;   



    
    labs_import(idat_driver).l = 'AMSRE';



%% ------------------------------
% ------ GOES data --------
% ------------------------------
idat_driver=idat_driver+1;

load(load_file_goes);
xdat_import(idat_driver).x = times_GOES_save + time_shift;
ydat_import(idat_driver).y = goes_LWP_mean;
    
labs_import(idat_driver).l = 'GOES';

sza = sun_pos(times_GOES_save,mean(LAT_val),mean(LON_val));
xdat_import(idat_driver).x(sza>65) = NaN;



%% -------- UM data -----------------------------------------------------

%--- Load and process the data
dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/longer_timeseries/';

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


%New runs as done for AGU
fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = 'No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1; %1000cm^{-3} 

%fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkx_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;

%fileUM{idat} = 'xkqkm_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkn_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;

fileUM{idat} = 'xkqkv_LWP_RWP_.pp.nc'; labs_UM(idat).l = '100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkw_LWP_RWP_.pp.nc'; labs_UM(idat).l = '10cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;


for idat_UM=1:length(fileUM)
    idat_driver = idat_driver+1;
    
    filename = [dirUM fileUM{idat_UM}];
    filename_rho = [dirUM fileUM_rho{idat_UM}];
    

%    lwp = get_LWP_RWP_UM('LWP',flag{idat_UM},filename,filename_rho);
    clear vars_in
    vars_in.var = 'LWP';
    vars_in.flag = flag{idat_UM};
    vars_in.file_lwp = filename;
    vars_in.file_rho = filename_rho;
    vars_in.pole_lat = pole_lat;
    vars_in.pole_lon = pole_lon;
    vars_in.time_in = [];

    [lwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    stime=size(lwp,1);
    [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
    lwp=permute(lwp,[2 3 1]);
    a=NaN*ones(size(lwp));
    a(iregion_lin)=0;
    lwp=lwp+a;
    lwp=permute(lwp,[3 1 2]);
    
vars_in.var = 'RWP';
    vars_in.flag = flag{idat_UM};
    vars_in.file_lwp = filename;
    vars_in.file_rho = filename_rho;
    vars_in.pole_lat = pole_lat;
    vars_in.pole_lon = pole_lon;
    vars_in.time_in = [];

    [rwp,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);

    stime=size(rwp,1);
    [iregion_lin,iregion_lin_edges,iregion_lin2D,iregion_lin2D_edges] = restrict_to_region_2D_lats(LAT_val,LON_val,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,stime);
    rwp=permute(rwp,[2 3 1]);
    a=NaN*ones(size(rwp));
    a(iregion_lin)=0;
    rwp=rwp+a;
    rwp=permute(rwp,[3 1 2]);    
    
    tlwp = lwp+rwp;
    
    ydat_import(idat_driver).y = 1e3*meanNoNan(tlwp(:,:),2);

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;   


%    xdat_import(idat_driver).x = time_matlab; %
    xdat_import(idat_driver).x = time_matlab + time_shift; %    
    
    labs_import(idat_driver).l = labs_UM(idat_UM).l;
end
%    xdat_import(idat_UM).x =



%% ---  Main script to do plots and save
savedir = savedir_driver;
%
DRIVER_lineplot_watervap


set(h(1).h,'linestyle','none');
set(h(1).h,'marker','*');
set(h(1).h,'markersize',20);
uistack(h(1).h,'top');


isave_plot=isave_plot_driver;
if isave_plot==1
    saveas_ps_fig_emf(gcf,[savename],'',0,1);
end
