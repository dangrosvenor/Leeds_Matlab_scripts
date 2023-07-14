% Timeseries plot for UM
isave_plot_driver=0;
savedir_driver='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';

LAT_val = [-23.5 -16.44]; LON_val = [-85.93 -78.08]; %Smaller region to but out model edges
LAT_val = [-24.5 -15.44]; LON_val = [-86.93 -77.08]; %GOES region for UM comparison xkqk 26thOct POC
LAT_val = [-21 -16.44]; LON_val = [-82 -78.08]; %Much smaller region in NW corner
LAT_val = [-24.5 -15.44]; LON_val = [-84 -77.08]; %Same region as for the AGU PDFs - Trying to match AMSRE and GOES domains
LAT_val = [-21 -15.44]; LON_val = [-86.93 -82]; %Same as for revised PDFs post-AGU. Revised smaller region in NW corner  

%--- run the file to set up the defaults
watervap_defaults
idat_driver=0;

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'LWP timeseries';
xlab='Time (UTC)';
ylab='Liquid Water Path (g m^{-2})';
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.



idate_ticks_fix=1;
iaxis_square=0; %switch to make axis square


%--- Load and process the data
dirUM='/home/disk/eos8/d.grosvenor/UM/26thOct_POC/';

clear fileUM xdat_import ydat_import flag

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
fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = '(xkqkf) Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkh_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkh) 100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkj_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkj) 400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
% fileUM{idat} = 'xkqkk_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkk)
% 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqko_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkl_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkq_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
fileUM{idat} = 'xkqkr_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkm_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkn_LWP_RWP_.pp.nc'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;


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
    
    ydat_import(idat_driver).y = 1e3*meanNoNan(lwp(:,:),2);

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;        
    xdat_import(idat_driver).x = time_matlab; %
    
    labs_import(idat_driver).l = labs_UM(idat_UM).l;
end
%    xdat_import(idat_UM).x =



%---  Main script to do plots and save
savedir = savedir_driver;
isave_plot=isave_plot_driver;
DRIVER_lineplot_watervap

