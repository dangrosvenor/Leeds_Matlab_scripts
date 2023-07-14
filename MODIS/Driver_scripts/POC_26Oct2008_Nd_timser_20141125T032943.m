% Timeseries plot for UM
savedir_driver ='/home/disk/eos1/d.grosvenor/modis_work/plots/UM/';
isave_driver=1;

%--- run the file to set up the defaults
watervap_defaults
idat_driver=0;

%--- set some options for this particular plot
graph=0; %graph choice in watervap
titlenam = 'Nd timeseries for columns with LWP.GT.10 g m^{-2} and z.LTE.2km';
xlab='Time (UTC)';
ylab='Domain mean of column max N_d (# mg^{-1})';
xlims=0;
xlimits=[0 100];

izlim=0;
zmin=1500;
zmax=3000;

lor=4; %1=top-right, 2=top-left 0=auto, 3=lower left, 4=lower right, -1=on the right of plot pane
nmark=-1; %-1 means that all points have markers. Otherwise only plot the number specified.

isave_plot=0;

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
%%fileUM{idat} = 'xkqkf_qL_qR_.pp.nc.mat'; labs_UM(idat).l = '(xkqkf) Old-mphys'; flag{idat}='load_mat'; fileUM_rho{idat} = 'xkqkf_rho_.pp.nc'; pole_lat=70; pole_lon=278; idat=idat+1;
%%fileUM{idat} = 'xkqkh_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkh) 100cm^{-3} RHcrit=0.8'; pole_lat=70; pole_lon=278; idat=idat+1;
% fileUM{idat} = 'xkqkj_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkj) 400cm^{-3}';  pole_lat=70; pole_lon=278;idat=idat+1;
% fileUM{idat} = 'xkqkk_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkk) 400cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkl'; labs_UM(idat).l = '(xkqkl) 1000cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqko'; labs_UM(idat).l = '(xkqko) 100cm^{-3} RHcrit=0.7'; pole_lat=70; pole_lon=278; idat=idat+1;
%fileUM{idat} = 'xkqkq_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkq) 100cm^{-3} No cloud-scheme'; pole_lat=70; pole_lon=278;idat=idat+1;
%fileUM{idat} = 'xkqkr_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkr) 1000cm^{-3} No cloud-scheme';pole_lat=70; pole_lon=278; idat=idat+1;
%%fileUM{idat} = 'xkqkm_Nd_.pp.nc'; labs_UM(idat).l = '(xkqkm) 1000cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;
fileUM{idat} = 'xkqkn'; labs_UM(idat).l = '(xkqkn) 100cm^{-3} RHcrit=0.7 AeroProc';pole_lat=70; pole_lon=278; idat=idat+1;

clear vars_UM vars_in

idat=1;
vars_UM{idat}='LWP_RWP'; idat=idat+1;
vars_UM{idat}='Nd'; idat=idat+1;
vars_UM{idat}='rho'; idat=idat+1;

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
    

    vars_in.var = 'Nd';
    vars_in.flag = flag{idat_UM};
%    vars_in.file_lwp = filename;
%    vars_in.file_rho = filename_rho;
    vars_in.pole_lat = pole_lat;
    vars_in.pole_lon = pole_lon;
    vars_in.time_in = [];

    [Nd,time_matlab,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM,it,daynum_timeseries3_UM,modisyear_timeseries3_UM,gcm_time_UTC_UM,gcm_time_matlab_UM] = get_LWP_RWP_UM(vars_in);
    

    
    ydat_import(idat_driver).y = 1e-6*Nd;

%    time=nc{'t'}(:);
%    t0_str=nc{'t'}.time_origin{1};
%    t0_str2=[t0_str(1:11) ' ' t0_str(13:17)];
%    xdat_import(idat_UM).x = datenum(t0_str2) + time;        
    xdat_import(idat_driver).x = time_matlab(1:length(ydat_import(idat_driver).y)); %
    
    labs_import(idat_driver).l = labs_UM(idat_UM).l;
end
%    xdat_import(idat_UM).x =



%---  Main script to do plots and save
isave_plot = isave_driver;
savedir = savedir_driver;

DRIVER_lineplot_watervap




                            

                            
                            

                            
                
                    
                    
                                    

                                    
                                    
                                

                          
                         
                            
                            
                         

        
