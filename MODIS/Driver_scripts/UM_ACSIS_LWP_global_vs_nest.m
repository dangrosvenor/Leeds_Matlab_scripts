%Pre-process using 
%   UM_quick_plot_global.m

%Add variables to :-
%   UM_var_defs.m

isave_plot_global_diff=1;
iplot_wind_arrows=0;

%If want to plot an outline of the nest on the global map
filename_nest = '/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/u-as103/BC__model_level01_time08.mat';
nest=load(filename_nest);
irotated_pole_box=1;

%% File for times for the loop
um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; %global

var_UM_DRIVER = 'LWP';
var_UM_DRIVER = 'accum_number_ukca';
var_UM_DRIVER = 'SW_down_surf';

var_UM = var_UM_DRIVER;
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
filename = [dirUM '/' run_type '_' var_UM '_ALL.mat'];
load(filename,'time_ALL');




for it_global_diff=1:length(time_ALL)
    
    time_round = time_ALL(it_global_diff);
    time_format_str = 'UTC';

    scrsz=get(0,'ScreenSize');
    posit=[9 60 scrsz(3) scrsz(4)];
    figure('position',posit);



%% global run
    subplot(1,3,1);
    var_UM = var_UM_DRIVER;

    um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; %global
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    filename = [dirUM '/' run_type '_' var_UM '_ALL.mat'];

    if it_global_diff==1
        dat_global = load(filename)
    end

    dat_modis = eval(['dat_global.' var_UM '_ALL{it_global_diff};']);
    dat_modis_global = dat_modis;
    gcm_Plat2D_edges_UM = dat_global.gcm_Plat2D_edges_UM;
    gcm_Plat2D_UM = dat_global.gcm_Plat2D_UM;
    gcm_Plon2D_edges_UM = dat_global.gcm_Plon2D_edges_UM;
    gcm_Plon2D_UM = dat_global.gcm_Plon2D_UM;


    %run plotting script
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands
    
    
%% Nest
    subplot(1,3,2);
    var_UM = var_UM_DRIVER;    

    um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umnsaa'; %nest
    %um_case='u-at459'; pole_lat=45.0; pole_lon=145.0; run_type = 'umglaa'; %global
    dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
    filename = [dirUM '/' run_type '_' var_UM '_ALL.mat'];

    if it_global_diff==1
        dat_nest = load(filename)
    end

    dat_modis = eval(['dat_nest.' var_UM '_ALL{it_global_diff};']);
    dat_modis_nest = dat_modis;
    gcm_Plat2D_edges_UM = dat_nest.gcm_Plat2D_edges_UM;
    gcm_Plat2D_UM = dat_nest.gcm_Plat2D_UM;
    gcm_Plon2D_edges_UM = dat_nest.gcm_Plon2D_edges_UM;
    gcm_Plon2D_UM = dat_nest.gcm_Plon2D_UM;

    

    %run plotting script
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands
    

%% Diff
    subplot(1,3,3);
    var_UM = [var_UM_DRIVER ' diff'];
    
    dat_modis = dat_modis_nest - dat_modis_global;   

    %run plotting script
    UM_ACSIS_LWP_global_vs_nest_quick_plot_commands
    

%% Save    
    savename=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS/global_vs_LWP_plots/' num2str(it_global_diff,'%02d') '_global_' titlenam_driver];

    if isave_plot_global_diff==1
        clear opts
%        opts.iplot_png=1;
        opts.iplot_eps=1;
        saveas_ps_fig_emf(gcf,[savename],'',0,1,0,'',[],0,opts);
        close(gcf);
    end



end