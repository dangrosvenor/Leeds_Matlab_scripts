
%% Load data

%save(save_time_lat_file,'dom_mean_LWP_PI','dom_mean_LWP_PD','-APPEND','-V7.3');

um_case='u-ba458'; model_name='CASIM-bug'; %Old CASIM with sed bug
%um_case='u-ba333'; model_name='UKCA';%UKCA PI run (volcano OFF)
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
save_time_lat_file = [dirUM '/save_time_lat.mat'];
dat_ukca=load(save_time_lat_file);

%Select PI run (volcano OFF)
um_case='u-bc309'; model_name2='CASIM'%New CASIM runs with fixed sed. 
dirUM = ['/home/disk/eos10/d.grosvenor/UM/12Nov2008_Boutle/' um_case];
save_time_lat_file = [dirUM '/save_time_lat.mat'];
dat_casim=load(save_time_lat_file);

load(save_time_lat_file); %load the lat,lon, etc.


%% Timeseries of dLWP
%save(save_time_lat_file,'dom_mean_LWP_PI','dom_mean_LWP_PD','-APPEND','-V7.3');
    time_mean_dLWP = meanNoNan(dom_mean_LWP_PD-dom_mean_LWP_PI,1);
    
%    [indirect_ALL_timser] = UM_make_regional_timeseries(indirect_ALL,nT,LAT_val_DRIVER2,LON_val_DRIVER2,gcm_Plat2D_UM,gcm_Plon2D_UM,gcm_Plat2D_edges_UM,gcm_Plon2D_edges_UM);
        
    time = dat_global.time_ALL(time_inds);
    
    figure
    plot(time,dom_mean_LWP_PD-dom_mean_LWP_PI,'linewidth',3);
    datetick('x','dd');
    
    %legend(leg_strs,'location','northwest');
    xlabel('Time');
    %ylabel('SW surface forcing (W m^{-2})');
    ylabel('\DeltaLWP (g m^{-2})');
    set(gca,'ylim',[-2 8])
    fontsize_figure(gcf,gca,18);
    grid on

  
%% Overall RR means
overall_mean_ukca_OFF = meanNoNan(dat_ukca.mean_RR_PI(:),1)
overall_mean_ukca_ON = meanNoNan(dat_ukca.mean_RR_PD(:),1)
overall_mean_casim_OFF = meanNoNan(dat_casim.mean_RR_PI(:),1)
overall_mean_casim_ON = meanNoNan(dat_casim.mean_RR_PD(:),1)

dRR_ukca = 100*(overall_mean_ukca_ON - overall_mean_ukca_OFF)./overall_mean_ukca_OFF
dRR_casim = 100*(overall_mean_casim_ON - overall_mean_casim_OFF) ./ overall_mean_casim_OFF

%% Calcluate zonal means

N_lat_bins=100;
N_lat_bins=50;

[zonal_OFF_ukca,lats,Nout_OFF]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_ukca.mean_RR_PI,[],N_lat_bins);
[zonal_ON_ukca,lats,Nout_ON]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_ukca.mean_RR_PD,[],N_lat_bins);
[zonal_OFF_casim,lats,Nout_OFF]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_casim.mean_RR_PI,[],N_lat_bins);
[zonal_ON_casim,lats,Nout_ON]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_casim.mean_RR_PD,[],N_lat_bins);

%Discard latitudes where there aren't many data poits (due to curved grid)
Nmin=min(Nout_OFF,Nout_ON);
thresh_Nmin = median(Nmin)*0.8;
i=find(Nmin<thresh_Nmin);
lats(i)=NaN;
        
    
%% Plot zonal (over all lons) mean - volcano ON minus OFF for UKCA and
%% CASIM

leg_str={''};

figure
plot(3600*(zonal_ON_ukca - zonal_OFF_ukca),lats,'linewidth',3);
hold on
leg_str{1}=[model_name];

%figure
plot(3600*(zonal_ON_casim - zonal_OFF_casim),lats,'r','linewidth',3);
leg_str{2}=[model_name2];


legend(leg_str);
increase_font_size_map_figures  %Need to put the labels after for some reason...
xlabel('\DeltaRR (mm day^{-1}) (ON minus OFF)');
ylabel('Latitude (degrees)');
grid on

set(gca,'xlim',[-0.005 0.005]);



%% Plot zonal (over all lons) mean - actual RRs (4 lines)
leg_str={''}; ileg=1;

figure
hold on
plot(3600*zonal_ON_ukca,lats,'b-','linewidth',3);
plot(3600*zonal_OFF_ukca,lats,'b--','linewidth',3);

%leg_str{ileg}='UKCA ON'; ileg=ileg+1;
%leg_str{ileg}='UKCA OFF'; ileg=ileg+1;

leg_str{ileg}=[model_name ' ON']; ileg=ileg+1;
leg_str{ileg}=[model_name ' OFF'];ileg=ileg+1;



%figure
plot(3600*zonal_ON_casim,lats,'r','linewidth',3);
plot(3600*zonal_OFF_casim,lats,'r--','linewidth',3);

%leg_str{ileg}='CASIM ON'; ileg=ileg+1;
%leg_str{ileg}='CASIM OFF'; ileg=ileg+1;

leg_str{ileg}=[model_name2 ' ON']; ileg=ileg+1;
leg_str{ileg}=[model_name2 ' OFF']; ileg=ileg+1;

legend(leg_str,'Location','SouthEast');
increase_font_size_map_figures
%Need to put the labels after for some reason...
xlabel('Surface rain rate (mm day^{-1})');
ylabel('Latitude (degrees)');
grid on

set(gca,'xlim',[0 0.25]);



%% LWP

%% Overall LWP means
overall_mean_LWP_ukca_OFF = meanNoNan(dat_ukca.mean_LWP_PI(:),1)
overall_mean_LWP_ukca_ON = meanNoNan(dat_ukca.mean_LWP_PD(:),1)
overall_mean_LWP_casim_OFF = meanNoNan(dat_casim.mean_LWP_PI(:),1)
overall_mean_LWP_casim_ON = meanNoNan(dat_casim.mean_LWP_PD(:),1)

dLWP_ukca = 100*(overall_mean_LWP_ukca_ON - overall_mean_LWP_ukca_OFF)./overall_mean_LWP_ukca_OFF
dLWP_casim = 100*(overall_mean_LWP_casim_ON - overall_mean_LWP_casim_OFF) ./ overall_mean_LWP_casim_OFF

%% Calcluate zonal means

N_lat_bins=100;
N_lat_bins=50;

[zonal_LWP_OFF_ukca,lats_LWP,Nout_LWP_OFF]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_ukca.mean_LWP_PI,[],N_lat_bins);
[zonal_LWP_ON_ukca,lats_LWP,Nout_LWP_ON]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_ukca.mean_LWP_PD,[],N_lat_bins);
[zonal_LWP_OFF_casim,lats_LWP,Nout_LWP_OFF]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_casim.mean_LWP_PI,[],N_lat_bins);
[zonal_LWP_ON_casim,lats_LWP,Nout_LWP_ON]=UM_calc_lat_mean(gcm_Plat2D_UM,gcm_Plon2D_UM,dat_casim.mean_LWP_PD,[],N_lat_bins);

%Discard latitudes where there aren't many data poits (due to curved grid)
Nmin_LWP=min(Nout_LWP_OFF,Nout_LWP_ON);
thresh_Nmin = median(Nmin_LWP)*0.8;
i=find(Nmin_LWP<thresh_Nmin);
lats_LWP(i)=NaN;
   


%% Plot zonal (over all lons) mean - actual LWPs (4 lines)
leg_str={''}; ileg=1;

figure
hold on
plot(zonal_LWP_ON_ukca,lats_LWP,'b-','linewidth',3);
plot(zonal_LWP_OFF_ukca,lats_LWP,'b--','linewidth',3);

%leg_str{ileg}='UKCA ON'; ileg=ileg+1;
%leg_str{ileg}='UKCA OFF'; ileg=ileg+1;

leg_str{ileg}=[model_name ' ON']; ileg=ileg+1;
leg_str{ileg}=[model_name ' OFF']; ileg=ileg+1;

%figure
plot(zonal_LWP_ON_casim,lats_LWP,'r','linewidth',3);
plot(zonal_LWP_OFF_casim,lats_LWP,'r--','linewidth',3);

%leg_str{ileg}='CASIM ON'; ileg=ileg+1;
%leg_str{ileg}='CASIM OFF'; ileg=ileg+1;

leg_str{ileg}=[model_name2 ' ON']; ileg=ileg+1;
leg_str{ileg}=[model_name2 ' OFF']; ileg=ileg+1;

legend(leg_str,'Location','SouthWest');
increase_font_size_map_figures
%Need to put the labels after for some reason...
xlabel('LWP (g m^{-2})');
ylabel('Latitude (degrees)');
grid on

set(gca,'xlim',[40 120]);
    
%% Plot zonal (over all lons) mean - volcano ON minus OFF for UKCA and
%% CASIM

leg_str={''};

figure
plot((zonal_LWP_ON_ukca - zonal_LWP_OFF_ukca),lats_LWP,'linewidth',3);
hold on
%leg_str{1}='UKCA';
leg_str{1}=[model_name];

%figure
plot((zonal_LWP_ON_casim - zonal_LWP_OFF_casim),lats_LWP,'r','linewidth',3);
%leg_str{2}='CASIM';
leg_str{2}=[model_name2];

legend(leg_str,'Location','SouthEast');
increase_font_size_map_figures  %Need to put the labels after for some reason...
xlabel('\DeltaLWP (g m^{-2}) (ON minus OFF)');
ylabel('Latitude (degrees)');
grid on

%set(gca,'xlim',[-0.005 0.005]);



%% Plot zonal (over all lons) mean - UKCA - CASIM

leg_str={''};

figure
plot((zonal_LWP_OFF_ukca - zonal_LWP_OFF_casim),lats_LWP,'linewidth',3);
hold on
leg_str{1}='Volcano OFF';


%figure
plot((zonal_LWP_ON_ukca - zonal_LWP_ON_casim),lats_LWP,'r','linewidth',3);
leg_str{2}='Volcano ON';

legend(leg_str,'Location','SouthEast');
increase_font_size_map_figures  %Need to put the labels after for some reason...
xlabel('\DeltaLWP (g m^{-2}) (UKCA minus CASIM');
ylabel('Latitude (degrees)');
grid on









