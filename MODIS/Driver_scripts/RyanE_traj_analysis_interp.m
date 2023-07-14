%Uses RyanE_traj_analysis_interp_FUNC to do the actual interp to
%trajectories.

%Normal data (not anomalies from 100 day running average)
%save_file='/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_TEST.mat';
%RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,time_out,LWP_PD_ALL, RWP_PD_ALL,...
%TLWP_PD_ALL,Conv_LWP_PD_ALL,Conv_RWP_PD_ALL,low_CF_PD_ALL,total_CF_PD_ALL,Nd_PD_ALL,max_cloud_height_in_cloud_LWC_IWC_PD)

var_list={'totCF_COSP_MODIS_PD_ALL','liqCF_COSP_MODIS_PD_ALL','tot_surf_precip_rate_ALL','BL_height_dthdz_ALL','BL_height_combined_ALL','qv700_PD_ALL','LTS_PD_ALL','omega700_PD_ALL','surface_temp_PD_ALL','U10_PD_ALL','V10_PD_ALL',...
    'LWP_PD_ALL','RWP_PD_ALL','TLWP_PD_ALL','Conv_LWP_PD_ALL','Conv_RWP_PD_ALL','low_CF_PD_ALL','total_CF_PD_ALL','Nd_PD_ALL','max_cloud_height_in_cloud_LWC_IWC_PD'};




%% Split the data into day and nighttime based on the local time. 
%Do this based on being within a certain tolerance of 01:30 and 13:30 local
%time. Should make this consistent with the trajectories.
%Set time tolerance within which they have to match the trajetory time
time_tol_traj = 3.5; %in hours

[LT] = local_times_from_UTC(HH_UM,gcm_Plon2D_UM); %This is a 3D array of the same size as the model data
inot_tol_LT_day = find(abs(LT-13.5)>time_tol_traj);
inot_tol_LT_night = find(abs(LT-1.5)>time_tol_traj); 

for ivar=1:length(var_list)    
    var = var_list{ivar};
    %Copy the full arrays to start
    eval_str = [var '_day = ' var ';']; eval(eval_str);
    eval_str = [var '_night = ' var ';']; eval(eval_str);
    
    %Mask out the non-day and non-night values
    eval_str = [var '_day(inot_tol_LT_day) = NaN;']; eval(eval_str);
    eval_str = [var '_night(inot_tol_LT_night) = NaN;']; eval(eval_str);
       
%    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_out;']; eval(eval_str2);     
%    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var ';']; eval(eval_str2);    
end  



%% No Nan restriction for 100-day running mean
%Anomalies from 100 day running average
%calc the anomolies
nthresh=0; %Require a minimum of X non-NaN points within the 100-day running mean before being accepted.

dat_str = '_100_0_nonNaN_DAY'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 

for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    %Calculate the 100-day running mean and the anomalies for var
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_day,nthresh);'];eval(eval_str);            
    %Add time data to the mean and anomaly structures
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    %Add the 100-day mean and anomaly data to the mean and anomaly structure
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);


%night
dat_str = '_100_0_nonNaN_NIGHT'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 

for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_night,nthresh);'];eval(eval_str);            
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);



%% Save the standard variables (not as running means or anomalies)
dat_str = '_standard_DAY'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}        
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_out;']; eval(eval_str2);     
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_day;']; eval(eval_str2);    
end
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

dat_str = '_standard_NIGHT'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}        
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_out;']; eval(eval_str2);     
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_night;']; eval(eval_str2);    
end
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);


%% Require a minimum of 20 non-NaN points within the 100-day running mean before being accepted.
%Anomalies from 100 day running average)
%calc the anomolies
nthresh=20; %Require a minimum of XX non-NaN points within the 100-day running mean before being accepted.

dat_str = '_100_20_nonNaN_DAY'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_day,nthresh);'];eval(eval_str);            
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);

%night
dat_str = '_100_20_nonNaN_NIGHT'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_night,nthresh);'];eval(eval_str);            
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);

irun_100=0;
if irun_100==1

%% Require a minimum of 100 (i.e. all) non-NaN points within the 100-day running mean before being accepted.
%Anomalies from 100 day running average)
%calc the anomolies
nthresh=100; %Require a minimum of >=80 non-NaN points within the 100-day running mean before being accepted.

dat_str = '_100_100_nonNaN_DAY'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_day,nthresh);'];eval(eval_str);            
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);


%night
dat_str = '_100_100_nonNaN_NIGHT'
eval(['clear Ryan_dat' dat_str ' Ryan_dat_anom' dat_str]) 
for ivar=1:length(var_list)    
    var = var_list{ivar}
    
    eval_str=['[time_100, ' var '_100,' var '_anom_100] = RyanE_traj_analysis_100_day_anomalies(time_out,' var '_night,nthresh);'];eval(eval_str);            
    eval_str2 = ['Ryan_dat' dat_str '.time_100 = time_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.time_100 = time_100;']; eval(eval_str2);   
    eval_str2 = ['Ryan_dat' dat_str '.' var ' = ' var '_100;']; eval(eval_str2);
    eval_str2 = ['Ryan_dat_anom' dat_str '.' var ' = ' var '_anom_100;']; eval(eval_str2);   
end

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat' dat_str ');']; eval(eval_str);

save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];
eval_str=['RyanE_traj_analysis_interp_FUNC(save_file,gcm_Plat2D_UM,gcm_Plon2D_UM,Ryan_dat_anom' dat_str ');']; eval(eval_str);


%% Plot histograms of the number of complete trajectories (no NaNs in whole trajectory) from loaded data

dat_str = '_100_100_nonNaN_NIGHT';
dat_str = '_100_80_nonNaN_DAY';
dat_str = '_100_0_nonNaN_DAY';
%dat_str = '_100_0_nonNaN_NIGHT';

% Running mean data
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
% Anomalies
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];

var='Nd';
%var='omega700';
%var='LTS';
var='CTH'; %Has slightly more NaN trajectories than Nd
var='TLWP_inc_Conv';
%var='low_CF';
var='total_CF';
var='qv700';
var='surface_temp';
var='U10';
%var='V10';

%standard=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_standard.mat');
%anom=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom_100_80_nonNaN.mat');

dat = load(save_file,['traj_' var]);
eval(['inan=isnan(dat.traj_' var ');']);
%inan=isnan(traj_low_CF);
num_nans = sum(inan,1);
iwhole = find(num_nans==0);
Ntraj = length(iwhole);

qh=ndhistc_run([num_nans(:)],[-0.5:1:7.5]);
figure
plot([0:7],qh,'bx-');
xlabel('Number of NaNs in each trajectory');
ylabel('Frequency');
title(var);

%% Combined day and night trajectories for test
var='Nd';
%var='omega700';
%var='LTS';
%var='CTH'; %Has slightly more NaN trajectories than Nd
%var='TLWP_inc_Conv';
%var='low_CF';
%var='total_CF';
%var='qv700';
%var='surface_temp';
%var='U10';
%var='V10';

dat_pre_str = '_100_0_nonNaN_';
%dat_pre_str = '_100_20_nonNaN_';


%day
dat_str = [dat_pre_str 'DAY'];


% Running mean data
%save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
% Anomalies
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];

%standard=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_standard.mat');
%anom=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom_100_80_nonNaN.mat');

dat = load(save_file,['traj_' var]);
eval(['iday=isnan(dat.traj_' var ');']);
eval(['day=dat.traj_' var ';']);


%night
dat_str = [dat_pre_str 'NIGHT'];

% Running mean data
%save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data' dat_str '.mat'];
% Anomalies
save_file=['/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom' dat_str '.mat'];

%standard=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_standard.mat');
%anom=load('/home/disk/eos10/d.grosvenor/UM/RyanE_traj_analysis/UKMET_trajectories_model_data_anom_100_80_nonNaN.mat');

dat = load(save_file,['traj_' var]);
eval(['inight=isnan(dat.traj_' var ');']);
eval(['night=dat.traj_' var ';']);

%combined number of NaNs
i2=iday+inight;
day(iday)=0;
night(inight)=0;
both=day+night;
both(i2==2)=NaN;


inan = isnan(both);
%inan=isnan(traj_low_CF);
num_nans = sum(inan,1);
iwhole = find(num_nans==0);
Ntraj = length(iwhole);

qh=ndhistc_run([num_nans(:)],[-0.5:1:7.5]);
figure
plot([0:7],qh,'bx-');
xlabel('Number of NaNs in each trajectory');
ylabel('Frequency');
title(var);

end