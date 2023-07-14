% 4th Jan, 2017 - Nd data from mock L3 data for Matt C - JJA 2008 only.
% Have saved the data from the MODIS* script described below here for quicker loading using less memory :-
% '/home/disk/eos8/d.grosvenor/mat_files_various/MattC_JJA2008_mockL3data/SAVED_JJA2008_CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65.mat'
% Load this, or...
% First of all run MODIS_multi_DAY_processL3L2.m to load in the mock L3
% data. Settings :-
%  multiL2L3_case = 'load L3 and concatenate';
%  multiL2L3_project='global';
%  daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/'];
%  direcs={'aqua'};  years=[2008]; days={[153:244]};  %JJA for 2008

isave=1;
%savedir = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/';
savedir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/';
savedir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_meanCTT_173_meanCTH_3.2km_SZA_65/';
save_name = 'Nd_re_tau_saved_JJA2008_no_CF_screening';

date_str = datestr(now,30);
savefile01 = [savedir save_name '_' date_str '.mat'];

date_tag = 'JJA_2008';
var_tags = {'Nd_16','Nd_21','Nd_37','re_16','re_21','re_37','tau'};

iscreen_seaice=1; %
iscreen_land=0;
iscreen_region=0; %Whether to apply the regional screening (same for all years, but different for each month - i.e. all Jans use the same one)
iscreen_sza=0; % This refers to L3 screening (i.e. based on daily averaged data)
%max_sza_thresh = 65;

iscreen_ndays=1;
%Exclude points below a threshold number of days over the whole season
thresh_ndays=5;



Date_text = 'For JJA 2008 (days 153 to 244)';
Notes =[...
    'Global droplet concentrations from MODIS using technique described in Grosvenor & Wood, ACP,2014, doi:10.5194/acp-14-7291-2014. '...
    'Data is given for the 1.6, 2.1 and 3.7um MODIS near_infrared bands. '...
    'Individual MODIS Aqua Collection 5.1 L2 swaths with 1km Nd data were averaged into 1x1 degree regions and then into daily data. Swaths with SZA >65 degrees are excluded. ' ...  
    'Daily data is then averaged into seasonal values. ' ...      
    'To be included in the seasonal average there must have been at least '  num2str(thresh_ndays) ' days of data. '...
    'Fill values are set to NaN. Standard deviations of the seasonal values are also given based on the inital daily data. '...
    'The total number of days of data going into each average are also given. '...
    'Cloud Top Height (CTH) is limited to <3.2km. CTH was calculated based on the 1x1 degree mean CTT and SST using the method described in Zuidema, 2009 (10.1175/2009jcli2708.1)'...
    '(Might want to update to 1km CTH at some point)'...
    'Data is also restricted to 1x1 deg regions with a liquid cloud fraction that is >80%.'...
    'Created with Matlab script MattC_mockL3_Nd_dataset_JJA2008.m'...
    ];





%% Combine NH and SH seaice (is only SH in the file)
%seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2007_all_NHemisphere_20150529T103046.mat';
seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2008_JJA_all_NHemisphere_20161222T222241.mat';
seaiceNH = load(seaice_fileload_1deg_NH);
seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2008_JJA_all_SHemisphere_20161221T095015.mat';
seaiceSH = load(seaice_fileload_1deg_SH);

%Do the time matching Load these to do the time matching (already loaded in
%MODIS* script)
%load(file01,'modisyear_timeseries3','daynum_timeseries3');
%seaice_match_times  %made a function version of this (below)
[seaice_time3_NH,seaice_max_time3_NH] = seaice_match_times_FUNC(seaiceNH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
[seaice_time3_SH,seaice_max_time3_SH] = seaice_match_times_FUNC(seaiceSH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
%So will now have NaN values in both SH and NH files. Will set the NaNs to zero for
%now and then make NaN the places where both were NaN.
inanNH = isnan(seaice_time3_NH); seaice_time3_NH(inanNH)=0;
inanSH = isnan(seaice_time3_SH); seaice_time3_SH(inanSH)=0;
inanNH_max = isnan(seaice_max_time3_NH); seaice_max_time3_NH(inanNH_max)=0;
inanSH_max = isnan(seaice_max_time3_SH); seaice_max_time3_SH(inanSH_max)=0;

iboth = inanNH & inanSH; %When both are NaN - need to use the logical since they are not just integers
iboth_max = inanNH_max & inanSH_max; %When both are NaN - need to use the logical since they are not just integers

seaice_time3 = seaice_time3_NH + seaice_time3_SH;
seaice_time3(iboth) = NaN;

seaice_max_time3 = seaice_max_time3_NH + seaice_max_time3_SH;
seaice_max_time3(iboth_max) = NaN;

%save(file01,'-APPEND','-V7.3','seaice_time3','seaice_max_time3');

%days of month funciton will be useful to give days required for each month
%- e.g. 32:60 for Feb.
% [days_required_for_mean,time_mean_str] = days_of_month(month_no)




clear idays Nd_16 Nd_21 Nd_37 sea_ice_max_tmp sea_ice_tmp re_16 re_21 re_37 tau


Nd_16 = N_time3_16;
Nd_21 = N_time3;
Nd_37 = N_time3_37;
re_16 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3;
re_21 = Cloud_Effective_Radius_Liquid_Mean.timeseries3;
re_37 = Cloud_Effective_Radius_16_Liquid_Mean.timeseries3;
tau = Cloud_Optical_Thickness_Liquid_Mean.timeseries3;

siz=size(Nd_16);

if iscreen_seaice==1
    
    sea_ice = seaice_time3;
    %sea_ice = seaice_max_time3;    %max over 2 week window
    
    %find points where there is any sea-ice using the threshold below
    sea_ice_thresh = 1e-5;
    %    sea_ice_thresh = 1e9;
    isea_high = find(sea_ice>sea_ice_thresh);
    %    isea_high = find(sea_ice_max>sea_ice_thresh);
    Nd_16(isea_high)=NaN;
    Nd_21(isea_high)=NaN;
    Nd_37(isea_high)=NaN;
    re_16(isea_high)=NaN;
    re_21(isea_high)=NaN;
    re_37(isea_high)=NaN;
    tau(isea_high)=NaN;
    
end


%% Make an overall mean for all days of the season

for ivar=1:length(var_tags)
    var_tag = var_tags{ivar};
    
    [temp,temp_Ndays,temp_std_dev] = eval(['meanNoNan(' var_tag ',3);']);
    
    if iscreen_ndays==1
        %Exclude data where there are fewer than thresh_nyears in the
        %climatology
        inan = find(temp_Ndays<thresh_ndays);
        temp(inan) = NaN;
        temp_std_dev(inan) = NaN;
        
    end
    
    eval_str = [var_tag '_' date_tag ' = temp;']; eval(eval_str);
    eval_str = [var_tag '_' date_tag '_Ndays = temp_Ndays;']; eval(eval_str);
    eval_str = [var_tag '_' date_tag '_std_dev = temp_std_dev;']; eval(eval_str);
    
    
    
    
    
    
    
    
    
    
    if iscreen_land==1  %This needs changing to work for mock L3
        
        landmask_load= load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
        landmask = repmat(flipdim(landmask_load.amsre_land_mask,1),[1 1 12]);
        %Smooth out the land mask to make it more aggressive - seemed to be
        %allowing some small bits of land otherwise.
        %Also need to do this for the sea ice
        
        %What is halo?
        
        fs=fspecial('average',mask_size); %Make the averaging filter
        landmask_tmp = add_halo(landmask_load.amsre_land_mask,halo);
        landmask_smooth = filter2(fs,landmask_tmp,'same');
        landmask_smooth = remove_halo(landmask_smooth,halo);
        landmask2 = repmat(flipdim(landmask_smooth,1),[1 1 12 15]);
        
        
        %Land mask is zero where there is land and NaN where not - so can just add to the
        %array to make the non-land regions NaN
        
        %Make the land regions NaN
        %    Nd_16_ocean = Nd_16 + landmask2;
        %    Nd_21_ocean= Nd_21 + landmask2;
        %    Nd_37_ocean = Nd_37 + landmask2;
        %    std_Nd_16_ocean = std_Nd_16 + landmask2;
        %    std_Nd_21_ocean= std_Nd_21 + landmask2;
        %    std_Nd_37_ocean = std_Nd_37 + landmask2;
        
        eval_str = ['Nd' var_tag ' = Nd' var_tag ' + landmask2;']; eval(eval_str);
        eval_str = ['Nd_Ndays' var_tag ' = Nd_Ndays' var_tag ' + landmask2;']; eval(eval_str);
        eval_str = ['Nd_std_dev' var_tag ' = Nd_std_dev' var_tag ' + landmask2;']; eval(eval_str);
        
    end
    
    
    
    
    clear vars_to_save
    i=1;
    vars_to_save{i}='Date_text'; i=i+1;
    vars_to_save{i}='Notes'; i=i+1;
    %vars_to_save{i}='Notes2'; i=i+1;
    vars_to_save{i}='MLAT'; i=i+1;
    vars_to_save{i}='MLON'; i=i+1;
    %         vars_to_save{i}='Nd_16_ocean'; i=i+1;
    %         vars_to_save{i}='Nd_21_ocean'; i=i+1;
    %         vars_to_save{i}='Nd_37_ocean'; i=i+1;
    %         vars_to_save{i}='std_Nd_16_ocean'; i=i+1;
    %         vars_to_save{i}='std_Nd_21_ocean'; i=i+1;
    %         vars_to_save{i}='std_Nd_37_ocean'; i=i+1;
    %         vars_to_save{i}='Nd_16'; i=i+1;
    %         vars_to_save{i}='Nd_21'; i=i+1;
    %         vars_to_save{i}='Nd_37'; i=i+1;
    %         vars_to_save{i}='std_Nd_16'; i=i+1;
    %         vars_to_save{i}='std_Nd_21'; i=i+1;
    %         vars_to_save{i}='std_Nd_37'; i=i+1;
    %         vars_to_save{i}='ndays_16'; i=i+1;
    %         vars_to_save{i}='ndays_21'; i=i+1;
    %         vars_to_save{i}='ndays_37'; i=i+1;
    
%     vars_to_save{i}=['Nd' var_tag '_clim_ANNUAL']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_DJF']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_MAM']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_JJA']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_SON']; i=i+1;
    
    vars_to_save{i}=[var_tag '_' date_tag]; i=i+1;
    
%     vars_to_save{i}=['Nd' var_tag '_clim_std_dev_ANNUAL']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_std_dev_DJF']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_std_dev_MAM']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_std_dev_JJA']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_std_dev_SON']; i=i+1;
    

      vars_to_save{i}=[var_tag '_' date_tag '_std_dev']; i=i+1;

%     vars_to_save{i}=['Nd' var_tag '_clim_Ndays_ANNUAL']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_Ndays_DJF']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_Ndays_MAM']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_Ndays_JJA']; i=i+1;
%     vars_to_save{i}=['Nd' var_tag '_clim_Ndays_SON']; i=i+1;
    
    vars_to_save{i}=[var_tag '_' date_tag '_Ndays']; i=i+1;
    
    %Save here
    if isave==1
        filename_savevars = savefile01;
        tag='';
        save_vars_mat
    end
    
end

if isave==1
    %Convert to netCDF
    mat2nc_Dan(filename_savevars,[filename_savevars '.nc']);
%    mat2nc_Dan_matlab_netcdf(filename_savevars,[filename_savevars '.nc']);
end
