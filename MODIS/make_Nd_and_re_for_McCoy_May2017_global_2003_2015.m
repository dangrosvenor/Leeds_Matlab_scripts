% This script automatically loads the data using MODIS_multi_DAY_processL3L2 with make_mockL3_variables_script_name = 'make_mockL3_variables_just_for_Nd_May2017';
% and multiL2L3_case = 'load L3 and concatenate'; for each year.
% Need to manually set the directory to load from in
% MODIS_multi_DAY_processL3L2.m (daily_averaged_files_loc2) and the action
% required (multiL2L3_case= 'load L3 and concatenate')

%Already have the Droplet_Number variables (made using 1km tau and reff).
%And 1x1 deg values (made in filtering_data_get) - N_time3 and N_time3_37
%McCoy would prefer Nd21 and 37 in separate files formatted as
%[lon,lat,time] with a time vector of size [N 3] (Y M D) and lat and lon
%(vectors I guess?)

iscreen_seaice=1;
iseaice_max_2week=1; %Whether to screen with max sea-ice over two week window or just daily values
iscreen_land=0;



CF_str = '0.8'; CF_str2 = '80';
CF_str = '0.0'; CF_str2 = '0'; %N.B. - don't currently have the cloud fraction for these...

file_dir=['/home/disk/eos15/d.grosvenor/mock_L3/CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];
file_dir2=['/home/disk/eos1/d.grosvenor/mock_L3/CF_' CF_str '_meanCTT_173_meanCTH_3.2km_SZA_65/'];

direcs={'aqua'};  years_multi=[2003:2015]; days_multi={[1:365],[1:366],[1:365],[1:365],[1:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365],[1:365]}; %,[1:365]};%

for iyear_Nd_dataset = 1:length(years_multi)
    
    years = years_multi(iyear_Nd_dataset); year_str = num2str(years);
    days = days_multi(iyear_Nd_dataset);
    %Run the script to load the data
    i_override_MODIS_multi=1;
    MODIS_multi_DAY_processL3L2
 

    
[date_str,date_num] = date_from_day_of_year_func(daynum_timeseries3,modisyear_timeseries3);
time = datevec(date_num); time=time(:,1:3);
lat=LAT;
lon=LON;

Notes='Produced using make_Nd_for_McCoy_May2017_global_2003_2015.m';

seaice_str = '';

%load and match sea-ice data
if iscreen_seaice==1
    
    seaice_fileload_1deg_NH = ['/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_' year_str '_all_NHemisphere_2018-07-18.mat'];
    seaiceNH = load(seaice_fileload_1deg_NH);
    seaice_fileload_1deg_SH = ['/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_' year_str '_all_SHemisphere_2018-07-18.mat'];
    seaiceSH = load(seaice_fileload_1deg_SH);
    
    seaice_str = '_screened_for_seaice_';
    

    
    %Do the time matching Load these to do the time matching (already loaded in
    %MODIS* script)
    %seaice_match_times  %made a function version of this (below)
    [seaice_time3_NH,seaice_max_time3_NH] = seaice_match_times_FUNC(seaiceNH,Cloud_Effective_Radius_37_Liquid_Mean,modisyear_timeseries3,daynum_timeseries3);
    [seaice_time3_SH,seaice_max_time3_SH] = seaice_match_times_FUNC(seaiceSH,Cloud_Effective_Radius_37_Liquid_Mean,modisyear_timeseries3,daynum_timeseries3);
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
    
    
%% Now screen data
    if iseaice_max_2week==0
        sea_ice = seaice_time3;
        seaice_max_2week_str = '';
    else
        sea_ice = seaice_max_time3;    %max over 2 week window
        seaice_max_2week_str = '_2week_max';
    end
    %find points where there is any sea-ice using the threshold below
    sea_ice_thresh = 1e-5;
    %    sea_ice_thresh = 1e9;
    isea_high = find(sea_ice>sea_ice_thresh);
    %    isea_high = find(sea_ice_max>sea_ice_thresh);
    
    N_time3(isea_high)=NaN;
    N_time3_37(isea_high)=NaN;
    Droplet_Number_Concentration.timeseries3(isea_high)=NaN;
    Droplet_Number_Concentration_37.timeseries3(isea_high)=NaN;
    Cloud_Optical_Thickness_Liquid_Mean.timeseries3(isea_high)=NaN;
    Cloud_Effective_Radius_Liquid_Mean.timeseries3(isea_high)=NaN;
    Cloud_Effective_Radius_37_Liquid_Mean.timeseries3(isea_high)=NaN;
    
    
end

if iscreen_land==1  %
    
    seaice_max_2week_str = [seaice_max_2week_str '_screened_for_land'];

    landmask_load= load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
%    landmask = repmat(flipdim(landmask_load.amsre_land_mask,1),[1 1 12]);
    %Smooth out the land mask to make it more aggressive - seemed to be
    %allowing some small bits of land otherwise.
    %Also need to do this for the sea ice

    %Halo adds some a buffer region to the edges, so we don't lose data
    %when smoothing
    
    mask_size=2;
    halo = mask_size*3;

    fs=fspecial('average',mask_size); %Make the averaging filter
    landmask_tmp = add_halo(landmask_load.amsre_land_mask,halo);
    landmask_smooth = filter2(fs,landmask_tmp,'same');
    landmask_smooth = remove_halo(landmask_smooth,halo);
    landmask2 = repmat(flipdim(landmask_smooth,1),[1 1 size(N_time3,3)]);


    %Land mask is zero where there is land and NaN where not - so can just add to the
    %array to make the non-land regions NaN
    
    vars = {'N_time3','N_time3_37','Droplet_Number_Concentration.timeseries3','Droplet_Number_Concentration_37.timeseries3','Cloud_Optical_Thickness_Liquid_Mean.timeseries3',...
        'Cloud_Effective_Radius_Liquid_Mean.timeseries3','Cloud_Effective_Radius_37_Liquid_Mean.timeseries3'};
    
    for ivar=1:length(vars)
        eval_str = [vars{ivar} ' = ' vars{ivar}  '+ landmask2;']; eval(eval_str);
    end
    
   

end



%2.1um
Nd = permute(N_time3,[2 1 3]);
Nd_1km = permute(Droplet_Number_Concentration.timeseries3,[2 1 3]);
tau_1km = permute(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,[2 1 3]);
re_1km = permute(Cloud_Effective_Radius_Liquid_Mean.timeseries3,[2 1 3]);
%cf_1km = permute(Cloud_Fraction_Liquid.timeseries3,[2 1 3]);
cf_1km = ['Need to add this for CF>0.0!'];
file_save = [file_dir2 'Nd_re_tau_cf_21_' num2str(time(1,1)) '_SZA_LT_65_CF_GT_' CF_str2 '_CTH_LT_3.2km' seaice_str seaice_max_2week_str '.mat'];

save(file_save,'-V7.3','Nd','tau_1km','re_1km','cf_1km','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
mat2nc_Dan(file_save,[file_save '.nc']);

%3.7um
Nd = permute(N_time3_37,[2 1 3]);
Nd_1km = permute(Droplet_Number_Concentration_37.timeseries3,[2 1 3]);
tau_1km = permute(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,[2 1 3]);
re_1km = permute(Cloud_Effective_Radius_37_Liquid_Mean.timeseries3,[2 1 3]);
%cf_1km = permute(Cloud_Fraction_Liquid.timeseries3,[2 1 3]);
file_save = [file_dir2 'Nd_re_tau_cf_37_' num2str(time(1,1)) '_SZA_LT_65_CF_GT_' CF_str2 '_CTH_LT_3.2km' seaice_str seaice_max_2week_str '.mat'];

save(file_save,'-V7.3','Nd','tau_1km','re_1km','cf_1km','time','lon','lat','Nd_1km','Notes','make_mockL3_variables_script_name','multiL2L3_case','multiL2L3_project','daily_averaged_files_loc2');
mat2nc_Dan(file_save,[file_save '.nc']);


end