% 10th Jan, 2017 - Nd data from mock L3 data for day 213 of 2005 to day 121
% of 2012 (7 full SH summers, but using global data)
% First of all run MODIS_multi_DAY_processL3L2.m to load in the mock L3
% data. Settings :-
%  multiL2L3_case = 'load L3 and concatenate';
%  multiL2L3_project='global';
%  direcs={'aqua'};  years=[2005:2012]; days={[213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121]};
%  daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/'];

iload_seaice=1;
isave=1;
%savedir = '/home/disk/eos1/d.grosvenor/saved_misc_mat_files/';
%savedir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/';
savedir = '/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_1_3.2km_SZA_65_meanTau10/';
save_name = 'Nd_saved_sea-ice-edge_Jan2017';

%default, changed below
mask_size=0; tag='_no_landmask';

%Can just set mask_size to =1 for no smoothing
mask_size=1;   tag='_landmask_no_smooth'; %Can just set mask_size to =1 for no smoothing
%mask_size=2;   tag='_landmask_smooth02'; %Can just set mask_size to =1 for no smoothing
%mask_size=4;   tag='_landmask_smooth04';

months=[9:12 1:4];
Nmonths=length(months);



date_str = datestr(now,30);
savefile01 = [savedir save_name tag '_' date_str '.mat'];

dataset_str = 'SH_summer_2005_to_2011';
var_tags = {'Nd_16','Nd_21','Nd_37'};
Date_text='SHemisphere summer (SONDJF) 2005 to 2011';
Notes='Monthly mean Nd from mock L3 using 1km Nd values.';

iscreen_seaice=1; %
iscreen_land=1;
iscreen_region=0; %Whether to apply the regional screening (same for all years, but different for each month - i.e. all Jans use the same one)
iscreen_sza=0; % This refers to L3 screening (i.e. based on daily averaged data)
%max_sza_thresh = 65;

iscreen_ndays=0;
%Exclude points below a threshold number of days over the whole season
thresh_ndays=5;

if iload_seaice==1
    %% Load and combine NH and SH seaice (is only SH in the file)
    %seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2007_all_NHemisphere_20150529T103046.mat';
    %seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2008_JJA_all_NHemisphere_20161222T222241.mat';
    %seaiceNH = load(seaice_fileload_1deg_NH);
    %seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2008_JJA_all_SHemisphere_20161221T095015.mat';
    %seaiceSH = load(seaice_fileload_1deg_SH);
    %seaice_fileload_1deg_NH = '/home/disk/eos8/d.grosvenor/sea_ice_data/north/saved_seaice_1deg_2000-2015_all_NHemisphere_20161223T184753.mat';
    %seaiceNH = load(seaice_fileload_1deg_NH);
    seaice_fileload_1deg_SH = '/home/disk/eos8/d.grosvenor/sea_ice_data/south/saved_seaice_1deg_2000-2015_all_SHemisphere_20170110T094057.mat';
    seaiceSH = load(seaice_fileload_1deg_SH);
    
    single_hemisphere='single'
    switch single_hemisphere
        case 'single'
            [seaice_time3,seaice_max_time3] = seaice_match_times_FUNC(seaiceSH,Cloud_Fraction_Liquid,modisyear_timeseries3,daynum_timeseries3);
        case 'combine'
            %Script to combine the two together - not needed for SO case
            sea_ice_combine_NH_SH
    end
    
    
    
end



clear idays Nd_16 Nd_21 Nd_37 sea_ice_max_tmp sea_ice_tmp


Nd_16 = N_time3_16;
Nd_21 = N_time3;
Nd_37 = N_time3_37;

siz=size(Nd_16);

if iscreen_seaice==1    
    
    %find points where there is any sea-ice using the threshold below
    sea_ice_thresh = 1e-5;
    %    sea_ice_thresh = 1e9;    
    isea_high = find(seaice_max_time3>sea_ice_thresh);
    %    isea_high = find(sea_ice_max>sea_ice_thresh); %max over 2 week window
    Nd_16(isea_high)=NaN;
    Nd_21(isea_high)=NaN;
    Nd_37(isea_high)=NaN;
    
    Notes = [Notes ' Seaice screening with a threshold of ' num2str(sea_ice_thresh) '.'];
    
end

if iscreen_land==1  %This needs changing to work for mock L3
    
    landmask_load= load('/home/disk/eos1/d.grosvenor/amsre_land_mask.mat');
    
    %landmask = repmat(flipdim(landmask_load.amsre_land_mask,1),[1 1 Nmonths]);
    %landmask = flipdim(landmask_load.amsre_land_mask,1);

    %Smooth out the land mask to make it more aggressive - seemed to be
    %allowing some small bits of land otherwise.
    %Also need to do this for the sea ice
    halo=8; 
    fs=fspecial('average',mask_size); %Make the averaging filter
    landmask_tmp = add_halo(landmask_load.amsre_land_mask,halo);
    landmask_smooth = filter2(fs,landmask_tmp,'same');
    landmask_smooth = remove_halo(landmask_smooth,halo);
    %landmask_smooth = repmat(flipdim(landmask_smooth,1),[1 1 Nmonths]);
    landmask_smooth = flipdim(landmask_smooth,1);
    
    
    
    
    %Land mask is zero where there is land and NaN where not - so can just add to the
    %array to make the non-land regions NaN
    
    %Make the land regions NaN
    %    Nd_16_ocean = Nd_16 + landmask2;
    %    Nd_21_ocean= Nd_21 + landmask2;
    %    Nd_37_ocean = Nd_37 + landmask2;
    %    std_Nd_16_ocean = std_Nd_16 + landmask2;
    %    std_Nd_21_ocean= std_Nd_21 + landmask2;
    %    std_Nd_37_ocean = std_Nd_37 + landmask2;
    
%     eval_str = ['Nd' var_tag ' = Nd' var_tag ' + landmask2;']; eval(eval_str);
%     eval_str = ['Nd_Ndays' var_tag ' = Nd_Ndays' var_tag ' + landmask2;']; eval(eval_str);
%     eval_str = ['Nd_std_dev' var_tag ' = Nd_std_dev' var_tag ' + landmask2;']; eval(eval_str);
%     

    
end






%% Make an overall mean for all days of the season

for ivar=1:length(var_tags)
    var_tag = var_tags{ivar};
    
    for im=1:length(months)
        month_no = months(im);
        %days of month funciton will be useful to give days required for each month
        %- e.g. 32:60 for Feb.        
        [days_required_for_mean,time_mean_str] = days_of_month(month_no);
        
        date_tag = [dataset_str '_' time_mean_str];
        
        time_inds=[];
        for it=1:length(days_required_for_mean)
            ifind = find(daynum_timeseries3_MODIS==days_required_for_mean(it));
            if length(ifind)>0
                time_inds = cat(2,time_inds,ifind);
            end
        end
        [temp,temp_Ndays,temp_std_dev] = eval(['meanNoNan(' var_tag '(:,:,time_inds),3);']);
        
        if iscreen_ndays==1
            %Exclude data where there are fewer than thresh_nyears in the
            %climatology
            inan = find(temp_Ndays<thresh_ndays);
            temp(inan) = NaN;
            temp_std_dev(inan) = NaN;            
        end
        
        if iscreen_land==1                       
            
            %Land mask is zero where there is land and NaN where not - so can just add to the
            %array to make the non-land regions NaN
            
            temp = temp + landmask_smooth;
            temp_std_dev = temp_std_dev + landmask_smooth;
            
            Notes = [Notes ' Landmask screening with mask_size = ' num2str(mask_size) '.'];
            
            
        end
        
        
        eval_str = [var_tag '_' date_tag ' = temp;']; eval(eval_str);
        eval_str = [var_tag '_' date_tag '_Ndays = temp_Ndays;']; eval(eval_str);
        eval_str = [var_tag '_' date_tag '_std_dev = temp_std_dev;']; eval(eval_str);
        
        
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
    
end

if isave==1
    %Convert to netCDF
    mat2nc_Dan_matlab_netcdf(filename_savevars,[filename_savevars '.nc']);
end
