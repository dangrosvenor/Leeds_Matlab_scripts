function [model_match,modis_match,d_match] = Hawaii_aod_satellite_model_match_func02(modis_aod,model_lat2D,model_lon2D,model_dat,time_out_var)

%% Load the model and modis data from the .mat files
%UM_base_dir = '/home/disk/eos15/d.grosvenor/UM/Hawaii/';
%um_case_PD = 'u-co295';  %volcano ON    

% Load the model AOD data from the .mat file
%aod_save_file = [UM_base_dir um_case_PD '/aod550_total.mat'];
%model_aod = load(aod_save_file);

% Load the MODIS AOD data
%filedir='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/Hawaii_Dec2020_ADVANCE/'    
%file_name_h5='MYD04_L2.A2020356.0055.061.2020357191130.hdf';
%save_name = 'Hawaii_all_L2_04_data.mat';
%modis_aod = load([filedir save_name]);
%modis_aod = load([modis_aod_file]);

%% Loop over the MODIS swaths (Aqua and Terra) and match to the nearest model time
% Just do this for the region of interest (Bottom left corner = 11.9674,
% -180, top right = 22.4980, -152.1999; these were the ones used for the
% MODIS AOD request for swaths).

lons_roi = [-180 -152.1999];
lats_roi = [11.9674 22.4980];
dt_tol = 31; %tolerance for time difference between MODIS and model in minutes.
dlat_init = 20/111; %smaller area to search for in model grid - will be the modis pixel +/- this number of degrees
dlon_init = 20/111; %smaller area to search for in model grid - will be the modis pixel +/- this number of degrees

nswaths = length(modis_aod.AOD_550_Dark_Target_Deep_Blue_Combined_all_times);
%loop over swaths
for isw=1:nswaths
    isw
    modis_time = minALL(modis_aod.scantime_matlab_all_times{isw}); %scantime gives one value per pixel - just take the min since the swath only
        %takes 5 minutes.
    %Find the nearest model output (in UTC)
    [min_val,it] = min(abs(modis_time - time_out_var)); %times are in days (since 1-Jan-0000)
    if min_val*24*60 < dt_tol
        igood = find( modis_aod.lat_all_times{isw} >= lats_roi(1) & modis_aod.lat_all_times{isw} < lats_roi(2) ...
            & modis_aod.lon_all_times{isw} >= lons_roi(1) & modis_aod.lon_all_times{isw} < lons_roi(2) ...
            & modis_aod.AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag_all_times{isw}==3);
        
        
        %Match each pixel to a model pixel.
        for ip=1:length(igood)
            i = igood(ip);
            lat_modis = modis_aod.lat_all_times{isw}(i);
            lon_modis = modis_aod.lon_all_times{isw}(i);
            %Narrow down the search for the nearest model point.
            lat_modis_init = [lat_modis - dlat_init lat_modis + dlat_init];
            lon_modis_init = [lon_modis - dlon_init lon_modis + dlon_init];
            i_init = find( model_lat2D >= lat_modis_init(1) & model_lat2D < lat_modis_init(2) ...
                & model_lon2D >= lon_modis_init(1) & model_lon2D < lon_modis_init(2) );
            
            %Find the nearest pixel and extract the model AOD
            if length(i_init)>0
                latA = repmat(lat_modis,[length(i_init) 1]);
                lonA = repmat(lon_modis,[length(i_init) 1]);
                latB = model_lat2D(i_init);
                lonB = model_lon2D(i_init);
                [d,aob]=distlatlon(latA,lonA,latB,lonB);
                
                [minval,i_match] = min(d);
                
                dat_it = model_dat(:,:,it);
                
                model_match{isw}(ip) = dat_it(i_init(i_match));
                modis_match{isw}(ip) = modis_aod.AOD_550_Dark_Target_Deep_Blue_Combined_all_times{isw}(i);
                d_match{isw}(ip) = minval;
            end
        end
        
    end
    
end
