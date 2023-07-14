lat_lon_inds_str = [':,:'];

try
    irun_from_modis_multi_DAY_processL3L2=1;

%  ** MODIS_process_global_L2_daily_with_screening ** -  this is generally run,
%  except for process action when uses this:-
%  ** modis_multi_dir_process **     

% for process action - make_MODIS_variables_01_mockL3  is where
% the data is added to the large array and
% filedir = the place where the L2 joint files are stored
% this is now set below.
% processed_daily_swath_loc is where the 20 swath processed data is stored 
    
%% JUST for 'load L3 and concatenate' and 'load, average and save' NOT process action
% However, for 'load, average and save' a subset of variables are first loaded for screening
%, which are are chosen in
%   MODIS_process_global_L2_daily_with_screening
% rather than using these files
  
% The process action uses make_mockL3_variables variables
% Choose the desired script below, or use make_mockL3_variables_Arctic_SZA_study - this contains the
%       load_special2 case. May get an error of undefined modis_var if have
%       selected load_special2 with the wrong
%       make_mockL3_variables_script_name
% N.B. - don't put .m at the end of the script names

make_mockL3_variables_script_name = 'make_mockL3_variables';
make_mockL3_variables_script_name = 'make_mockL3_variables_Arctic_SZA_study'; %inc. Cloud_Fraction_Liquid2, etc.
make_mockL3_variables_script_name = 'make_mockL3_variables_Arctic_SZA_study_without_CTT_ice_liq_etc'; %For VOCALS 2007 Jan 2015 tests
%make_mockL3_variables_script_name = 'make_mockL3_variables_DanMcCoy';
%make_mockL3_variables_script_name = 'make_mockL3_variables_just_for_Nd_May2017';
              % *** see make_mockL3_variables_DanMcCoy_saved for the vars
            % saved for Daniel. ***
%make_mockL3_variables_script_name = 'make_mockL3_variables_mockL3_basic'; 
    %Basic variables for the earlier mockL3 files
    
%make_mockL3_variables_script_name = 'make_mockL3_variables_select'; %e.g. just the date for e.g. POLDER comparison
%make_mockL3_variables_script_name = 'make_mockL3_variables_Puijo_CDP'; %e.g. just the date

% - N.B. the save_or_load variable determines some extra vars to be loaded
% this is set in :-
%   MODIS_process_global_L2_daily_with_screening.m


post_process_L2_and_L3 = 1; %flag to determine whether to post process the data or not after 
%'load L3 and concatenate' and 'load processed L2 and concatenate'. e.g.
%create sst array, run filtering_data_get etc. - may not want to if have only
%loaded a few variables

    iload_amsre_mat = 1; %flag to tell it load amsre data for the time range - saved in a mat file
       %N.B. - not used for the process action.
       
    ignore_amsre = 0; %Flag to say to not use amsre data for load, average* action.   
       
       %to switch off sst matching set ino_sst_matching=1 in the command window before running
                    %N.B. will need to do this every time as it is reset
                    %upon running
                    
                    
%% -- Default flags for the process action (will be overwritten below) - whether to apply the water path
% confidence screening in open_L2_MODIS_file_01 and whether to limit re16
% and re37 to <20um as done for re21 when the water path confidence screening is
% done
    i_apply_water_path_confidence = 1;                    
    i_limit_16_37_to_20_um = 0;  %N.B. will only want to apply this when doing water path screening
    i_make_wavelengths_consistent = 0; %Makes pixels NaN for all wavelengths if any wavelength has NaN value
        % So, all wavelengths will be using the same pixels.
%The above are the defaults that may get overwritten below for the
%regional switch (multiL2L3_project switch)

                    
%% choose the options for what action to do
    
multiL2L3_case = 'process';  %process L2 joint files into the 20 swath arrays (one for each day)
  %at the moment this only considers liquid pixels, so ice info is lost
  %(except for ice cloud fraction). Saves .mat files in
  %processed_daily_swath_loc = ['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];
  %One for each day, all in the same dir. Only the most recent ones are
  %used for the other multiL2L3_case actions (below). Are
  %separate directories for aqua (if it gets included) and terra.
  % Note - joint L2 swaths with no daylight in them (no SZA<81.4) are not
  % processed by MODIS_process_multiple_L2_files and so will not be
  % included in the big arrays.

%multiL2L3_case = 'load, average and save';  %loads in a series of the daily files produced by the process action
% and does screening to produce average values for each day. These are saved in
  %directories that indicate the screening that was done.
  %N.B. - during a first pass variables are loaded just for screening, which are chosen in
  % MODIS_process_global_L2_daily_with_screening
  % rather than using the make_mockL3_variables_script_name file. Then
  % later the latter file is used for the main variables.
  
%multiL2L3_case = 'load L3 and concatenate'; %stiches together several of
  %the above daily average (mock L3) files to make a large array ready for plotting
  
%multiL2L3_case = 'save concatenated mockL3 data'; %save the above files to a .mat file
  

%multiL2L3_case = 'load processed L2 and concatenate'; %loads in processed L2 files and
  %cats them together without doing daily averaging (e.g. as done previously
  %for the SZA work - allows multiple screenings on the fly. But uses lots of memory
  %for multiple years. Might be ok if we restrict the region.

  
%% Choose the region to work on  
  
%need to set this properly too in order to choose the correct folder for
%the files
     multiL2L3_project='Arctic_summer_JointL2';
     multiL2L3_project='Arctic_summer_JointL2_box_only';  %have more years for the box only region  
     multiL2L3_project='global';
%     multiL2L3_project='Southern_Ocean_box_only'; %for sub-sets of the global data set the global_region below 
%     multiL2L3_project='Southern_Ocean_all_lats'; %N.B. - this should be all lons, not all lats!
%    multiL2L3_project='Southern_Ocean_all_lons_no_water_path_confidence';
%    multiL2L3_project='Southern_Ocean_all_lons_water_path_confidence_ON';    
%    multiL2L3_project='Southern_Ocean_RITS94_GASSP_region';    
    
     %for sub-sets of the global data set the global_region below 
     %along with corresponding LAT_val and LON_val values
%     multiL2L3_project='Tasmania (Boers)';
%     multiL2L3_project='VOCALS';
   % --- N.B. set MODIS_filetype_choose in later on when processing to either L2 or
   % jointL2 to say which one we want to look at
%      multiL2L3_project='Puijo_Sami'; %CDP matches
%      multiL2L3_project='Puijo_Sami_L3'; %L3 CDP matches      
%      multiL2L3_project='Puijo_Sami_jointL2';      
%      multiL2L3_project='Puijo_Sami_DMPS_L2';
%     multiL2L3_project='VOCALS_JointL2_box_only';  %have processed aqua 2007 for the box region only
     %(whereas have terra 2007 globally)
     
%     multiL2L3_project = 'VOCALS_JointL2_box_only_no_confidence_screen'; %NOT applying water path confidence =very good (allows marginal, good and verygood, but not bad)
%     multiL2L3_project = 'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37' %Limiting re16 and re37 to 20um as happens for re2.1 (with "very good" water path confidence screening)
%     multiL2L3_project = 'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37_consistent_pixels' %Limiting re16 and re37 to 20um as happens for re2.1 (with "very good" water path confidence screening)
%     multiL2L3_project = 'Iceland_30W30E_40N90N'; Npix_Iceland = 10;
%                                                    %think I set this as 0.25deg lat and lon instead
     
     %Set the directory of where to store and load mock L3 (can be global
     %or regional)
%      switch multiL2L3_project
%          case 'global'
%              dir_for_mockL3 = '';
%          otherwise
%              dir_for_mockL3 = multiL2L3_project;
%      end

     switch multiL2L3_project
         case {'Southern_Ocean_box_only','Tasmania (Boers)','VOCALS','Southern_Ocean_all_lats'}
             %this is if we have selected the global option above, but only want to load
             %a smaller region
             global_region=multiL2L3_project;
             multiL2L3_project='global'; %need to set this back to global to give the right directory
         case {'Puijo_Sami_L3'}
             %this is if we have selected the global option above, but only want to load
             %a smaller region
             global_region=multiL2L3_project;
             multiL2L3_project='global_L3'; %need to set this back to global to give the right directory
             
         otherwise
             global_region='global';
     end
     
     

     
% ***  only needed for 'process' action  ***





ialternate_file_structure=0; %set this flag to one for the Puijo station MODIS data where all the files are
%located for a whole year, rather than day by day.

%whether to only open files that match certain times
time_restrict=0;  dtime_restrict=1; %in hours: can be +/- this many hours

MODIS_filetype_choose = 'L2 Joint';
                        
                        
clear MLAT MLON

switch multiL2L3_case

%% Choose the screening - for everything except process action
%   (switch multiL2L3_case) 
    case {'load, average and save','load L3 and concatenate','load processed L2 and concatenate','save concatenated mockL3 data'}
        %this directory is where to save to for 'load, average and save' and where
        %to load from for *** 'load L3 and concatenate' *** - select the required
        %screening (so, this dir does not need to be set for 'load processed L2 and concatenate', but enters here because it sets the 
        % lat and lon of the sub-region, if using)
        
        % Based upon the screening case requested below, the acutal screening to do is chosen in
        %    MODIS_process_global_L2_daily_with_screening   
        %   Add a case in the above for the regional directory that you
        %   are dealing with.
        
        % N.B. - have done L3 data for the VOCALS region for Aqua for one
        % year, but don't have global - so have created additional entries
        % for 'no screening' case in a special folder
        
    %no screening case:-
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_minCTT_173_ice_allowed_SZA_0-90_thresh_NP_0/']; %i.e. no screening        
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/no_screening/']; %no screening           
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37/no_screening/']; %no screening VOCALS tests (re<20 for 1.6 and 3.7um)                  
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/VOCALS_JointL2_box_only_no_confidence_screen/no_screening/']; %no screening VOCALS tests with no water path confidence screening                          
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37_consistent_pixels/no_screening/']; %no screening VOCALS tests (re<20 for 1.6 and 3.7um and making the pixels consistent - processed for 2005 to end 2011)     
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Southern_Ocean_all_lons_no_water_path_confidence/CF_0.8_minCTT_173_ice_allowed_SZA_65/']; %

    %screening cases:-
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_273_SZA_65/'];
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_268_SZA_65/']; 
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_268_ice_allowed_SZA_65/'];           
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.6_meanCTT_273_ice_allowed_SZA_65/'];                   
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.6_meanCTT_268_ice_allowed_SZA_65/'];        
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_273_ice_allowed_SZA_65/'];          
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_268_ice_allowed_SZA_65/'];        
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_278_ice_allowed_SZA_65/'];   
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_280_ice_allowed_SZA_65/'];           
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_284_ice_allowed_SZA_65/'];                   
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_273_ice_allowed_SZA_65/'];                           
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65_tau.GT.5/'];    
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_268_ice_allowed_SZA_65_tau.GT.5/'];          
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65_tau.LT.5/'];          
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_0-90_tau.GT.5/'];         
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_SZA_75-85/'];         
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65-85/'];                 
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_0-90/']; 

%% most used global ones I think
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_65/'];  %Old data given to Daniel McCoy
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_minCTT_173_ice_allowed_SZA_65/'];          
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65_meanTau10/']; %One used for sea-ice edge stuff, Jan 2017 (added mean Tau screening)
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_1_3.2km_SZA_65_meanTau10/']; %One used for sea-ice edge stuff, Jan 2017 (added mean Tau screening)


%Used for Dan McCoy 2003-2015 dataset and Matt C for JJA 2008 comparison in Nd review paper :-
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_173_meanCTH_3.2km_SZA_65/'];  
        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_meanCTT_173_meanCTH_3.2km_SZA_65/'];          

        % Added for ACSIS metric to calculate the fraction of days with Sc (low cloud
        % >80%). Needed for the denominator of the fraction in order for
        % clear skies (no cloud) to be counted
        %daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.0_SZA_65/'];

%% ...end of most used ones
% Others :-
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_55/']; 
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_173_ice_allowed_SZA_45/']; 
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_263_no_iceCF_SZA_65/'];                        
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_minCTT_268_meanCTH_3km_SZA_65/'];   

% -- 2014 SZA paper ---
        %processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.6_CTT_273_SZA_65/'];
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.99_meanCTT_268_ice_allowed_SZA_65/'];   %re-done -- think have this globally. At least do the Southern Ocean
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.99_meanCTT_268_ice_allowed_SZA_0-90/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.99_meanCTT_173_ice_allowed_SZA_0-90/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.8_meanCTT_173_ice_allowed_SZA_0-90/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.8_meanCTT_173_ice_allowed_SZA_0-90_minfrac_CF_0.5/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.8_meanCTT_173_ice_allowed_SZA_0-90_minfrac_NpNd_0.5/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.8_meanCTT_173_ice_allowed_SZA_0-90_minfrac_CF_0.5_minfrac_NpNd_0.5/']; %re-done
%         daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/SZA_paper_screening_CF_0.1_meanCTT_173_ice_allowed_SZA_0-90_minfrac_CF_0.0_minfrac_NpNd_0.0/']; %re-done

   %8th April, 2014 - screening for Southern Ocean data
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/CF_0.8_meanCTT_263_meanCTH_3.2km_SZA_65/'];  

% Iceland 2014 with Anja
%           daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Iceland_30W30E_40N90N/CF_0.8_minCTT_173_ice_allowed_SZA_65/'];  
%           daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Iceland_30W30E_40N90N/CF_0.0_minCTT_173_ice_allowed_SZA_65/'];          

%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Iceland_30W30E_40N90N/CF_0.8_meanCTT_173_meanCTH_0.5_3.2km_SZA_65_meanTau10/']; %One used for sea-ice edge stuff, Jan 2017 (added mean Tau screening)       
%        daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L2/averaged_daily_mockL3/Iceland_30W30E_40N90N/CF_0.0_meanCTT_173_meanCTH_0.5_3.2km_SZA_65_meanTau10/']; %One used for sea-ice edge stuff, Jan 2017 (added mean Tau screening)       
        
%default value - leave like this
        lat_lon_inds_str = ':,:'; %a  string like '1:12,20:30' to select
        
        switch multiL2L3_project
            case {'global','global_L3'}
                
                switch global_region  % --- N.B. - this does not apply to the process action ---
                    case 'global'
                        LAT_val = []; LON_val = []; %select this for all of the region
                        %or specify LAT_val and LON_val for a portion of the globe
                    case 'Southern_Ocean_box_only'
                        LAT_val = [-60 -45]; LON_val = [50 100]; %Southern Ocean
%                        LAT_val = [-70 -60]; LON_val = [-140 -60]; %Southern Ocean - looking further south
                          %at region where get less sea-ice
                    case {'Southern_Ocean_all_lats'}
                        LAT_val = [-60 -45]; LON_val = [-180 180]; %Southern Ocean  
                        LAT_val = [-60 -40]; LON_val = [-180 180]; %Southern Ocean   
                        LAT_val = [-60 -40]; LON_val = [-180 180]; %Southern Ocean - did this for Daniel McCoy Science paper 2014
                        LAT_val = [-90 -20]; LON_val = [-180 180]; %Southern Ocean - trying a larger region for 2015 analysis                        
%                        LAT_val = [-45 -40]; LON_val = [-180 180]; %Southern Ocean                           
                            %N.B. searches on MLON_AMSRE and MLON
                            %MLON(end)=179.5 and searches for
                            %<LON_val(end), so make sure LON_val(end) is
                            %above 179.5 if want all points included.
                            
                    case 'Southern_Ocean_all_lons_no_water_path_confidence'
                        %Test to see if can load a smaller region of the
                        %20-70S, all lons data to save memory etc. - for
                        %SEASCAPE proposal CTH histograms.
                        % Region shown in box in proposal - not sure if the
                        % aircraft can cover all of this or not?
                         LAT_val = [-70 -60]; LON_val = [-20 0];
                            
                    case 'Tasmania (Boers)'
                        LAT_val = [-50 -40]; LON_val = [120 150]; %Tasmania, Boers, OJRMS, 1998 - Southern Ocean
                        LAT_val = [-60 -40]; LON_val = [40 144]; %Tasmania, Boers, OJRMS, 1998 - Southern Ocean                        
                    case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2','Puijo_Sami_L3'}  %at 62.91N, 27.66 E
                        LAT_val = [62 63]; LON_val = [27 28]; %Puijo_Sami : 62.91N, 27.66 E
%                        LAT_val = [62 64]; LON_val = [26 29]; %Puijo_Sami : 62.91N, 27.66 E   
                        
                       
                                              
                    case {'VOCALS','VOCALS_JointL2_box_only'}
%                        lats_multi = {[-21.5 -18.5],[-21.5 -18.5],[-21.5 -18.5]};
%                        lons_multi = {[-76.25 -71.25], [-81.25 -76.25], [-87.25 -81.25]};
%                        LAT_val = [-30.5 -19.5]; LON_val = [-73.25 -72.25]; %VOCALS longitude bin, one lat
                        LAT_val = [-30.5 -19.5]; LON_val = [-77.25 -76.25]; %VOCALS longitude bin, one lat                        
%                        LAT_val = [-60 -40]; LON_val = [40 144]; %                                               
                        LAT_val = [-40 0]; LON_val = [-140 -68]; %
                        
                    case {'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37','VOCALS_JointL2_box_only_no_confidence_screen'}
                        LAT_val=[-50 20]; LON_val=[-160 -60]; %VOCALS limiting to the box
                        
                   

                end  %switch global_region
                
          %Cases where the processed files are for a limited region only - N.B. - this not required for the process action  
          %Just for the load_average_save, and load L3 and concatenate actions (I think)
            case {'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37','VOCALS_JointL2_box_only_no_confidence_screen','VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37_consistent_pixels'}
                LAT_val=[-50 20]; LON_val=[-160 -60]; %VOCALS limiting to the box
            case {'Arctic_summer_JointL2_box','Arctic_summer_JointL2_box_only'}
                LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer
            case {'Southern_Ocean_box_only'}
                LAT_val = [-60 -45]; LON_val = [50:100]; %Southern Ocean 
            case {'Southern_Ocean_RITS94_GASSP_region'}
                LAT_val = [-80 -50]; LON_val = [-140:-60]; %Southern Ocean                 
            case {'VOCALS_JointL2_box_only'}
                LAT_val = [-40 0]; LON_val = [-140 -68]; %
            case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2'}
                LAT_val = []; LON_val = []; %    
            case {'Southern_Ocean_all_lons_no_water_path_confidence','Southern_Ocean_all_lons_water_path_confidence_ON'}
                LAT_val = [-90 -20]; LON_val = [-180 180]; %Southern Ocean - trying a larger region for 2015 analysis
            case {'Iceland_30W30E_40N90N'}    %N.B. - this not required just for the process action
                LAT_val = [40 90]; LON_val = [-30 30]; %               
        end

%% Choose dates to do daily averaging for or to load for L3 and concatenate
%   Does not apply to process action
        switch multiL2L3_project  % NOTE - many regions have multiL2L3_project changed to 'global' above by this stage
            case {'global','Southern_Ocean_box_only'}
                
              if exist('i_override_MODIS_multi')~=1 %(==1 for a variable exising in workspace)
                %        direcs={'terra'};
                %        direcs={'aqua'};

                %        years=[2007:2010];
                                
                
%                direcs={'terra'}; years=[2006 2007]; days = {[306:365],[1:19 22:210]};                 
%                direcs={'terra'}; years=[2006 2007]; days = {[306:365],[1:19 22:243]}; %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333days_2006_2007.mat'; %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333_2007.mat';

%                direcs={'terra'}; years=[2007]; days = {[211:211]}; amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333days_2006_2007.mat'; %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333_2007.mat';

%this is the latest one:-
                direcs={'terra'}; years=[2006 2007]; days = {[335:365],[1:19 22:365]}; 
                direcs={'terra'}; years=[2006 2007]; days = {[335:365],[1:19 22:365]}; 
                direcs={'terra'}; years=[2006 2007]; days = {[335:365],[1:19 22:334]}; 


                                
%                direcs={'aqua'}; years=[2007]; days = {[1:365]};                 

                %direcs={'terra'}; years=[2007]; days = {[333:365]};                 
%                direcs={'terra'}; years=[2006 2007]; days = {[335:365],[1:19 22:243]};                 
                
%                direcs={'terra'}; years=[2007]; days = {[1:19 22:365]}; 
                               
                
%                direcs={'terra'}; years=[2007]; days = {[1:19 22:182]};
%                direcs={'terra'}; years=[2007]; days = {[182:365]};
                
%                direcs={'aqua'}; years=[2007]; days = {[354]};  %20th Dec   
%                direcs={'terra'}; years=[2007]; days = {[354]}; %20th Dec
%                direcs={'aqua'}; years=[2007]; days = {[171]}; %20th June                                              
%                direcs={'terra'}; years=[2007]; days = {[171]}; %20th June                              

%                direcs={'terra'}; years=[2007]; days = {[1:19 22:25]};     
                
%                direcs={'terra'}; years=[2007]; days = {[333:365]};
                
                %shorter period for test
%                direcs={'terra'}; years=[2006]; days = {[335]};  
                
                
                
%                direcs={'terra'}; years=[2007]; days = {[211:244]};                      
%                direcs={'terra'}; years=[2006 2007]; days = {[335:365],[1:19 22:365]};                      
                %day 335 of 2007 (&2006) is 1st Dec

%                direcs={'terra'};  years=[2007]; days={[152:243]}  %JJA 2007 for Koren region study                
%                direcs={'aqua'};  years=[2006 2007]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)
                %For processing from 20swath to daily :-
                direcs={'aqua'};  years=[2005:2012]; days={[213:365],[1:121 213:365],[1:365],[1:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121]};  %processing just when have SZA<65 data (01Aug06 - 30Apr07)
                direcs={'aqua'};  years=[2008:2012]; days={[122:212 317:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365]};  %                
                %Matt C JJA 2008  (Jan 2017)
                direcs={'aqua'};  years=[2008]; days={[153:244]};  %JJA for 2008 for Matt C
%                direcs={'aqua'};  years=[2008]; days={[1:365]};  % Whole of  2008
                %direcs={'aqua'};  years=[2008]; days={[1:152 245:366]};  % 
                direcs={'aqua'};  years=[2006]; days={[1:365]};  % Whole of  2006
                %Sea-ice edge Jan 2017
                %direcs={'aqua'};  years=[2005:2012]; days={[213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121]};
                
                %Missed out day 122 of 2003 somehow
                direcs={'aqua'};  years=[2003:2014]; days={[123:365],[1:366],[1:365],[1:365],[1:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365]}; %,[1:365]};%    
                
               
                %Had to change variable script since *** not available for
                %day 131 of 2007
direcs={'aqua'};  years=[2007:2014]; days={[131:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365]}; %,[1:365]};%                                                

                %Had to take out Cloud_Top_Temperature_Day_Standard_Deviation.timeseries3 variable script since not available for
                %day 211 of 2007. And a few other variables too.
direcs={'aqua'};  years=[2007:2014]; days={[211:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365]}; %,[1:365]};%                                                
                %no data for day 366 of 2008
direcs={'aqua'};  years=[2009:2015]; days={[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365],[1:365]}; %,[1:365]};%                                                

direcs={'aqua'};  years=[2008]; days={[366]};

% For loading
direcs={'aqua'};  years=[2003]; days={[1:365]};
 direcs={'aqua'};  years=[2004]; days={[1:366]};
%  direcs={'aqua'};  years=[2005]; days={[1:365]};
%  direcs={'aqua'};  years=[2006]; days={[1:365]};
  direcs={'aqua'};  years=[2007]; days={[1:365]};
%  direcs={'aqua'};  years=[2008]; days={[1:366]};
%  direcs={'aqua'};  years=[2009]; days={[1:365]};
  
% direcs={'aqua'};  years=[2010]; days={[1:166 168:365]};
 %direcs={'aqua'};  years=[2011]; days={[1:365]};
 %direcs={'aqua'};  years=[2012]; days={[1:366]};
 %direcs={'aqua'};  years=[2013]; days={[1:365]};
 %direcs={'aqua'};  years=[2014]; days={[1:365]};
 %direcs={'aqua'};  years=[2015]; days={[1:365]};
 

% direcs={'aqua'};  years=[2003:2015]; days={[1:365],[1:366],[1:365],[1:365],[1:365],[1:366],[1:365],[1:166 168:365],[1:365],[1:366],[1:365],[1:365],[1:365]}; %,[1:365]};%

 
%  direcs={'aqua'};  years=[2010]; days={[30]};
  direcs={'aqua'};  years=[2009 2010]; days={[87:365],[1:88]};  %ACSIS 2009-2010 global run
  %direcs={'aqua'};  years=[2014]; days={[242:253]};  %Iceland, Aug30-Sep10th, 2014
  
  direcs={'aqua'};  years=[2014]; days={[321]};  %Rosenfeld test case from ppt 17th Nov, 2014
  %direcs={'aqua'};  years=[2014]; days={[1:365]};  %Rosenfeld test case from ppt 17th Nov, 2014
  
                %(max extent of processed data so far)
                
                %        direcs={'terra'}; years=[2006 2007]; days = {[336:365],[182:210]}; %data that forgot to process
                %        direcs={'terra'}; years=[2007]; days = {[211:243]}; %new data
%                        direcs={'terra'}; years=[2006 2007]; days = {[306:365],[1:19 22:210]}; %max extent of processed data so far


                        
                %days determines the size of the pre-allocated array, so set it properly
                %        days = {[306:365],[1:59]}; %original
                %        days = {[336:365],[1:59]};
                %        days = {[22:59]};
                %        days = {[306:365],[1:19 22:32]}; %
                %        days = {[306:335]}; %
                %        days = {[34:59]}; %

                %        days = {[33:59]}; %days 20-21 seem to be missing

                %        days = {[306:335],[1:19 22:181]};
                
        end
                
                %amsre_ssts_global_2006-2010.mat contains 5 day smoothed
                %data from 2006-2010 (including Dec 2005 and Jan 2007 to
                %avoid smoothing problems at end and beginning)
%                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
                
%Trying out the HADSST dataset (interpolated from weekly data) that runs
%from 1990 to Oct 2015. This was made using sst_make_amsre_equivalent_from_HADSST.m
                %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/HADSST_matched_to_amsre_ssts_global_1990-Oct2015.mat';               
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/HADSST_matched_to_amsre_ssts_global_1990-Apr2017.mat';
                

                
             case {'global_L3'}
                 %still need to finish this read in method
                 
                direcs={'aqua','terra'}; years=[2006:2010]; days = {[1:366]}; amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333days_2006_2007.mat'; %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_333_2007.mat'; 

                
                %amsre_ssts_global_2006-2010.mat contains 5 day smoothed
                %data from 2006-2010 (including Dec 2005 and Jan 2007 to
                %avoid smoothing problems at end and beginning)
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
                
                
                daily_averaged_files_loc2=['/home/disk/eos8/d.grosvenor/saved_data_L3/'];
                


            case {'Arctic_summer_JointL2','Arctic_summer_JointL2_box_only'}
                %for the SZA box
                direcs={'terra','aqua'}; years=[2007:2010]; days = {[164:181]};  
%                direcs={'aqua'}; years=[2007:2010]; days = {[164:181]};                  
%                direcs={'aqua'}; years=[2010]; days = {[179]};                  
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_box_2007_2010.mat';
%                years=[2007]; days = {[164:165]};

            case {'VOCALS_JointL2_box_only'}
%                direcs={'aqua'}; years=[2007]; days = {[1:181 211:365]};                                 
                direcs={'aqua'}; years=[2008]; days = {[1:366]};   
                
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_box_2007_2010.mat';
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
                
            case {'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37','VOCALS_JointL2_box_only_no_confidence_screen','VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37_consistent_pixels'}
               %N.B. this does not appy to process action

                direcs={'aqua'}; years=[2005]; days = {[1:365]};                
                direcs={'aqua'}; years=[2006]; days = {[1:365]};  
                direcs={'aqua'}; years=[2007]; days = {[1:365]};  
                direcs={'aqua'}; years=[2008]; days = {[1:366]};                  
                direcs={'aqua'}; years=[2009]; days = {[1:365]};                   
                direcs={'aqua'}; years=[2010]; days = {[1:365]};                 
                direcs={'aqua'}; years=[2011]; days = {[1:365]};                  
                direcs={'aqua'}; years=[2012]; days = {[1:366]};                     
                
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_box_2007_2010.mat';
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_VOCALS_AQUAonly_2006_to_2010_20131025T044307.mat';

%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_raw_amsre_2002-2005.mat';
%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_Aqua_2002_to_2006_20150225T031906.mat';
                                
            case {'Southern_Ocean_all_lons_no_water_path_confidence','Southern_Ocean_all_lons_water_path_confidence_ON'}                        
                % Not for the process action
                
                
                direcs={'aqua'};  years=[2013 2014]; days={[213:365],[1:121]}  %
               
                %For loading L3
                direcs={'aqua'};  years=[2005:2012]; days={[213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
               
               amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
               
               
            case 'Southern_Ocean_RITS94_GASSP_region'   
                % Not for the process action                
                direcs={'aqua'};  years=[2013 2014]; days={[213:365],[1:121]}  %
               
                %For loading L3
                direcs={'aqua'};  years=[2005:2012]; days={[213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121 213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
               
%               amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
%Trying out the HADSST dataset (interpolated from weekly data) that runs
%from 1990 to Oct 2015. This was made using sst_make_amsre_equivalent_from_HADSST.m
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/HADSST_matched_to_amsre_ssts_global_1990-Oct2015.mat';               
               
            case {'Iceland_30W30E_40N90N'}
               %N.B. this does not appy to process action                 
                direcs={'aqua'}; years=[2014]; days = {[245:245]};   
                direcs={'aqua'}; years=[2014]; days={[242:253]};  %a
                
                ignore_amsre=0;
                
%                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_180_360_box_2007_2010.mat';
               

%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_VOCALS_AQUAonly_2006_to_2010_20131025T044307.mat';

%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_raw_amsre_2002-2005.mat';
%                amsre_matfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_Aqua_2002_to_2006_20150225T031906.mat';

                %amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/HADSST_matched_to_amsre_ssts_global_1990-Oct2015.mat';
                                
             
                                
            case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2'}   %load_L2 etc.
                direcs={'aqua','terra'}; years=[2006:2012]; days = {[1]};  
                direcs={'aqua','terra'}; years=[2006:2011]; days = {[1]};                  
%                direcs={'aqua','terra'}; years=[2006:2007]; days = {[1]};                  
%                direcs={'terra'}; years=[2006:2007]; days = {[1]};                    
                
%                iload_amsre_mat = 0; %flag to tell it load amsre data for the time range - saved in a mat file
                amsre_matfile = '/home/disk/eos8/d.grosvenor/saved_data_L2/amsre_ssts_global_2006-2010.mat';
                
                 %this determines which files are loaded - e.g. 11
                        %km boxes, 5 km, 2min time difference, or 15mins,
                        %etc.
                        % the script looks for files that have the
                        % following keywords in their filename
                        Puijo_keyword = 'L2_Puijo_matches_2013'; %30 seconds, 11x11km
                        Puijo_keyword = 'L2_Puijo_matches_2mins2013'; %2mins, 11x11km                        
%                        Puijo_keyword = 'L2_Puijo_matches_15mins2013'; %15mins 11x11km                                                                        
                        Puijo_keyword = 'L2_Puijo_matches_15mins_5km2013'; %15mins 5x5km
                        Puijo_keyword = 'L2_Puijo_matches_2mins_11km_assumeUTC2013'; %
                        Puijo_keyword = 'L2_Puijo_matches_2mins_11km_single_liquid_layers_only'; %with the number of pixels with single liquid layers
                        Puijo_keyword = 'L2_Puijo_matches_2mins_11km_UTC_minus_10mins'; %with the number of pixels with single liquid layers
                        Puijo_keyword = 'L2_Puijo_matches_2mins_5km_single_liquid_layers_only'; %with the number of pixels with single liquid layers                        
                        Puijo_keyword = 'L2_Puijo_matches_2mins_5km'; %with the number of pixels with single liquid layers                        
                        %The above was (accidentally) used for the data processed from the 0.5x1
                        %degree L2-joint data too
                        Puijo_keyword = 'L2_Puijo_matches_30mins_5km_DMPS';
%                        Puijo_keyword =
%                        'L2_Puijo_matches_30mins_11km_DMPS';
        end

        processed_daily_swath_loc=''; %for concatenation don't need the processed swaths (the 16 orbits)


       switch multiL2L3_case
           case {'load, average and save','load processed L2 and concatenate'}
               processed_daily_swath_loc = ['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];
       end

%% (switch multiL2L3_case - last case) For the process cases - choose lat,lon and the files to process       
    case 'process'

        switch multiL2L3_project
            case 'Arctic_summer_JointL2_box_only'
                time_restrict_case = 'choose here';
                override_mockL3_options2=1;
% LAT_val = [72 75]; LON_val = [-3:48]; %Arctic summer                
                LAT_min=71; LAT_max=76; LON_min=-4; LON_max=49; %Arctic summer files limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        
                years=[2007:2010];
%                years=[2010];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[164:181]};
 %              days = {[179:179]};
 
             case 'VOCALS_JointL2_box_only'
                 time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-40; LAT_max=0; LON_min=-140; LON_max=-68; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        
                years=[2007];
%                years=[2010];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[1:365]};
                
              case 'VOCALS_JointL2_box_only2'
                 time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-50; LAT_max=20; LON_min=-160; LON_max=-60; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        
                years=[2007];
%                years=[2010];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[1:365]};
                
            case 'VOCALS_JointL2_box_only_no_confidence_screen'
                %For the PROCESS action
                %whether to apply the water path
                % confidence screening in open_L2_MODIS_file_01 and whether to limit re16
                % and re37 to <20um as done for re21 when the water path confidence screening is
                % done
                i_apply_water_path_confidence = 0;
                i_limit_16_37_to_20_um = 0;  %N.B. will only want to apply this when doing water path screening
                
                filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-50; LAT_max=20; LON_min=-160; LON_max=-60; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        


                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly  
                
               years=[2008]; days={[1:351 355:361 363 365]}  %missing:- 352:354, 362, 364 and 366
               years=[2008]; days={[324:324]}  %all there for Aqua 2008 (1-366)    
               
               years=[2006]; days={[1:365]};  %all there for Aqua 2006 (1-365)
               years=[2005]; days={[1:365]};  %all there for Aqua 2005 (1-365)               
               years=[2009]; days={[1:365]};  %all there for Aqua 2005 (1-365)                 
               years=[2010]; days={[1:365]};  %all there for Aqua 2005 (1-365)                                
               years=[2011]; days={[1:365]};  %all there for Aqua 2005 (1-365)                                
               years=[2012]; days={[1:366]};  %all there for Aqua 2005 (1-365)                 
               years=[2012]; days={[366]};  %all there for Aqua 2005 (1-365)                                
               
            case 'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37'
                %whether to apply the water path
                % confidence screening in open_L2_MODIS_file_01 and whether to limit re16
                % and re37 to <20um as done for re21 when the water path confidence screening is
                % done
                i_apply_water_path_confidence = 1;
                i_limit_16_37_to_20_um = 1;  %N.B. will only want to apply this when doing water path screening
                
                filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-50; LAT_max=20; LON_min=-160; LON_max=-60; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        
                years=[2007];
%                years=[2010];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[1:365]};   
                
            case 'VOCALS_JointL2_box_only_confidence_screen_20um_limit_re16_37_consistent_pixels'
                % apply the water path
                % confidence screening in open_L2_MODIS_file_01 and limit re16
                % and re37 to <20um as done for re21 when the water path confidence screening is
                % done
                % Plus, discard all pixels where don't have a retrieval
                % (after screening) for ALL wavelengths
                i_apply_water_path_confidence = 1;
                i_limit_16_37_to_20_um = 1;  %N.B. will only want to apply this when doing water path screening
                i_make_wavelengths_consistent = 1;
                
                filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-50; LAT_max=20; LON_min=-160; LON_max=-60; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        
%                years=[2007];
%                years=[2010];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
%                days = {[1:365]};    
                
%               direcs={'aqua'};  years=[2009 2010 2010 2011 2011 2012]; days={[213:365],[1:121],[213:365],[1:121],[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                                                                             
               direcs={'aqua'};  years=[2005 2006 2007 2008 2009 2010 2011]; days={[1:365],[1:365],[1:365],[1:366],[1:365],[1:365],[1:365]}  %
                
                
            case 'Iceland_30W30E_40N90N'
                %For the PROCESS action
   
                filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=40; LAT_max=90; LON_min=-30; LON_max=30; %VOCALS limiting to the box
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                        


                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly  
                
%               years=[2008]; days={[1:351 355:361 363 365]}  %missing:- 352:354, 362, 364 and 366
                             
               years=[2014]; days={[245:245]};  %all there for Aqua 2005 (1-365)                                 
               years=[2014]; days={[242:253]};  %all there for Aqua 2005 (1-365)
                                
                
            case 'Arctic_summer_JointL2'
                time_restrict_case = 'choose here';
                override_mockL3_options2=1;
                LAT_min=40; LAT_max=90; LON_min=-180; LON_max=180; %Arctic summer 2007 files

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

                
                direcs={'terra','aqua'};
%                direcs={'aqua'};  
                direcs={'terra'};                  
                        
                years=[2006:2007];
                years=[2009];

                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[164:181]};
                
            case 'Southern_Ocean_all_lons_no_water_path_confidence' 
                % Process action
                %whether to apply the water path
                % confidence screening in open_L2_MODIS_file_01 and whether to limit re16
                % and re37 to <20um as done for re21 when the water path confidence screening is
                % done
                i_apply_water_path_confidence = 0;
                i_limit_16_37_to_20_um = 0;  %N.B. will only want to apply this when doing water path screening
                
                 filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-90; LAT_max=-20; LON_min=-180; LON_max=180; 
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

%                direcs={'terra','aqua'};
 
%                direcs={'terra'};                  
                        


                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[1:365]};   
                
               years=[2008]; days={[1:351 355:361 363 365]}  %missing:- 352:354, 362, 364 and 366
               years=[2008]; days={[324:324]}  %all there for Aqua 2008 (1-366)    
               
               years=[2006]; days={[1:365]}  %all there for Aqua 2006 (1-365)
               direcs={'aqua'};  years=[2007]; days={[1:365]}  %all there for Aqua 2007 (1-365)                    
               
               direcs={'aqua'};  years=[2006 2007]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)
%               direcs={'aqua'};  years=[2007 2008]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
%               direcs={'aqua'};  years=[2007 2008]; days={[339:365],[1:121]}  %must have crashed or something after doing day 338.
               direcs={'aqua'};  years=[2008 2009]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
               direcs={'aqua'};  years=[2009 2010]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2011 2012 2012]; days={[213:365],[1:121],[213:365]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2013 2013 2014 2014]; days={[1:121],[213:365],[1:121],[213:365]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2010 2011]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)    
               direcs={'aqua'};  years=[2005 2006]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)    
               
            case 'Southern_Ocean_all_lons_water_path_confidence_ON'
                % Process action
                %whether to apply the water path
                % confidence screening in open_L2_MODIS_file_01 and whether to limit re16
                % and re37 to <20um as done for re21 when the water path confidence screening is
                % done
                i_apply_water_path_confidence = 1;
                i_limit_16_37_to_20_um = 0;  %N.B. will only want to apply this when doing water path screening
                
                 filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1; 
                LAT_min=-90; LAT_max=-20; LON_min=-180; LON_max=180; 
                %to save processing time

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];

%                direcs={'terra','aqua'};
 
%                direcs={'terra'};                  
                        


                %days determines the size of the pre-allocated array (only when concatenating)
                %, so set it properly
                days = {[1:365]};   
                
               years=[2008]; days={[1:351 355:361 363 365]}  %missing:- 352:354, 362, 364 and 366
               years=[2008]; days={[324:324]}  %all there for Aqua 2008 (1-366)    
               
               years=[2006]; days={[1:365]}  %all there for Aqua 2006 (1-365)
               direcs={'aqua'};  years=[2007]; days={[1:365]}  %all there for Aqua 2007 (1-365)                    
               
               direcs={'aqua'};  years=[2006 2007]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)
%               direcs={'aqua'};  years=[2007 2008]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
%               direcs={'aqua'};  years=[2007 2008]; days={[339:365],[1:121]}  %must have crashed or something after doing day 338.
               direcs={'aqua'};  years=[2008 2009]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)               
               direcs={'aqua'};  years=[2009 2010]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2011 2012 2012]; days={[213:365],[1:121],[213:365]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2006 2013 2014 2014]; days={[213:365],[1:121],[213:365],[1:121],[213:365]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                              
               direcs={'aqua'};  years=[2010 2011]; days={[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)    
               direcs={'aqua'};  years=[2006 2007 2007 2008 2008 2009]; days={[213:365],[1:121],[213:365],[1:121],[213:365],[1:121]}  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                                             
               
            case 'global'
                filedir='/home/disk/eos5/d.grosvenor/joint_L2/';
                
                time_restrict_case = 'choose here';
                override_mockL3_options2=1;
                LAT_min=-90; LAT_max=90; LON_min=-180; LON_max=180; %Global                

                ioverride_multi_dir=1;
                processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];
                
                direcs={'terra','aqua'};
                direcs={'aqua'};  
%                direcs={'terra'};                  
                
%                years=[2006:2007]; days = {[306:365],[1:19 22:181]}; %
                %days 20-21 seem to be missing
%                years=[2007];
                
                
%                days = {[136:181]};
%                days = {[182:210]}; 
                
                direcs={'terra'};  years=[2007]; days = {[182:210]}; %
                 direcs={'aqua'};  years=[2006]; days = {[182:210]}; %               
%                direcs={'aqua'};  years=[2007]; days = {[211:243]}; %
%                direcs={'terra'};  years=[2006 2007]; days = {[333:365],[1:19 22:332]}; %data processed so far
                
%                direcs={'terra'};  years=[2006 2007]; days = {[333:365],[1:19 22:332]}; %data processed so far                
%                direcs={'aqua'};  years=[2007 2008]; days = {[1:365],[1:316]}; %started 5th March 2013
%                direcs={'terra'};  years=[2007]; days = {[1:19 22:365]}; %       
%                direcs={'terra'};  years=[2007]; days = {[333:365]}; %  
%                direcs={'terra'};  years=[2007]; days = {[1:1]}; % 
%                direcs={'aqua'};  years=[2007]; days = {[354]}; % 20th Dec                
%                direcs={'aqua'};  years=[2007]; days = {[171]}; % 20th June                                
%                direcs={'aqua'}; years=[2007]; days = {[131:210 244:365]};

               direcs={'aqua'};  years=[2006 2007 2007 2008 2008 2009]; days={[213:365],[1:121],[213:365],[1:121],[213:365],[1:121]};  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                                             
               direcs={'aqua'};  years=[2009 2010 2010 2011 2011 2012]; days={[213:365],[1:121],[213:365],[1:121],[213:365],[1:121]};  %processing just when have SZA<65 data (01Aug06 - 30Apr07)                                                             
               direcs={'aqua'};  years=[2005 2005 2006]; days={[21:121],[213:365],[1:121]};
               direcs={'aqua'};  years=[2008 2008]; days={[153:212],[122:152]};                 
               direcs={'aqua'};  years=[2006]; days={[122:212]};         
               direcs={'aqua'};  years=[2005 2009 2010 2012]; days={[163:212],[122:212],[122:212],[122:366]};               
               
               direcs={'aqua'};  years=[2003 2004 2015]; days={[1:365],[1:366],[1:365]};                                
               direcs={'aqua'};  years=[2015 2011]; days={[20:365],[122:212]};                                               
%               direcs={'aqua'};  years=[2013 2014]; days={[18:365],[1:365]};                                
%               direcs={'aqua'};  years=[2003 2004]; days={[10:365],[1:366]};

            direcs={'aqua'};  years=[2009 2013]; days={[122:212],[18:365]}; %Challenger             
            direcs={'aqua'};  years=[2010 2014]; days={[137:212],[1:365]}; %Olympus
            direcs={'aqua'};  years=[2010 2011]; days={[172:212],[155:212]}; %Challenger (to help out Olympus)           
            direcs={'aqua'};  years=[2011 2015]; days={[130:212],[20:365]}; %Olympus 
            direcs={'aqua'};  years=[2012]; days={[133:366]}; %Challenger
            %direcs={'aqua'};  years=[2004]; days={[1:366]}; %Challenger
            %direcs={'aqua'};  years=[2003]; days={[10:365]}; %        
            %direcs={'aqua'};  years=[2013]; days={[18:365]}; %Challenger batch 
            %direcs={'aqua'};  years=[2014]; days={[1:365]}; %Olympus batch
            %direcs={'aqua'};  years=[2015]; days={[20:365]}; %Olympus batch  - stopped 
            direcs={'aqua'};  years=[2011]; days={[170:212]}; %Challenger - DONE
direcs={'aqua'};  years=[2010]; days={[168:212]}; %Challenger - was an issue with 167 (out of range subscript)?
%direcs={'aqua'};  years=[2010]; days={[167]}; %Challenger - was an issue with 167 (out of range subscript)?

direcs={'aqua'};  years=[2010]; days={[1:365]}; %
direcs={'aqua'};  years=[2009]; days={[1:365]}; %

%direcs={'aqua'};  years=[2012]; days={[326:366]}; %olympus direct
%direcs={'aqua'};  years=[2013]; days={[280:365]}; %olympus direct
%direcs={'aqua'};  years=[2014]; days={[162:365]}; %Challenger direct
%direcs={'aqua'};  years=[2014]; days={[15:159]}; %Challenger direct
%direcs={'aqua'};  years=[2015]; days={[15:159]}; %Challenger direct
%direcs={'aqua'};  years=[2015]; days={[100:192]}; %Challenger direct
%direcs={'aqua'};  years=[2015]; days={[198:365]}; %Challenger direct
%direcs={'aqua'};  years=[2003]; days={[42:100]}; %Challenger direct             
%direcs={'aqua'};  years=[2003]; days={[311:365]}; %Challenger direct     
%direcs={'aqua'};  years=[2008]; days={[366]}; %Challenger direct   


%2003 - 11:100; 101:188; 189:277; 278:365

               %McCoy wants 2003-2015
               % After the ones above will need full years for 2003-2004 &
               % 2013-2015
               % Processing 2003-2004; 2015; 2013-2014; 2005-2012
               
               %2005 done. as is 2006, 2007, 2008. Need 2009-2013
               
            case {'Puijo_Sami','Puijo_Sami_DMPS_L2'}
                
               
                                 
                MODIS_filetype_choose = 'L2'
%                MODIS_filetype_choose = 'L2 Joint'     


                disp('Make sure that MODIS_filetype_choose is correct !! ');
                
                switch MODIS_filetype_choose
                    case 'L2'
                        filedir='/home/disk/eos8/d.grosvenor/MOD_L2/Finland_Sami/'; %for L2 CDP
                        filedir='/home/disk/eos8/d.grosvenor/MOD_L2/Finland_Sami/DMPS_matches/'; %for L2 DMPS Nacc                        
                        processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '/'];  
                    case 'L2 Joint'
                        filedir='/home/disk/eos5/d.grosvenor/joint_L2/Puijo/'; %for jointL2
                        processed_daily_swath_loc=['/home/disk/eos8/d.grosvenor/saved_data_L2/' multiL2L3_project '_jointL2/'];
                end
                


                Puijo_box_size = '0.5 x 1 deg (lat x lon)';
%                Puijo_box_size = 'NxNkm';
                
                 Npix_Puijo = 5; %NxN km
                 Npix_Puijo = 11;

                
                %whether to only open files that match certain times
                time_restrict=1;  dtime_restrict=0.5; %in hours: can be +/- this many hours
%                time_restrict=1;  dtime_restrict=2/60; %for matching the CDP 5mins averages    
%                time_restrict=1;  dtime_restrict=15/60; %for matching the CDP 5mins averages                    
                use_wind_Irshad = 0; %flag to say whether to use the wind speed and direction or not
                time_restrict_case = 'Puijo station';
                
                override_mockL3_options2=1;
                LAT_min = 62; LAT_max = 64; 
                LON_min = 26; LON_max = 29; %Puijo_Sami : 62.91N, 27.66 E

                ialternate_file_structure=1;                
                ioverride_multi_dir=1;
                               
                
%The Puijo data is not in separate dirs for each day. So, set days to =1
%and load all together
                direcs={'Puijo/terra','Puijo/aqua'};  years=[2006:2012]; days = {1}; %
                
                
                direcs={'terra','aqua'};  years=[2008:2012]; days = {1}; %                
                %                direcs={'terra'};  years=[2006]; days = {1}; %
%                direcs={'aqua'};  years=[2009:2012]; days = {1}; %
                 

                
        end

end

%% End of option setting


switch multiL2L3_case
    case 'save concatenated mockL3 data'
        %don't do anything except save at the end
    
    otherwise

if length(days)<length(years)
    if length(days)>1
        disp('days and years mismatch!! Set days to correspond to each year');
        return
    elseif length(days)==1
        ain='n';
        while strcmp(ain,'y')~=1
            ain = input('\n**** days{} is only one set, expanding to all years - enter ''y'' to continue : ','s');
        end
        for idays_expand=2:length(years)
            days{idays_expand}=days{1};
        end
    end
end

fprintf(1,['\n---' multiL2L3_case ' action selected! ---']);
if exist('i_override_MODIS_multi')~=1
    ain='n';
%    while strcmp(ain,'y')~=1
%        ain = input('\n**** Don''t U forget to list the variables that you want to save \nin make_mockL3_variables.m - enter ''y'' to continue : ','s');
%    end
end


clear daynum_timeseries3 aqua_terra_timeseries3 modisyear_timeseries3


switch multiL2L3_case
    case {'load, average and save','load processed L2 and concatenate'}
%load/use SSTs for everything except 'load L3 and concatenate'
%taken it off for the 'process' action too


if iload_amsre_mat==1  %if we have saved the data already
    %otherwise need to run multi_read_amsre_daily to load the amsre
    %data - loading avoids the lengthy read in of the AMSRE data from the
    %individiual AMSRE files
    load(amsre_matfile); %loads in day_amsre, gcm_Plat2D_AMSRE (etc), year_amsre, 
    %month_amsre, sst_amsre_smooth
  
end

%pick out the correct days for those that we want to process from the sst
%data (uses sst_amsre_smooth, day_amsre, year_amsre etc).
    iday_all_L2L3=0;
    clear modisyear_timeseries3_block daynum_timeseries3_block
    %loop through all years and days to create daynum_timeseries3 and
    %modisyear_timeseries3 for amsre sst processing file
    for idirec=1:length(direcs)

        for iyear_L2L3=1:length(years)
            days_L2L3 = days{iyear_L2L3};
            for iday_L2L3=1:length(days_L2L3)
                iday_all_L2L3=iday_all_L2L3+1;
                modisyear_timeseries3_block(iday_all_L2L3) = years(iyear_L2L3);
                daynum_timeseries3_block(iday_all_L2L3) = days_L2L3(iday_L2L3);
            end
        end

    end
    %N.B. - these daynum_timeseries and modisyear_timeseries values may
    %change below (e.g. may be replicated for all orbits to make bigger
    %arrays

    % --- sort out the SST from AMSRE for this one day -------
    ioverride_block_amsre=1;
    amsre_block_average_time %this uses modisyear_timeseries3_block for the size
      %N.B. - this is also run later (after loading) for the 'load processed L2 and concatenate','load L3 and concatenate'
      %processes in order to create an SST field the same size as the data loaded in
      %(Ndays*20).
    % --------------------------------------------------------

sst_amsre_time3_MODIS_multi_day = sst_amsre_time3;
             
end             
             

day_range_catL2=0;  %this is the index count for all data (inc. both terra and aqua directories)
iday_all_L2L3=0;
%% directory loop (e.g. aqua and terra)
for idirec=1:length(direcs)
    %where to load from
     savedir_var = [processed_daily_swath_loc direcs{idirec} '/'];
     

     switch multiL2L3_case
         case 'process'
             %need to make one directory before it will make the next
%              if exist([processed_daily_swath_loc])~=7 %if not a directory
%                  eval(['!mkdir ' processed_daily_swath_loc]); %make it
%              end
%              if exist([savedir_var])~=7 %if not a directory
%                  eval(['!mkdir ' savedir_var]); %make it
%              end
             
             make_dir(savedir_var); %funciton that makes the directory (if it doesn't exist)

         case {'load, average and save','load L3 and concatenate','load processed L2 and concatenate'}

             daily_averaged_files_loc = [processed_daily_swath_loc direcs{idirec} '/'];

             switch multiL2L3_case
                 case 'load, average and save'
                     if exist([daily_averaged_files_loc2])~=7 %if not a directory
                         eval(['!mkdir ' daily_averaged_files_loc2]); %make it
                     end
             end
          
     end   %switch multiL2L3_case

     Ndays_L2L3 = length([days{:}]); %this gives the total number of days in the days array


%% year loop     
for iyear_L2L3=1:length(years)
    days_L2L3 = days{iyear_L2L3};
      
%% day loop
    for iday_L2L3=1:length(days_L2L3)
        iday_all_L2L3=iday_all_L2L3+1;
        day_str = num2str(days_L2L3(iday_L2L3));
     
        %add zeros to the front until we have 3 digits
        Lstr=length(day_str);
        for istr=Lstr:2
            day_str=['0' day_str];
        end
        
        year_str=num2str(years(iyear_L2L3));
        
        if ialternate_file_structure==1
            filedir2 = [direcs{idirec} '/' year_str '/']            
        else
            filedir2 = [direcs{idirec} '/' year_str '/' day_str '/']
        end
        
        switch multiL2L3_case
            case {'load, average and save','load L3 and concatenate','load processed L2 and concatenate'}
                switch multiL2L3_case
                    case 'load L3 and concatenate'
                        savedir_var = [daily_averaged_files_loc2 direcs{idirec} '/'];
                end
                
                %get the saved variables that fit the required day
                switch multiL2L3_project
                    case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2'}
                        saved_files = dir([savedir_var '*' year_str '*' Puijo_keyword '*.mat']);
                    otherwise
                        saved_files = dir([savedir_var '*' year_str '_' day_str '*make mock L3 data daily*.mat']); %think this get mangled?
                        saved_files = dir([savedir_var '*' year_str '_' day_str '*.mat']);                        
                end
                
                
                %find the most recent version
                %list of the modification dates in Matlab datenum format
                % -- had to switch to a different method here due to some
                % symbolic links being used - maybe something that Marc did
                % - that led to some older files being given e.g. 2038
                % timestamps...!
                   %file_times = [saved_files.datenum];
                % -- will use the timestamp in the filename instead
                clear file_times
                for ifile_time=1:length(saved_files)
                    time_str = saved_files(ifile_time).name(end-18:end-4);  
                    file_times(ifile_time) = datenum(time_str,'yyyymmddTHHMMSS');                    
                end
                
                [times_sorted,isort] = sort(file_times);
                %pick the last one (latest file)
                load_file = saved_files(isort(end)).name;

                L2L3_load_file = [savedir_var load_file];     
                
%%% --------------------------------------------------------------------------                                
                MODIS_process_global_L2_daily_with_screening
%%% --------------------------------------------------------------------------                
                
            case 'process'
                ioverride_multi_dir=1;

                % action='store data etc';
                % action='draw and save plots';
                % action='make mock L3 data';
                action='make mock L3 data daily';                                             

                idir=1;
                modis_dir_multi{idir}=filedir2; MODIS_filetype = MODIS_filetype_choose; idir=idir+1;

%----------------------------------------------------------------------
                MODIS_multi_dir_process %for the 'process' action
%----------------------------------------------------------------------                
        end
        
        switch multiL2L3_case
            case 'load processed L2 and concatenate'
                modisyear_timeseries3(day_range_catL2)=years(iyear_L2L3);
            otherwise
                %moved this to inside MODIS_process_global_L2_daily_with_screening
                %modisyear_timeseries3(iday_all_L2L3)=years(iyear_L2L3);
        end

    end  %for iday_L2L3=1:length(days_L2L3)

end %for iyear_L2L3=1:length(years)

end  %for idirec=1:length(direcs)


end  %switch multiL2L3_case



switch multiL2L3_case
   
        
    case {'load processed L2 and concatenate','load L3 and concatenate'}
 %       daynum_timeseries3=[1:length(Cloud_Fraction_Liquid.timeseries3)];
        aqua_terra_timeseries3=upper(direcs);
        
        LAT=MLAT(ilat);
        LON=MLON(ilon);
        
        MLAT=LAT;
        MLON=LON;
        %LAT_min=40; LAT_max=90; LON_min=-180; LON_max=180;
        mod_data_type='timeseries3';
        
        
        if post_process_L2_and_L3==1
            switch multiL2L3_case

               %load ssts for L3 concatenate as didn't do earlier
                case {'load L3 and concatenate'}
                    if iload_amsre_mat==1  %if we have saved the data already
                        %otherwise need to run multi_read_amsre_daily to load the amsre
                        %data - loading avoids the lengthy read in of the AMSRE data from the
                        %individiual AMSRE files
                        load(amsre_matfile); %loads in day_amsre, gcm_Plat2D_AMSRE (etc), year_amsre,
                        %month_amsre, sst_amsre_smooth
                    end
            end
            
            switch multiL2L3_project
                case {'Puijo_Sami','Puijo_Sami_jointL2','Puijo_Sami_DMPS_L2'}
                    %don't do this for Puijo (is over land and only have one location)
                    sst_amsre_time3 = NaN*ones(size(Cloud_Fraction_Liquid.timeseries3));
                otherwise  %to switch off sst matching set ino_sst_matching=1 in the command window before running
                    %N.B. will need to do this every time as it is reset
                    %upon running
                    if ~exist('ino_sst_matching') | ino_sst_matching==0
                        amsre_block_average_time %this uses modisyear_timeseries3 for the size
                    end
            end

            filtering_data_get

            %        mod_data_type_bk = mod_data_type;
            %        mod_data_type = 'timeseries3';

            [N_time3_16]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_16_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');
            [N_time3_37]=MODIS_justN_func(Cloud_Optical_Thickness_Liquid_Mean.timeseries3,Cloud_Effective_Radius_37_Liquid_Mean.timeseries3*1e-6,Wflag,0,Cloud_Top_Temperature_Day_Mean.timeseries3,'N');

            %        mod_data_type = mod_data_type_bk;

        end
        
    case 'save concatenated mockL3 data' %just save mock L3 data in the directory of the screened data
        savedir_var= daily_averaged_files_loc2;
        %    savevarname = [savedir_var 'timeseries3_TProfiles_CTT_CTP_data_' dataset_modis '_' datestr(now,30)];
        savevarname = [savedir_var 'mockL3_saved_data_' datestr(now,30)];
        
        %Add N_time3, W_time3 and H_time3 to the list to save
        imod_var=length(modis_var)+1;
        modis_var{imod_var} = 'N_time3_16'; imod_var=imod_var+1;
        modis_var{imod_var} = 'N_time3_37'; imod_var=imod_var+1;        
        modis_var{imod_var} = 'N_time3'; imod_var=imod_var+1;
        modis_var{imod_var} = 'W_time3'; imod_var=imod_var+1;
        modis_var{imod_var} = 'H_time3'; imod_var=imod_var+1;

        MLAT = MLAT_full(ilat_loaded);
        MLON = MLON_full(ilon_loaded);        
        
       % --------------------------------------
              save_MODIS_L3_data
       % --------------------------------------

end

switch multiL2L3_case
    case 'process'
    otherwise
        if exist('daynum_timeseries3')            
            daynum_timeseries3_MODIS = daynum_timeseries3;
        end
        gcm_str_select = 'MODIS';
        gcm_str = 'MODIS';

        LAT_val_loaded = LAT_val;
        LON_val_loaded = LON_val;

        %Here we find the incides of the locations laoded in for the full lat lon
        %arrays - this routine can load in only select areas, so the ilat and ilon
        %are indices for those rather than for the full global lat lon
        MLAT_full=[-89.5:89.5];
        MLON_full=[-179.5:179.5];
        
        if exist('LAT')

            clear ilat_loaded ilon_loaded
            for i=1:length(LAT)                
                 ifind = find(MLAT_full == LAT(i));
                 if length(ifind)>0
                     ilat_loaded(i) = ifind;
                 else
                     ilat_loaded(i) = 0;
                 end
            end
            for i=1:length(LON)
                ifind = find(MLON_full == LON(i));
                if length(ifind)>0
                    ilon_loaded(i) = ifind;
                else
                    ilon_loaded(i) = 0;
                end
            end

        end

end

clear irun_from_modis_multi_DAY_processL3L2 ino_sst_matching i_override_MODIS_multi
catch multi_DAY_error
    clear irun_from_modis_multi_DAY_processL3L2 ino_sst_matching i_override_MODIS_multi
    rethrow(multi_DAY_error)
end




