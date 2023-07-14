%read and average multiple modis files from multiple directories
%Steps for running / things to make sure of:-
% 1) *** Don't forget to list the variables that you want to save in
% make_mockL3_variables.m
% 2) Set the required lat/lon boundaries in mockL3_from_L2
% 3) Set time_restrict below - =0 if don't want to choose files at certain times/locations  
% 4) Set the directory of files required below - uncomment all not required
% as otherwise it processes multiples ones

%runs this  - MODIS_process_multiple_L2_files

try

if ~exist('ioverride_multi_dir') | ioverride_multi_dir==0
    
    savedir_var=['/home/disk/eos1/d.grosvenor/modis_work/saved_data_L2/'];
    savedir_var=['/home/disk/eos8/d.grosvenor/saved_data_L2/Arctic_summer_JointL2/'];
    savedir_var=['/home/disk/eos8/d.grosvenor/AP_Feb_2010_daytime/06Feb2010_flight99/aqua/images/'];
    
    action='store data etc';
    action='draw and save plots';
    %action='make mock L3 data';
    %action='make mock L3 data daily';

    
    ain='n';
    while strcmp(ain,'y')~=1
        ain = input('\n**** Don''t U forget to list the variables that you want to save \nin make_mockL3_variables.m - enter ''y'' to continue : ','s');
    end

    filedir='/home/disk/eos8/d.grosvenor/';
    filedir='';    

    %modis_dir='MPACE_L2_17days/aqua/';
    %modis_dir='MPACE_L2_17days/';

    %make sure to comment out directories that we don't want to process
    % *** THIS WILL PROCESS ALL THE SELECTED DIRECTORIES!! ***
    idir=1;
    %modis_dir_multi{idir}='MPACE_JointL2_17days_all_lats/terra/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='MPACE_JointL2_17days_all_lats/aqua/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='MPACE_L2_17days_all_lats/terra/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='MPACE_L2_17days_all_lats/aqua/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='summer_2004_17days_all_lats/aqua/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='summer_2004_17days_all_lats/terra/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/TERRA/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='VOCALS/10-40S_100-60W_12Oct-17Nov_2008/AQUA/'; MODIS_filetype = 'L2'; idir=idir+1;
    %modis_dir_multi{idir}='AP_Dec_2009_daytime/AQUA/'; MODIS_filetype = 'L2'; idir=idir+1; %Antarctic summer swaths
    %modis_dir_multi{idir}='summer_2004_17days_all_lats/aqua/'; MODIS_filetype = 'L2'; idir=idir+1; %Arctic summer (2004) swaths
    %days 164-181
    %modis_dir_multi{idir}='Arctic_summerL2_2007_20W-60E_70-80N/terra/'; MODIS_filetype = 'L2'; idir=idir+1; %New Arctic summer (2007) swaths
    %modis_dir_multi{idir}='Arctic_summerL2_2007_20W-60E_70-80N/aqua/'; MODIS_filetype = 'L2'; idir=idir+1; %New Arctic summer (2007) swaths
    %modis_dir_multi{idir}='Arctic_summerL2_2007_20W-60E_70-80N/Joint_5km_files/terra/'; MODIS_filetype = 'L2 Joint'; idir=idir+1; %New Arctic summer (2007) swaths - Joint 5 km data (subsampled)
    %modis_dir_multi{idir}='MOD_L2/2012_06_05_hole_Australia/'; MODIS_filetype = 'L2'; idir=idir+1; %"Hole near Tasmania"
    %modis_dir_multi{idir}='joint_L2/terra/2007/001/'; MODIS_filetype = 'L2 Joint'; idir=idir+1; %day 001 of 2007 L2 joint files

    %modis_dir_multi{idir}='AP_Feb_2010_daytime/AQUA/'; MODIS_filetype = 'L2'; idir=idir+1; %Antarctic Flight 104 12th Feb 2010 swaths
    %modis_dir_multi{idir}='AP_Feb_2010_daytime/TERRA/'; MODIS_filetype = 'L2'; idir=idir+1; %have all TERRA, flight 104 with AQUA

%    modis_dir_multi{idir}='AP_Feb_2010_daytime/06Feb2010_flight99/aqua/'; MODIS_filetype = 'L2'; idir=idir+1; %Flight 99 aqua swaths
    
%    modis_dir_multi{idir}='MOD_L2/Iceland_2014/terra/249/';  MODIS_filetype = 'L2'; idir=idir+1; 6th Sep
    %modis_dir_multi{idir}='MOD_L2/Iceland_2014/aqua/';  MODIS_filetype = 'L2'; idir=idir+1;  
%    modis_dir_multi{idir}='MOD_L2/Iceland_2014/terra/Sep2014/';  MODIS_filetype = 'L2'; idir=idir+1; 
%    modis_dir_multi{idir}='MOD_L2/Iceland_2014/terra/01-07Sep2014/';  MODIS_filetype = 'L2'; idir=idir+1;     
    %modis_dir_multi{idir}='MOD_L2/Iceland_2014/terra/01-07Sep2014/stripped_down_for_nicer_plot/';  MODIS_filetype = 'L2'; idir=idir+1;      
%modis_dir_multi{idir}='/home/disk/eos15/d.grosvenor/eos8/Koike_Japan/';  MODIS_filetype = 'L2'; idir=idir+1;        
modis_dir_multi{idir}='/home/disk/eos15/d.grosvenor/eos8/MOD_L2/ship_tracks_UK_Spain_NdReview/';  MODIS_filetype = 'L2_C6'; idir=idir+1;     



%    modis_dir_multi{idir}='MOD_L2/Azores_22Nov_2009/';  MODIS_filetype = 'L2'; idir=idir+1;       
    % *** make sure to set the lats and lons required for the locations to be stored
    %     in MODIS_average...
    % ***************************************************************************

    
    %whether to only open files that match certain times
    time_restrict=0;  dtime_restrict=1; %in hours: can be +/- this many hours
    
    time_restrict_case = 'choose here';
    %time_restrict_case = 'VOCALS cloud segments (Chris T)';

    if time_restrict==0
        time_restrict_case = 'choose here';
    end


end

isave_MODIS_data=1;  %flag to say whether to save all of the data gathered or not
climits=[0 500];



if time_restrict==1
%     ain='n';
%     while strcmp(ain,'y')~=1
%         ain = input('\n**** time_restrict==1 - is this correct? \nEnter ''y'' to continue : ','s');
%     end


disp('*** WARNING time_restrict==1 !! ***')

end




day_only=0; %flag to stop the read in if the time is nighttime (estimated on time - for VOCALS only)
composite_single_days=1; %flag to plot all the images for each day on one plot instead of all of the 
%swaths individually.
iclose_figures=0; %flag to close the figures once saved

% ---------------------------------------------

for idir=1:length(modis_dir_multi)
    modis_dir=modis_dir_multi{idir};
    override_flags_openL2 = 1; %flag to say we want to override the flags in open_L2_MODIS_file_01
    MODIS_process_multiple_L2_files
end
% ---------------------------------------------

clear ioverride_multi_dir
catch multi_error
    clear ioverride_multi_dir
    rethrow(multi_error)
end