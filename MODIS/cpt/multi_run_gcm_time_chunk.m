try
%make sure that the COSP flag is set below   
%
%runs:-
%   am3_multi_load_time_chunks - runs it multiple times if required
%     set the model run name below
%   am3_choose_load_file - set the files/dir to proces here

%NOTE - after running this copy the filename for the .mat file that is
%printed at the end in load_catfiles_timechunks

% gcm_list_of_saved_variables   - variables to save are set in this

% Also runs 
%     am3_choose_load_file  -- this contains the directories for the
%         raw GCM data. Select what years and files to process here.

% scripts for reading the data for the different models are repeated below
% from am3_multi_load_time_chunks.m   :-

%  AM3 COSP
%             iread_cosp=1;
%             read_AM3T_cosp_savemem
%  AM3 no COSP
%             iread_cosp=0;
%             read_AM3T_cosp_savemem       
%  CAM no COSP or COSP
%             read_cosp_cam
%             read_cam
%
%  This script processes each file into several .mat files (split by time chunks) and then the filenames
%  of these are stored in am3_filepaths_chunk. This is a structure with
%  each entry the filename of the "chunk"
%  am3_filepaths_chunk.y2010_CAM5_COSP_10_chunk_5
%  ans =
%  saved_gcm_processed_data_CAM5_COSP_CAM5_new_1deg_COSP_cam5_1_26_AMIP_1deg.cam.h1.2010-10-11-64800.nc_20130725T065734.mat
%  am3_filepaths_chunk, along with variables containing the years and
%  months of the chunks are then saved in another .mat file -
%  gcm_saved_filename_file.
%  One of these is created for each year of data.

temp_gz_file = '/home/disk/eos5/d.grosvenor/temp_gz';
temp_gz_inst = [temp_gz_file '.nc'];
temp_gz_grid = [temp_gz_file '_grid.nc'];
temp_gz_prepostLWP = [temp_gz_file '_prepostLWP.nc'];

tempval=input('Make sure that the COSP flag is set correctly!! Press enter to continue...');

ioverride_multi_chunk=1;
% am3_choose_load_file is where to choose the required file
cosp_flag=1; %flag to save the COSP variables.
ino_dbz_cfads=1; %if set to one then DOESN'T read in the radar CFADS (cloudsat)
%  - not available for the post diurnal CAMCLUBB runs (inc. prepost LWP)
ilwcAPBP=1; %set to one if we have the pre and post mphys LWP (zero otherwise)
  %Not sure if this does anything at the mo

nT = 10; %choose the number of time indices to load and  process at one time
%can now choose nT=1 if required (didn't work previously)- at least for
%CAM5
%nT = 20;

%% Select the lat lon range to extract - N.B. - only applies to CAM models at
%present. For AM3 just one of the tiles is loaded - see below
       %VOCALS region:-
    lat_range = [-40 10];
    lon_range = [-140 -50]+360;

%For AM3 select one tile to process    
%all have same lat range +44.2 to -44.2
%except two polar sterographic plots
%1= 305.6-34.2 (centred on lon=0)
%2= 35.8 - 124.2
%3 = north pole spherical 90 to 36 in lat approx 
%4 =  125.8-214.2 in lon
%5 is the best for the VOCALS-CAPT region, 215.8-304.2 (or -144.2 to -55.8)
%6 = South pole -36 to -90 in lat
    tile_no = '5';  %N.B. may not have data for the other tiles for a lot of runs
    
    
filedir_gcm_load = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/'; %Directory to save in
%- need another dir in here for master .mat files

if cosp_flag==0
    cosp_flag4D=0;
    iread_cosp=0;
else
    cosp_flag4D=1;
    iread_cosp=1;
end

%chunk_model is used to determine the style of filename for the output -
%the year is obtained from the filename later in this script, depending on its location in
%the filename string. Will need to set a new type if the filename
%convention is changed.
%  Is also saved in gcm_month_cases{i} variable, which is one of the variables saved in the
%  final "master" file for each year (gcm_saved_filename_file), along with
%  am3_filepaths_chunk, etc.

chunk_model='AM3'; %choose am3_dataset below
%chunk_model='AM3_CLUBB';
%chunk_model='CAM5_COSP';
%chunk_model='CAM5';
%chunk_model='CAM5_CLUBB';
%chunk_model='CAM5_CLUBB_COSP';  %CLUBBv1
%chunk_model='CAM5_CLUBBv2_COSP';
%chunk_model='AM3_CLUBBv2'; %CLUBBv2, 50km
%chunk_model='AM3_CLUBBv2_COSP_100km'; %CLUBBv2, 100km
%chunk_model='AM3_CLUBBv2_COSP_200km'; %CLUBBv2, 200km (/home/disk/eos8/d.grosvenor/CPT/AM3_COSPCLUBB/AM3_CLUBB_200km_COSP_20130715)
%chunk_model='AM3_CLUBBv2_COSP_50km'; %CLUBBv2, 200km (/home/disk/eos8/d.grosvenor/CPT/AM3_COSPCLUBB/AM3_CLUBB_50km_COSP_20130715)
%chunk_model = 'AM3_CLUBBv1_1deg';
%chunk_model = 'AM3_CLUBBv1_2deg'; %no COSP
%chunk_model = 'AM3_100km_20130110';
chunk_model = 'CAM5_prepostLWP'; %
%chunk_model = 'CAMCLUBBv2_prepostLWP';

switch chunk_model
    case 'AM3'
        am3_dataset = 'AM3new_2deg';
        am3_dataset = 'AM3new_0.5deg';
        suffix = '';
    case 'AM3_CLUBB'
        am3_dataset = 'AM3_CLUBB';  
        suffix = '';
    case 'CAM5_CLUBB'
        am3_dataset = 'CAM5_CLUBB_new_1deg';
         suffix = '';
    case 'CAM5_CLUBB_COSP'
        am3_dataset = 'CAM5_CLUBB_new_1deg';  %original CLUBB output on which I revealed the diurnal issue
        suffix = '';
    case 'CAM5_CLUBBv2_COSP';
        am3_dataset = 'CAMCLUBB_COSP_postDiurnal'; %revised runs with attempt to fix diurnal issue (Jan 2013)
        %Changed the auto/accretion balance
        suffix = '_LON_150w_to_55w_LAT_45s_to_45n'; %this is needed because the variable names in the netCDF changed
    case 'AM3_CLUBBv2';
        am3_dataset = 'AM3_CLUBB_postDiurnal'; %revised runs with attempt to fix diurnal issue (Jan 2013)
        %Changed the auto/accretion balance        
    case 'AM3_CLUBBv2_COSP_100km'
        am3_dataset = 'AM3_CLUBB_COSP_prepostLWP'; %same as postdiurnal except now with COSP
    case 'AM3_CLUBBv1_1deg'
        am3_dataset = 'AM3_CLUBBv1_1deg';
    case 'AM3_CLUBBv1_2deg'
        am3_dataset = 'AM3_CLUBBv1_2deg';        
    case 'AM3_100km_20130110'; %AM3 base 100km res, w/COSP 2006-2010
        am3_dataset = 'AM3_100km_20130110';
    case 'AM3_CLUBBv2_COSP_200km'
        am3_dataset = 'AM3_CLUBBv2_2deg';
    case 'AM3_CLUBBv2_COSP_50km'
        am3_dataset = 'AM3_CLUBBv2_50km';        
    case 'CAM5'
        am3_dataset = 'CAM5_new_1deg';
        %am3_dataset = 'CAM5_old_2deg';
         suffix = '';
    case 'CAM5_COSP'
        am3_dataset = 'CAM5_new_1deg_COSP';  
        suffix = '';
    case 'CAM5_prepostLWP'
        am3_dataset = chunk_model;  %could just make this default so don't have set every time
        suffix = '';
    case 'CAMCLUBBv2_prepostLWP'
        am3_dataset = chunk_model;  %could just make this default so don't have set every time
        suffix = '';    
        
end

%% New 'master' list of gcm_saved_filename_file strings - just have to
%% copy the name of this to the loading routine now
datestrnow_file_list = datestr(now,30); %datestamp to create a unique file
save_file_list_filename = [filedir_gcm_load 'save_file_list_timechunks_' am3_dataset '_' datestrnow_file_list '.mat'];
clear load_file load_file_year


iread_normal_vars=1; 

% ------  file is chosen here --------- 
    clear files_multi
    ioverride_am3choose=1;
    am3_choose_load_file
% -------------------------------------

if ~exist('files_multi')
    nfiles_multi = 1;
else
    nfiles_multi = length(files_multi);
end

clear am3_filepaths_chunk
itam32=1; %counter for gcm_year_cases etc
clear gcm_year_cases gcm_model_cases gcm_month_cases

isave_file = 1;
file_year_old='';

if nfiles_multi==0
    error(' *** No files - nfiles_multi=0 ***');
end
for imulti_run = 1:nfiles_multi
    
    if exist('files_multi')
        nc_inst_file = files_multi(imulti_run).name;
    end
    
    if exist('files_multi_grid')
        nc_grid_file = files_multi_grid(imulti_run).name;
    end
    
    if exist('files_multi_prepostLWP')
        nc_prepostLWP_file = files_multi_prepostLWP(imulti_run).name;
    end
    
    
    
    % -------  Uzip (if  required) and open the NetCDF files --------------------
       
    nc_filepath =  [nc_dir nc_inst_file];  
    unzipped=0;
    if strcmp(nc_inst_file(end-2:end),'.gz')==1
        if exist(temp_gz_inst)==2
            eval(['!rm -f ' temp_gz_inst]);
        end
        eval(['!gunzip -c ' nc_dir nc_inst_file ' > ' temp_gz_inst]);
%        nc_inst_file = nc_inst_file(1:end-3); %remove the .gz
        nc_filepath = temp_gz_inst ;
        unzipped=1;
    end
    nc_inst = netcdf([nc_filepath],'nowrite');

       
    nc_filepath_grid = [nc_dir_grid nc_grid_file];
    unzipped_grid=0;
    if exist('nc_grid_file')
        if strcmp(nc_grid_file(end-2:end),'.gz')==1
            if exist([temp_gz_grid])==2
                eval(['!rm -f ' temp_gz_grid]);
            end
            eval(['!gunzip -c ' nc_dir_grid nc_grid_file ' > ' temp_gz_grid]);
%            nc_grid_file = nc_grid_file(1:end-3); %remove the .gz
            nc_filepath_grid = [temp_gz_grid];
            unzipped_grid=1;
        end
        nc_grid = netcdf([nc_filepath_grid],'nowrite');
    end
    
    nc_filepath_prepostLWP = [nc_dir_prepostLWP nc_prepostLWP_file];
    unzipped_prepostLWP=0;
    if exist('nc_prepostLWP_file')
        if strcmp(nc_prepostLWP_file(end-2:end),'.gz')==1
            if exist([temp_gz_prepostLWP])==2
                eval(['!rm -f ' temp_gz_prepostLWP]);
            end
            eval(['!gunzip -c ' nc_dir_prepostLWP nc_prepostLWP_file ' > ' temp_gz_prepostLWP]);
%            nc_grid_file = nc_grid_file(1:end-3); %remove the .gz
            nc_filepath_prepostLWP = [temp_gz_prepostLWP];
            unzipped_prepostLWP=1;
        end
        nc_prepostLWP = netcdf([nc_filepath_prepostLWP],'nowrite');
    end
    
    
    
    % --------------------------------------------------
    
    
%---------  Get the year and month from the filename ---------------    
    if length(strfind(chunk_model,'AM3'))>0
        file_year = nc_inst_file(1:4);
        file_month = nc_inst_file(5:8);
    elseif length(strfind(chunk_model,'CAM'))>0
        istr=strfind(nc_inst_file,'.cam.h')
        file_year = nc_inst_file(istr+8:istr+8+3);
        file_month = nc_inst_file(istr+8+5:istr+8+5+1);
    end
% ------------------------------------------------------------------
    
    %(just saving the .mat with the filenames here)
    %save am3_filepaths_chunk, etc. if we are now processing a new year -
    %before am3_multi_load_time_chunks. N.B. on this is done again after
    %the loop, so still gets saved even if we just do one year!!
    %Won't be executed first time around, so doesn't matter if the variables
    %don't exist then
    if strcmp(file_year,file_year_old)==0 & imulti_run>1
        save(gcm_saved_filename_file,'am3_filepaths_chunk','gcm_year_cases','gcm_month_cases','gcm_model_cases');
        % - need to enter the gcm_saved_filename_file filepath in load_gcm_processed_data

%        gcm_saved_filename_file
        load_file{isave_file}=gcm_saved_filename_file
        load_file_year{isave_file}=str2num(file_year_old); %are saving for the previous year
        isave_file = isave_file + 1;
        

        if exist(save_file_list_filename)==2
            save(save_file_list_filename,'load_file','load_file_year','-APPEND');
        else
            save(save_file_list_filename,'load_file','load_file_year');
        end
        
        clear am3_filepaths_chunk gcm_year_cases gcm_model_cases gcm_month_cases
        itam32=1; %reset this counter for gcm_year_cases etc
    end
    file_year_old=file_year;




    % ----------------- run am3_multi_load_time_chunks ---------------
    %this loads, processes and saves the data. Also adds to am3_filepaths_chunk
    %uses the nc_inst open netCDF file
    am3_multi_load_time_chunks
    % ----------------------------------------------------------------

    %re-zip if unzipped the file (saves filespace, but is slow and perhaps risky RE file loss)
    if unzipped==1
%        eval(['!gzip ' nc_dir nc_inst_file]);
        eval(['!rm -f ' temp_gz_inst]);        
    end
    if unzipped_grid==1
%        eval(['!gzip ' nc_dir_grid nc_grid_file]);
        eval(['!rm -f ' temp_gz_grid]);
    end
    if unzipped_prepostLWP==1
%        eval(['!gzip ' nc_dir_grid nc_grid_file]);
        eval(['!rm -f ' temp_gz_prepostLWP]);
    end
    
end

%do a final save of gcm_year_cases, etc and am3_filepaths_chunk etc.
save(gcm_saved_filename_file,'am3_filepaths_chunk','gcm_year_cases','gcm_month_cases','gcm_model_cases');
% - need to enter the gcm_saved_filename_file filepath in load_gcm_processed_data

load_file{isave_file}=gcm_saved_filename_file;
load_file_year{isave_file}=str2num(file_year); %saving for the current year at the end
isave_file = isave_file + 1;


if exist(save_file_list_filename)==3
    save(save_file_list_filename,'load_file','load_file_year','-APPEND');
else
    save(save_file_list_filename,'load_file','load_file_year');
end

fprintf(1,['\nMaster file is :- ' save_file_list_filename '\n']);



%Display the name on the screen - will want to save these somewhere to
%automate the process
%gcm_saved_filename_file


clear ioverride_multi_chunk
catch error_multi
    clear ioverride_multi_chunk
    rethrow(error_multi)
end


