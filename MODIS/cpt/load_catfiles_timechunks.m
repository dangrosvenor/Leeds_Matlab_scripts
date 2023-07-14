
if ~exist('ioverride_load_gcm') | ioverride_load_gcm==0
%am3_dataset2 = 'old'; %AM3

%am3_dataset2 = '0.5deg'; %AM3
%am3_dataset2 = 'CAMCLUBB';
%am3_dataset2 = 'AM3CLUBB'; % CLUBBv1, 50km - don't have COSP for this
%am3_dataset2 = 'CAM5_2deg';
%am3_dataset2 = 'CAM5_1deg';

%am3_dataset2 = 'CAMCLUBB_COSP_postDiurnal'; %(CLUBBv2) Now have newer
                                             %files with the pre and post mphys LWP
%am3_dataset2 = 'AM3CLUBB_postDiurnal'; %(CLUBBv2, 50km resolution)
%am3_dataset2 = 'AM3_CLUBBv2_COSP_100km'; %CLUBB w/ COSP
%am3_dataset2 = 'AM3_CLUBBv1_1deg'; %CLUBBv2?? Not v1 ! at 100km, with no COSP
%am3_dataset2 = 'AM3_100km_20130110'; %AM3-base runs at 100km, with COSP

am3_dataset2 = '2deg'; %AM3
%am3_dataset2 = 'AM3_CLUBBv1_2deg'; %AM3CLUBBv1 at 100km, with no COSP
%am3_dataset2 = 'AM3_CLUBBv2_200km'; %CLUBBv2, 2 deg, w/ COSP

%am3_dataset2 = 'CAM5_prepostLWP';
%am3_dataset2 = 'CAM_CLUBB_COSP';  %(CLUBBv1)
%am3_dataset2 = 'CAMCLUBBv2_prepostLWP';


end





clear load_file loadchunks_years
iload_file=1;

am3_dataset = am3_dataset2;
%have decided to make am3_dataset be a long descriptive name for choosing
%between model runs. But in some cases it can be removed when saving to
%avoid variable names that are too long.

switch am3_dataset2
    case 'AM3CLUBB' %AM3 CLUBBv1, 50km - 2007 only, no COSP
        cosp_flag=0;
        
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_AM3_CLUBB_20120813T174243.mat'];
%        iload_file=iload_file+1;
        
%with LTS
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_AM3_CLUBB_20130220T143345.mat'];
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_AM3_CLUBB_20130403T110023.mat'];        
        iload_file=iload_file+1;
        
               
        
        gcm_strs = 'AM3'; %N.B. this is just either CAM5 or AM3 I think (don't need to specify CLUBB).
        

        
        
    
    case 'CAM5_2deg'
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y000101_timechunks_CAM5_old_2deg_20120705T152636.mat'];
        iload_file=iload_file+1;
        
        gcm_strs = 'CAM5';
    
    case 'CAM5_1deg'
        cosp_flag=0;
        %1 degree h1 files (no COSP)
        %2006
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_new_1deg_20120704T150018.mat'];
%        iload_file=iload_file+1;
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_new_1deg_20120705T091336.mat'];
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_new_1deg_20120830T155031.mat'];        
        
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20120914T101616.mat'];                
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20120924T135329.mat'];                


       %1 degree h1 & h2 files (h2 are COSP)
        %2008
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20120914T101616.mat'];                      
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20120925T131212.mat'];                              
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20121130T184711.mat'];          
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20121130T215512.mat'];                               
 %with LTS        
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20130220T175506.mat'];                               
%with LTS & qv       
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20130221T182814.mat'];
%with LTS & qv and SW, LW etc.   
%       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20130530T142618.mat'];
       
       %processed all years
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_new_1deg_COSP_20130724T232906.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_new_1deg_COSP_20130725T011354.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_new_1deg_COSP_20130725T024617.mat']; 
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200911_timechunks_CAM5_new_1deg_COSP_20130725T052918.mat']; 
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y201012_timechunks_CAM5_new_1deg_COSP_20130725T070242.mat']; 
       iload_file=iload_file+1;

       

        
        gcm_strs = 'CAM5';
        
    case 'CAM_CLUBB_COSP' %aka CLUBBv1
        cosp_flag=1;
        
       %1 degree h1 & h2 files (h2 are COSP)
        %2008
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20120926T105745.mat'];                              
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20121201T091301.mat'];                                      

        %with LTS and qv
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20130221T183629.mat'];
        %also with SW fluxes
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20130530T121052.mat'];

        %2006 - 2010 CLUBBv1 files
%          load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_CLUBB_new_1deg_20130725T140358.mat'];
%          iload_file = iload_file + 1;  
%          load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_CLUBB_new_1deg_20130725T160434.mat'];
%          iload_file = iload_file + 1;  
%          load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20130725T171854.mat'];
%          iload_file = iload_file + 1;  
%          load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200911_timechunks_CAM5_CLUBB_new_1deg_20130725T191130.mat'];
%          iload_file = iload_file + 1;  
%          load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y201012_timechunks_CAM5_CLUBB_new_1deg_20130725T205101.mat'];

loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_CAM5_CLUBB_new_1deg_20130828T162035.mat';
loadchunks_years = [2006:2010]; %2006-2010 available - enter range as [2006:2010]
         


         gcm_strs = 'CAM5_CLUBB';
         am3_dataset=''; %to avoid the name being too long
         
    case 'CAMCLUBB_COSP_postDiurnal' %aka CLUBBv2
        cosp_flag=1;
        
         %1 degree h1 & h2 files (h2 are COSP)
        %2008
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAMCLUBB_COSP_postDiurnal_20130123T195159.mat'];                                              
%  Jan 30th - files with the RWP added
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAMCLUBB_COSP_postDiurnal_20130129T191410.mat'];                                                      
%with LTS added
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAMCLUBB_COSP_postDiurnal_20130220T120244.mat'];                                                      
%with qv at 700mb added
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAMCLUBB_COSP_postDiurnal_20130221T144418.mat'];                                                              
%with LTS, qv and SW LW fluxes
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAMCLUBB_COSP_postDiurnal_20130529T183813.mat'];
        
%with LTS, qv and SW LW fluxes
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_CLUBB_new_1deg_20130725T140358.mat'];
        iload_file = iload_file + 1;  
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_CLUBB_new_1deg_20130725T160434.mat'];
        iload_file = iload_file + 1;  
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200812_timechunks_CAM5_CLUBB_new_1deg_20130725T171854.mat'];
        iload_file = iload_file + 1;  
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200911_timechunks_CAM5_CLUBB_new_1deg_20130725T191130.mat'];
        iload_file = iload_file + 1;  
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y201012_timechunks_CAM5_CLUBB_new_1deg_20130725T205101.mat'];
        iload_file = iload_file + 1;  
        





         gcm_strs = 'CAM5_CLUBB';
         am3_dataset=''; %to avoid the name being too long
         
         
    case 'CAM5_prepostLWP' %1deg with pre and post LWC and other new variables
        cosp_flag=1;
        ilwcAPBP=1;
        
%        loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_filenames_dir/save_file_list_timechunks_CAM5_prepostLWP_20130801T115322.mat';

        %select which years to load
%        loadchunks_years = [2006:2009]; %2006-2010 available
        
%        iload_file = 1;
%load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_prepostLWP_20130731T192018.mat'];
%iload_file = iload_file + 1;


%loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_CAM5_prepostLWP_20130801T204826.mat';
%loadchunks_years = [2006:2010]; %2006-2010 available - enter range as [2006:2010]


%Newly processed file with new prepost LWP variables from Pete -
%also have added Nd in the ISCCP low, mid and high height ranges :-
loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_CAM5_prepostLWP_20140210T070452.mat';
loadchunks_years = [2006:2010]; %2006-2010 available - enter range as [2006:2010]

    case  'CAMCLUBBv2_prepostLWP' %1deg with pre and post LWP and other new variables
        cosp_flag=1;
        ilwcAPBP=1;

%        loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_CAMCLUBBv2_prepostLWP_20131024T005341.mat';
        %new version with changed diags for pre and post mphys LWP (first outputs from Pete were wrong) :-
        loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_CAMCLUBBv2_prepostLWP_20140210T022739.mat';
        loadchunks_years = [2006:2010]; %2006-2010 available - enter range as [2006:2010]
                        
    case 'AM3CLUBB_postDiurnal'
        cosp_flag=0;
        %2008 - still need to process the others
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3_CLUBB_postDiurnal_20130128T161942.mat'];                                              
  %with LTS
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3_CLUBB_postDiurnal_20130220T141732.mat'];                                                              
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3_CLUBB_postDiurnal_20130322T155336.mat'];
         gcm_strs = 'AM3_CLUBBv2';
         am3_dataset=''; %to avoid the name being too long
        
    case 'AM3_CLUBBv2_200km'
        
        cosp_flag=1;

        loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_AM3_CLUBBv2_2deg_20131010T143137.mat';
        loadchunks_years = [2006:2010]; %2006-2010 available - enter range as [2006:2010]

        gcm_strs = 'AM3_CLUBBv2';
        am3_dataset=''; %to avoid the name being too long
        
case 'AM3_CLUBBv2_50km'
        
        cosp_flag=1;

        loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_AM3_CLUBBv2_50km_20131013T014756.mat';
        loadchunks_years = [2007]; %2006-2010 available - enter range as [2006:2010]

        gcm_strs = 'AM3_CLUBBv2';
        am3_dataset=''; %to avoid the name being too long        
         
         
    case 'CAMCLUBB'
        cosp_flag=1;
        %1 degree h1 files, CLUBB (w/ COSP)
        %2006
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_CLUBB_new_1deg_20120704T151120.mat'];
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200612_timechunks_CAM5_CLUBB_new_1deg_20120705T084900.mat'];
%        iload_file=iload_file+1;
%        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_CLUBB_new_1deg_20120705T125511.mat'];
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y200711_timechunks_CAM5_CLUBB_new_1deg_20120830T134150.mat'];
        
        iload_file=iload_file+1;
        
        gcm_strs = 'CAM5_CLUBB';
        
    case '0.5deg'
        cosp_flag=1;
        
        %      load(['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_20120622T112004.mat']); %path to a mat file with these stored for a particular half-year
        %0.5 degree (no CFADs)
        %load(['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_20120622T141135.mat']); %path to a mat file with these stored for a particular half-year

%0.5 degree (w/ CFADs, no heights)
%2007_01
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_new_0.5deg_20120627T223738.mat'];
% iload_file=iload_file+1;
% %2007_07
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_new_0.5deg_20120627T075031.mat'];
% iload_file=iload_file+1;
% %2008_01
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080101_timechunks_new_0.5deg_20120627T190758.mat'];
% iload_file=iload_file+1;
% %2008_07
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_new_0.5deg_20120628T132931.mat'];
% iload_file=iload_file+1;
% 
% %2008_07
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_0.5deg_20130225T110023.mat'];
% iload_file=iload_file+1;

%2008_01 and 07  - how does this 01 and 07 work again??!. Think that now I
%just process the whole year together (that is what I did below). Will keep
%each year separate to save memory when loading, though.
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_0.5deg_20130225T142541.mat'];
% iload_file=iload_file+1;

%post diurnal with LTS (inc. LTS1000) & qv700 included
%  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_0.5deg_20130405T151118.mat'];  %with LTS included
%  iload_file=iload_file+1;

%with COSP fields (Nov 2013)
loadchunks_filename = '/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/save_file_list_timechunks_AM3new_0.5deg_20131112T030110.mat';
loadchunks_years = [2007:2008]; %2007-2008 available - enter range as e.g. [2006:2010]
        





gcm_str_select='AM3';

    case '2deg'
        cosp_flag=1;
%2 degree (no heights)
% 2007_01
%load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_new_2deg_20120626T161638.mat'];
%iload_file=iload_file+1;
% 2007_07
%load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_new_2deg_20120627T162226.mat'];
%iload_file=iload_file+1;

%ignore 2008_01 now as is combined with 2008_07 below
%%% %2008_01
%%% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080101_timechunks_new_2deg_20120628T104316.mat'];
%%% iload_file=iload_file+1;

% %2008_01 & 2008_07 - with 0.05 g/m3 LWC thresh
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20120913T101150.mat'];  %old - y20080701_timechunks_new_2deg_20120628T100834.mat
% iload_file=iload_file+1;  

% %2008_01 & 2008_07 - WITHOUT any LWC thresh
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20120913T142436.mat'];  %old - y20080701_timechunks_new_2deg_20120628T100834.mat
% iload_file=iload_file+1;


 
% %2009_01
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090101_timechunks_new_2deg_20120628T105053.mat'];
% iload_file=iload_file+1;
% %2009_07
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090701_timechunks_new_2deg_20120628T110928.mat'];
% iload_file=iload_file+1;
% %2010_01
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100101_timechunks_new_2deg_20120628T113302.mat'];
% iload_file=iload_file+1;
% %2010_07
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100701_timechunks_new_2deg_20120628T124647.mat'];
% iload_file=iload_file+1;


% ?2007_07?
%load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_new_2deg_20120627T142818.mat'];
%iload_file=iload_file+1;


%post-diurnal
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20120926T150301.mat'];  %old - y20080701_timechunks_new_2deg_20120628T100834.mat
% iload_file=iload_file+1; 

% %post-diurnal
%  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20130222T175732.mat'];  
%  iload_file=iload_file+1; 

 %post-diurnal with LTS and LTS1000
 %add other years here
% load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20130410T175356.mat'];  
% iload_file=iload_file+1; 

  %post-diurnal with LTS and LTS1000 - all years 2006-2010
 load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20060701_timechunks_AM3new_2deg_20130724T212440.mat'];  
 iload_file=iload_file+1; 
  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_AM3new_2deg_20130724T213606.mat'];  
 iload_file=iload_file+1; 
  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3new_2deg_20130724T214806.mat'];  
 iload_file=iload_file+1; 
  %load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090701_timechunks_AM3new_2deg_20130724T220009.mat'];  
  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090701_timechunks_AM3new_2deg_20130725T140332.mat'];
 iload_file=iload_file+1; 
%  load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100701_timechunks_AM3new_2deg_20130724T221723.mat'];  
   load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100701_timechunks_AM3new_2deg_20130725T145430.mat'];  
 iload_file=iload_file+1; 
 
gcm_str_select='AM3';

    case 'AM3_CLUBBv2_COSP_100km'
        cosp_flag=1;
 %July 2013 runs when got the pre and post LWP for CAM (but not AM3)
   %The major change for AM3 was that there is now COSP for CLUBB
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20060101_timechunks_AM3_CLUBB_COSP_prepostLWP_20130718T161141.mat'];  
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_AM3_CLUBB_COSP_prepostLWP_20130718T165306.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080101_timechunks_AM3_CLUBB_COSP_prepostLWP_20130718T173618.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090101_timechunks_AM3_CLUBB_COSP_prepostLWP_20130718T181859.mat'];  
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100101_timechunks_AM3_CLUBB_COSP_prepostLWP_20130718T190052.mat'];
       iload_file=iload_file+1;       

        gcm_str_select='AM3';
        
    case 'AM3_CLUBBv1.9_1deg'
        cosp_flag=0;
        %CLUBB runs at 1 degree - no COSP - these are later than v1 and so
        %are prob more like v2
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20060101_timechunks_AM3_CLUBBv1_1deg_20130726T004708.mat'];
        iload_file=iload_file+1;
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_AM3_CLUBBv1_1deg_20130726T005926.mat'];
        iload_file=iload_file+1; 
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080101_timechunks_AM3_CLUBBv1_1deg_20130726T011116.mat'];
        iload_file=iload_file+1; 
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090101_timechunks_AM3_CLUBBv1_1deg_20130726T012246.mat'];
        iload_file=iload_file+1; 
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100101_timechunks_AM3_CLUBBv1_1deg_20130726T013430.mat'];
        iload_file=iload_file+1;

        gcm_str_select='AM3';
        
case 'AM3_CLUBBv1_2deg'
        cosp_flag=0;
        %CLUBBv1 at 2 degree - no COSP        
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20060701_timechunks_AM3_CLUBBv1_2deg_20130726T122524.mat'];
        iload_file=iload_file+1;
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070701_timechunks_AM3_CLUBBv1_2deg_20130726T123053.mat'];
        iload_file=iload_file+1;
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080701_timechunks_AM3_CLUBBv1_2deg_20130726T123620.mat'];
        iload_file=iload_file+1;
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090701_timechunks_AM3_CLUBBv1_2deg_20130726T124200.mat'];
        iload_file=iload_file+1;
        load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100701_timechunks_AM3_CLUBBv1_2deg_20130726T124735.mat'];
        iload_file=iload_file+1;


        gcm_str_select='AM3';
        

    case 'AM3_100km_20130110'
        cosp_flag=1;
       %100 km base runs - presumably ran in Jan 2013 (with COSP) 
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20060101_timechunks_AM3_100km_20130110_20130726T005114.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20070101_timechunks_AM3_100km_20130110_20130726T015355.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20080101_timechunks_AM3_100km_20130110_20130726T024751.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20090101_timechunks_AM3_100km_20130110_20130726T034558.mat'];
       iload_file=iload_file+1;
       load_file{iload_file}=['/home/disk/eos8/d.grosvenor/CPT/saved_processed_gcm_data/y20100101_timechunks_AM3_100km_20130110_20130726T044100.mat'];
       iload_file=iload_file+1;
       
       gcm_str_select='AM3';
end

%% end of selections 


if exist('loadchunks_years')==1 %if are using the new style of loading with a master file
    load(loadchunks_filename);
    load_file_loaded = load_file;
    clear load_file
    for iload_chunk_year = 1:length(loadchunks_years);
       iy = find([load_file_year{:}] ==  loadchunks_years(iload_chunk_year));
       load_file{iload_chunk_year} = load_file_loaded{iy};
    end
end

am3_filepaths_chunk=[];
gcm_year_cases=[];
gcm_month_cases=[];
gcm_model_cases=[];

for iload_file=1:length(load_file)
    [am3_filepaths_chunk2,gcm_year_cases2,gcm_month_cases2,gcm_model_cases2] = load_chunk_file_info(load_file{iload_file});

    %    if iload_file>1


    names=fieldnames(am3_filepaths_chunk2);
    L2=length(names);

    for iam3_chunk=1:L2
        %add the new names to the old structure
        eval_str=['am3_filepaths_chunk.' names{iam3_chunk} '=am3_filepaths_chunk2.' names{iam3_chunk} ';'];
        eval(eval_str);
    end
    %concatenate these
    gcm_year_cases= cat(2,gcm_year_cases,gcm_year_cases2);
    gcm_month_cases= cat(2,gcm_month_cases,gcm_month_cases2);
    gcm_model_cases= cat(2,gcm_model_cases,gcm_model_cases2);
end
    





