%load amsre data
%N.B. AMSR-E is onboard the AQUA satellite

% The files below all contain sst_amsre and lwp_amsre for global data.
% But the matching to MODIS datasets has only been done for certain
% regions. E.g. for the sst_amsre_time3 variables.
% So if want to match for another region can load the global data and run
% amsre_block_average_time

% ------------------------------------------------
load_time3_only=1; % Flag to tell it to only load the time3 variables to save memory
% ------------------------------------------------

amsre_case ='Global 2008 Aqua only';
%amsre_case ='Global 2008 Aqua and Terra';  %Note Terra refers to SST matching - don't use for LWP
  %as just uses the same value for a given day from Aqua
%amsre_case ='VOCALS 2006-2010 Aqua and Terra';
%amsre_case ='VOCALS 2006-2010 Aqua only';  %For comparisons with GCMs for LWP just use Aqua
   %Could consider checking the (daytime)_ results with Terra as has a
   %slightly different overpass time (10:30 vs 13:30)
%amsre_case ='Southern Ocean 2006-2010 Aqua and Terra';

switch amsre_case
    case 'Global 2008 Aqua only'
        %global daily data for 2008 only (but with Dec 2007 and Jan 2009)
        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2007_to_2009_20130430T150542.mat';
        
    case 'Global 2008 Aqua and Terra'
        %global daily data for 2008 only (but with Dec 2007 and Jan 2009)
        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2007_to_2009_Aqua_Terra_global.mat';        

    case 'VOCALS 2006-2010 Aqua and Terra';

        % --------------- AMRE data made coincident with both Aqua and Terra (SSTs
        % - prob not a good idea for LWP)
        % ------VOCALS only region 2006-2010. lat_restrict = [-50 20]; lon_restrict =
        %[-160 -60]; N.B. uses edges in amsre_block_average_time to be consistent
        %with load_saved_modis_vars

        %amsre_loadfile =
        %'/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2005_to_2011_20130812T220913.mat';
        %amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_VOCALS_2006_to_2010_20130814T152756.mat';

        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_2006_to_2010_20130814T105849.mat';

    case 'VOCALS 2006-2010 Aqua only'
        % --------------- Aqua only
        %VOCALS region 2006-2010 for Aqua only (for CALIPSO  and AMSRE comparisons)
        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_VOCALS_AQUAonly_2006_to_2010_20131025T044307.mat';


    case 'Southern Ocean 2006-2010 Aqua and Terra'
        % All longitudes, 40-60S.
        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_SouthernOcean_AquaTerra_2006_to_2010_20140313T073945.mat';
        
    case 'Raw daily global data from multi_read_amsre_daily 2002-2005'
        % To save loading the files in again - result of running multi_read_amsre_daily
        amsre_loadfile = '/home/disk/eos5/d.grosvenor/AMSRE/amsre_daily_global_raw_amsre_2002-2005.mat';
end

if load_time3_only==1
    exclude={'lwp_amsre','sst_amsre'};
    vars_avail = whos('-FILE',amsre_loadfile);
    for ivars=1:length(vars_avail)
        if length(strmatch(vars_avail(ivars).name,exclude{1},'exact'))==0 & length(strmatch(vars_avail(ivars).name,exclude{1},'exact'))==0
            if length(strmatch(vars_avail(ivars).name,'lwp_amsre_VOCALS','exact'))>0 | length(strmatch(vars_avail(ivars).name,'sst_amsre_VOCALS','exact'))>0
                itime3 = strfind(vars_avail(ivars).name,'_VOCALS');
                var_str = vars_avail(ivars).name(1:itime3-1);

                eval([var_str ' = load_mat_to_var(amsre_loadfile,vars_avail(ivars).name);']);

            else
                load(amsre_loadfile,vars_avail(ivars).name);
            end
        end
    end
    gcm_Plat2D_AMSRE = Plat2D_AMSRE_time3;
    gcm_Plon2D_AMSRE = Plon2D_AMSRE_time3;    
    gcm_Plat2D_edges_AMSRE = Plat2D_AMSRE_time3_edges;
    gcm_Plon2D_edges_AMSRE = Plon2D_AMSRE_time3_edges;    
    LAT_AMSRE = minALL(Plat2D_AMSRE_time3) :1: maxALL(Plat2D_AMSRE_time3);
    LON_AMSRE = minALL(Plon2D_AMSRE_time3) :1: maxALL(Plon2D_AMSRE_time3);    
else
    load(amsre_loadfile);
end

if exist('daynum_timeseries3')
    daynum_timeseries3_time3 = daynum_timeseries3;
end

if ~exist('daynum_timeseries3_AMSRE') | ~ exist('gcm_time_UTC_AMSRE') 
    nT = size(lwp_amsre,3);
    daynum_timeseries3_AMSRE = 1:nT;
    gcm_time_UTC_AMSRE=zeros([1 nT]);
end

% whos('-FILE',amsre_loadfile)
%   Name                   Size                     Bytes  Class     Attributes
%
%   day_amsre              1x428                     3424  double
%   lwp_amsre              4-D                  443750400  double
%   lwp_amsre_time3        4-D                  379468800  double
%   month_amsre            1x428                     3424  double
%   sst_amsre            180x360x428            221875200  double
%   sst_amsre_time3      180x360x366            189734400  double
%   year_amsre             1x428                     3424  double
