% ACSIS_Robson_paper_plot_timeseries_RUN_multi_generic.m
% This script can be run several times from here :- ACSIS_Robson_paper_load_data_plot_MULTI.m

% savedir_date=['/home/disk/eos1/d.grosvenor/modis_work/ACSIS_Nd_trends/plots_' datestr(now,30) '/'];
% eval(['!mkdir ' savedir_date]); 

%Need to run 
%ACSIS_Robson_paper_load_data  - sets var_ukesm in there
%before this to load in the data. Also sets savedir_date. And sorts out the
%obs data (gcm_area_ etc.)


% --- If running for the 1870-1970 and 1971-2014 trends then best to set
% obs_str='' since otherwise the trend periods get changed automatically
% below
% obs_str is changed automatically in
% ACSIS_Robson_paper_plot_timeseries_generic, so set obs_str each time

%obs_str=''; %If want to manually set the trend years
try

    if ~exist('ioverride_vals')
        ioverride_vals=0;
    end
    
    ino_tit = 1; %whether to not plot a lower title
    
    if ioverride_vals==0        
        iload_DAMIP = 0; %only actually does this if iadd_DAMIP==1 too
        iadd_DAMIP = 0;
        
        iload_AerChemMIP=0; %only actually does this if iadd_AerChemMIP=1 too
        iadd_AerChemMIP=0;
        
        i_CMIP6_multimodel=0; %Runs ACSIS_Robson_paper_load_data_generic_multi_model
        iload_multi_model=0;
        
        iadd_HADGEM=1; %
        iload_HADGEM=1; %Currently done in same was as CMIP6 - may have to change if plot rsut for CMIP6.
        
        iadd_amip=0;
        iload_amip=0;
    end
    
    iplot_lin = 1; %whether to plot straight lines for the linear trends - NOTE - have to set this and iplot_trend, the latter
    %sets which trend line to plot.

season='Annual';
%season='DJF';
%season='MAM';
%season='JJA';
%season='SON';

%% Choose region
% Box for regional averages
%box_region_DRIVER = '1'; %'US outflow'  US small box

% Uses this script to select the region :-
% ACSIS_Robson_paper_choose_regional_box2
regions={'1','3','8','4'};
regions={'4'};
%regions={'20'}; %Full NA region used in 1st ACSIS paper
%regions={'10'};
%regions={'11'}; %latest US outflow region
%regions={'12'}; %US mainland and east coast for emissions.
%regions={'13'}; %NA region up to 50S instead of 60S to avoid sea-ice region.
%regions={'0'}; %Global
%regions={'00'}; %Global -60 to 60
%regions={'01'}; %Global -55 to 60
%regions={'02'}; %Global -50 to 55
%regions={'03'}; %Southern Ocean -60 to -40
%regions={'14'};
%regions={'17'};
%regions={'21'}; %Same lons as region 4, but further south to avoid sea-ice region
%regions={'22'}; %Sea-ice region around Newfoundland
%regions={'23'}; %US high clear-sky bias region
%regions={'24'}; %US Jon Robson AMOC work region.

%regions={'0','00','4','24'};

%box_region_DRIVER = '3'; %'Europe outflow' SW of Spain
%box_region_DRIVER = '4'; %'N. Atlantic basin' All NA region
%box_region_DRIVER = '5'; %UK region
%box_region_DRIVER = '6'; %West of UK region
%box_region_DRIVER = '7'; %20 deg further west of there
%box_region_DRIVER = '8'; %'Mid-Atlantic' middle of Atlantic at US box latitude
%box_region_DRIVER = '9'; %20 deg further west of there

iplot_individual_ens=0; %Set this if plotting the individual ensemble memebers
    % *** USE ACSIS_Robson_paper_plot_timeseries_RUN_multi_ENS_DeepC.m ***

%Can choose different periods for the trends
% yr_start_trend_box = [1851 1970 2001];
% yr_end_trend_box = [1970 2014 2014];
% yr_start_trend_box = [1985]; yr_end_trend_box = [2014]; %DEEP-C comparison
% yr_start_trend_box = [1983]; yr_end_trend_box = [2009]; %PATMOS / ISCCP comparison
% %yr_start_trend_box = [2003]; %For MODIS comparison
% 
% %yr_start_trend_box = [1980]; yr_end_trend_box = [2005];
% 
% iplot_trend=[1 1 0]; %whether to plot the straight trend lines

itrend_box_whisker=2; %default - may be overridden below

% iplot_trend=[1 1]; %whether to plot the straight trend lines
% %iplot_trend_amip=[0 1]; %whether to plot the straight trend lines

if ioverride_vals==0
    
    switch var_ukesm
        case {'clt'}
            obs_str = 'ISCCP';
            obs_str = 'PATMOSx';
            obs_str = 'ISCCP PATMOSx'; %both obs sets
            %obs_str = 'none';
            %obs_str='';
            
        otherwise
            %Use what was set in ACSIS_Robson_paper_load_data.m
    end
    
end

%set in ACSIS_Robson_paper_load_data.m, but not properly
switch obs_str
    case 'SW calc'
        %set in SW script
%         iadd_trend_str=1; %whether to add the trend value to the legend
%         iadd_trend_str=0;
%         itrend_lab=iconstant_trend;
        

    otherwise
        itrend_lab=1; %which trend to put in the legend label
        yr_start_plot = 1900; %where to start the plot from - where have data from more like - plot_year_lims is the xlim range
        yr_start_plot = 1950; %where to start the plot from
        yr_start_plot = 1850; %where to start the plot from        
        %yr_start_plot = 1979; %where to start the plot from
        %yr_start_plot = 1984; %where to start the plot from
        %yr_start_plot = 2000; %where to start the plot from
        %yr_start_plot = 1960;
        
        yr_end_plot = 2018;
        
        iadd_trend_str=0; %whether to add the trend value to the legend
        plot_year_lims = [yr_start_plot-2 yr_end_plot];
        plot_year_lims = [yr_start_plot-7 yr_end_plot]; %Matt C plot
        

end

switch obs_str
    case 'AMIP'
        %yr_start_trend_box = [1870 2003]; yr_end_trend_box = [1970 2014]; %MODIS comparison
        yr_start_trend_box = [1850 1971 1985]; yr_end_trend_box = [1970 2014 2014]; %MODIS comparison
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [2003]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one        
        itrend_box_whisker=3; 
        
    case 'MODIS'
        %yr_start_trend_box = [1870 2003]; yr_end_trend_box = [1970 2014]; %MODIS comparison
        yr_start_trend_box = [1850 1971 2003]; yr_end_trend_box = [1970 2014 2014]; %MODIS comparison
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [2003]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one        
        itrend_box_whisker=3; 
        
    case {'DEEP-C','Deep-C'}     
        yr_start_trend_box = [1850 1971 1985]; yr_end_trend_box = [1970 2014 2014]; %DEEP-C comparison
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [1985]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one
        itrend_amip = 3; %Which trend(s) to plot on the timeseries
        %for UKESM above
        itrend_box_whisker=3; %trend to use to compare to the obs
                
    case {'PATMOSx','ISCCP','ISCCP PATMOSx'}
        yr_start_trend_box = [1850 1971 1983]; yr_end_trend_box = [1970 2014 2009]; %PATMOS / ISCCP comparison
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [1983]; yr_end_trend_box_AMIP = [2009]; %PATMOS / ISCCP comparison
        itrend_box_whisker=3; %trend to use to compare to the obs
        
    case 'SW calc'        
        %yr_start_trend_box = [1870 1985]; yr_end_trend_box = [1970 2014]; %DEEP-C comparison
        %iplot_trend=[1 1]; %whether to plot the straight trend lines
        %Now set in ACSIS_Robson_paper_offline_SW_calcs_Sep2020.m
        
    case 'MAC'
        yr_start_trend_box = [1850 1971 1988]; yr_end_trend_box = [1970 2014 2014]; 
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines       
        yr_start_trend_box_AMIP = [1988]; yr_end_trend_box_AMIP = [2014];
        itrend_box_whisker=3; %trend to use to compare to the obs
        
    case 'MODIS_AOD'
        yr_start_trend_box = [1850 1971 2003]; yr_end_trend_box = [1970 2014 2014]; %MODIS comparison
        iplot_trend=[0 0 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [2003]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one        
        itrend_box_whisker=3; 
        
    case '1850-1970 1971-2014 plot'
        yr_start_trend_box = [1850 1971]; yr_end_trend_box = [1970 2014]; %
        iplot_trend=[1 1]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [1985]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one
        
    case '1850-1970 1971-2014 no plot'
        yr_start_trend_box = [1850 1971]; yr_end_trend_box = [1970 2014]; %
        iplot_trend=[0 0]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [1985]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one
        
    %Might want to also add some of these cases to the switch statement
    %below (if aren't actually using obs).
    
    otherwise
        yr_start_trend_box = [1850 1971 1985]; yr_end_trend_box = [1970 2014 2014]; %
        iplot_trend=[1 1 0]; %whether to plot the straight trend lines
        yr_start_trend_box_AMIP = [1985]; yr_end_trend_box_AMIP = [2014]; %and AMIP - need to set both this and one
        %for UKESM above  
        itrend_box_whisker=3; 
end



%yr_start_trend_box_nudged = yr_start_trend_box_AMIP;
%yr_end_trend_box_nudged = yr_end_trend_box_AMIP;
yr_start_trend_box_nudged = 1982;
yr_end_trend_box_nudged = 2014; %run hasn't finished yet for u-by844


switch season
    case 'DJF'
        yr_start_trend_box = yr_start_trend_box + 1; %start a year later for now - although could do DJF of Y2000 for CERES
end




%yr_title = 1965; %Year to plot the title at.
yr_title = 1980; %Year to plot the title at.
yr_title = yr_start_plot + (2014 - yr_start_plot) *0.6;

yr_PI_bar = 1847.5;
yr_PI_bar = 2016;





%Whether to include the main set of obs - note that 02 is actually the main
%obs and 01 the secondary - should fix this!
inc_obs_02=1;
switch var_ukesm
    case {'SO2_low_anthropogenic_emissions','dust_od550','clwvi','lwpic'}
        iplot_obs02=0;  
        inc_obs_02=0;
    otherwise
        iplot_obs02=1;
end

%For obs where we request mulitple variables here we select which ones to
%include.
switch obs_str                
    case {'ISCCP PATMOSx'}
        obs_str1 = 'PATMOSx';
        obs_str2 = 'ISCCP';   
        inc_obs_01 = 1;
        inc_obs_02 = 1;
    otherwise
        obs_str2 = obs_str; 
        inc_obs_01 = 0; 
end

switch obs_str
    case {'none','','1850-1970 1971-2014 plot','1850-1970 1971-2014 no plot'}
        iplot_obs01=0;  
        inc_obs_01=0;
        iplot_obs02=0;  
        inc_obs_02=0;
end


% switch var_ukesm
%     case {'ts','SO2_low_anthropogenic_emissions','dust_od550','SWTOA Calc'}
%         iplot_PI=0; %Looks like this is not used anymore - use no_PI (set in ACSIS_Robson_paper_load_data.m instead).        
%     otherwise
%         iplot_PI=1;
% end  
 
%Whether to remove the xlabel and xticklabels for the nolab plot -
%generally for the bottom plot in a subfigure
switch var_ukesm
    %case {'lwpic'}
        %iremove_xlabs = 0;        
    case {'SWTOA Calc'}
        %set in offline script
    otherwise
        iremove_xlabs = 1;       
end 

switch var_ukesm
    case {'SWTOA Calc'}
        %set in offline script
    otherwise        
        ipad_legend=0;
        %Pad the first entry to avoid box cutting into the text when have
    %superscripts
        iplot_box_whisker=1;
        period_str='';
end 



clear stats_regions trend_dat_box
for ibox=1:length(regions)
    
    box_region_DRIVER = regions{ibox};
    %ACSIS_Robson_paper_plot_timeseries_deepc
    ACSIS_Robson_paper_plot_timeseries_generic
    
end

%iscreen_land=0;

%%

%clear the override and also clear if there is an error in the code
clear ioverride_vals

catch timser_error
    clear ioverride_vals
    rethrow(timser_error);
end
